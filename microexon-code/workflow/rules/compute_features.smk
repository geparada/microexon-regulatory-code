import pathlib

species_constraint = r"(\w{2,6}\d{1,2})|(\[.+\])"

rule features_seqweaver:
  """Run Seqweaver"""
  input:
    event_file=config["input_file_schema"],
    genome_file=ancient("resources/genomes/{species}.2bit"),
    pytorch_model_file=ancient("resources/seqweaver_model/human.seqweaver.model.cpu.nn.net"),
    feature_name_file=ancient("resources/seqweaver_model/human.feature.colnames"),
  output: "results/features/{task}/features.{prefix}.{species}.seqweaver.csv.gz"
  params:
    use_cuda=config['cuda'],
    batch_size=128,
  wildcard_constraints:
    species=species_constraint
  conda: "../envs/pytorch.yml"
  script: "../scripts/compute_pytorch_features.py"

rule features_length:
  """Computes several length features"""
  input:
    exons=config["input_file_schema"],
    exon_length_stats=ancient("resources/exon_length_stats.hg38.npz")
  output: "results/features/{task}/features.{prefix}.{species}.length.csv.gz"
  wildcard_constraints:
    species=species_constraint
  conda: "../envs/base_env.yml"
  script: "../scripts/compute_length_features.py"


rule features_ss_strength:
  """Compute MaxEntScan scores for splice sites"""
  input:
    exons=config["input_file_schema"],
    genome=ancient("resources/genomes/{species}.2bit")
  output: "results/features/{task}/features.{prefix}.{species}.ss_strength.csv.gz"
  conda: "../envs/base_env.yml"
  wildcard_constraints:
    species=species_constraint
  script: "../scripts/compute_ss_strength_features.py"


rule features_gc_content:
  """Compute GC content in the upstream and downstream introns and the exon"""
  input:
    exons=config["input_file_schema"],
    genome=ancient("resources/genomes/{species}.2bit")
  output: "results/features/{task}/features.{prefix}.{species}.gc_content.csv.gz"
  params:
    intron_window=100
  wildcard_constraints:
    species=species_constraint
  conda: "../envs/base_env.yml"
  script: "../scripts/compute_nucleotide_content.py"


rule features_fimo_fasta_inputs:
  """Generate FASTA files as input to the FIMO tool"""
  input:
    config["input_file_schema"],
    ancient("resources/genomes/{species}.2bit")
  output: temporary("results/features/{task}/tmp/features.{prefix}.{species}.FIMO.input.fa")
  wildcard_constraints:
    species=species_constraint
  conda: "../envs/base_env.yml"
  script: "../scripts/generate_fimo_fasta_inputs.py"


rule features_cisbp:
  """Use FIMO to compute matches of the CisBP motifs"""
  input:
    rules.features_fimo_fasta_inputs.output[0],
    ancient("resources/CisBP-RNA-Homo_sapiens.meme")
  output: "results/features/{task}/features.{prefix}.{species}.CisBP.csv.gz"
  params:
    max_stored_scores=10000000,
    thresh=.5,
    script_path=pathlib.Path(workflow.basedir) / "scripts/aggregate_fimo_output.py"
  conda: "../envs/meme.yml"
  wildcard_constraints:
    species=species_constraint
  shell:
    "fimo "
    "--norc "
    "--max-stored-scores {params.max_stored_scores} "
    "--thresh {params.thresh} "
    "--no-qvalue "
    "--text "
    "--skip-matched-sequence "
    "{input[1]} {input[0]} | "
    "python {params.script_path} {output}"


rule features_manual_rbp_motifs:
  """Compute matches of manually defined RBP motifs"""
  input:
    exons=config["input_file_schema"],
    genome=ancient("resources/genomes/{species}.2bit")
  output: "results/features/{task}/features.{prefix}.{species}.manual_rbp_binding.csv.gz"
  wildcard_constraints:
    species=species_constraint
  conda: "../envs/base_env.yml"
  script: "../scripts/compute_manual_rbp_binding_features.py"


rule features_prepare_svm_bp_fasta_inputs:
  """Create a FASTA file with sequences to use as input for SVM-BP algorithm"""
  input:
    config["input_file_schema"],
    ancient("resources/genomes/{species}.2bit")
  output: temp("results/features/{task}/tmp/features.{prefix}.{species}.svm_bp_inputs.fa")
  params:
    seq_len=200
  wildcard_constraints:
    species=species_constraint
  conda: "../envs/base_env.yml"
  script: "../scripts/generate_svm_bp_fasta_inputs.py"


rule features_run_svm_bp:
  """Run the SVM-BP algorithm. This requires the use of the --use-conda flag"""
  input: rules.features_prepare_svm_bp_fasta_inputs.output[0]
  output: temporary("results/features/{task}/tmp/features.{prefix}.{species}.svm_bp.tab.gz")
  params:
    svm_bp_path=pathlib.Path(workflow.basedir).parent / "svm-bpfinder/svm_bpfinder.py",
    svm_light_path=pathlib.Path(workflow.basedir).parent / config["svm_light_path"]
  wildcard_constraints:
    species=species_constraint
  conda:
    # This script needs to run in Python 2.7
    "../envs/py27_env.yml"
  shadow: "shallow"
  shell: "python {params.svm_bp_path} -i {input} -s Hsap --svm-classify-path={params.svm_light_path} | gzip > {output}"


rule features_svm_bp:
  """Process the output of the SVM-BP algorithm"""
  input:
    svm_bp_scores=rules.features_run_svm_bp.output[0]
  output: "results/features/{task}/features.{prefix}.{species}.svm_bp_scores.csv.gz"
  wildcard_constraints:
    species=species_constraint
  conda: "../envs/base_env.yml"
  script: "../scripts/compute_bp_features.py"


rule features_cossmo:
  """Run the COSSMO algorithm and extract hidden features, logits, and PSI values"""
  input:
    event_file=config["input_file_schema"],
    acceptor_model_path=ancient("resources/cossmo_saved_models/acceptor_100nt_lstm_1_fold0"),
    donor_model_path=ancient("resources/cossmo_saved_models/donor_100nt_lstm_1_fold0"),
    genome_2bit_file=ancient("resources/genomes/{species}.2bit")
  output:
    output_file="results/features/{task}/features.{prefix}.{species}.cossmo.csv.gz"
  wildcard_constraints:
    species=species_constraint
  conda: "../envs/tensorflow.yml"
  script: "../scripts/compute_cossmo_features.py"


rule features_nsr100:
  """Computes known nSR100 binding motif features"""
  input:
       config["input_file_schema"],
       ancient("resources/genomes/{species}.2bit")
  output: "results/features/{task}/features.{prefix}.{species}.nsr100.csv.gz"
  params:
    bin_width=3,
    upstream_window=60
  wildcard_constraints:
    species=species_constraint
  conda: "../envs/base_env.yml"
  script: "../scripts/compute_nsr100_features.py"


rule features_kmers:
  input:
    config["input_file_schema"],
    ancient("resources/genomes/{species}.2bit"),
    ancient("resources/enriched_kmers.csv.gz")
  output: "results/features/{task}/features.{prefix}.{species}.kmers.csv.gz"
  params:
    upstream_window=100,
    bin_width=20
  wildcard_constraints:
    species=species_constraint
  conda: "../envs/base_env.yml"
  script: "../scripts/compute_kmer_features.py"
