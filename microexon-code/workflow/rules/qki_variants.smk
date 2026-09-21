rule generate_qki_variants:
  input:
    "resources/inputs/wt/MicEvents.{genome}.tab.gz",
    "resources/genomes/{genome}.2bit"
  output:
    "resources/inputs/qki/qki_{motif}.{genome}.tab.gz"
  params:
    window_size=150,
    event_label_filter="hi"
  wildcard_constraints:
    motif="[ACGTWSMKRYBDHVN]+"
  conda:
    "../envs/base_env.yml"
  script:
    "../scripts/generate_qki_variants.py"

rule split_qki_input_files:
  input:
    rules.generate_qki_variants.output[0]
  output:
    expand("resources/inputs/qki/qki_split_{split:02d}_{{motif}}.{{genome}}.tab.gz",
      split=range(10))
  conda:
    "../envs/base_env.yml"
  script:
    "../scripts/split_tsv.py"
