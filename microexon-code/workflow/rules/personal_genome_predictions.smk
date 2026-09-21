VCF_FILTER = (
  "'FILTER=\"PASS\" &" 
  "  FORMAT/DP>=10 &" 
  "  ((GT=\"het\" & (TYPE=\"snp\" | TYPE=\"mnp\") & GQ>=99 & AD[:1]/(FORMAT/DP)>=0.3 & AD[:1]/(FORMAT/DP)<0.8) | "
  "   (GT=\"het\" & TYPE=\"indel\" & GQ>=90 & AD[:1]/(FORMAT/DP)>=0.3 & AD[:1]/(FORMAT/DP)<0.8) | "
  "   (GT=\"AA\" & GQ>=25 & AD[:1]/(FORMAT/DP)>=0.8))'"
)


rule tabix:
  input:
    vcf="resources/personal_genomes/vcf/{dataset}/{sample}.vcf.gz",
  output:
    "resources/personal_genomes/vcf/{dataset}/{sample}.vcf.gz.tbi",
  conda: "../envs/bcftools.yml"
  shell:
    "tabix -p vcf {input.vcf}"

rule consensus:
  input:
    g="resources/genomes/hg38.fa.gz",
    vcf="resources/personal_genomes/vcf/{dataset}/{sample}.vcf.gz",
    index="resources/personal_genomes/vcf/{dataset}/{sample}.vcf.gz.tbi"
  output:
    twobit=temp("resources/genomes/[{dataset}--{sample}].2bit"),
    chain="results/personal_genomes/consensus/{dataset}/chain/{sample}.chain"
  log:
    "results/personal_genomes/consensus/{dataset}/logs/{sample}.txt",
  conda: "../envs/bcftools.yml"
  shell:
    "bcftools consensus "
    "  {input.vcf} "
    "  -f {input.g} "
    "  -H A "
    "  -o $TMPDIR/{wildcards.sample}.fa " 
    "  -c {output.chain} " 
    f"  -i  {VCF_FILTER} "
    "2> {log} && "
    "faToTwoBit $TMPDIR/{wildcards.sample}.fa {output.twobit} && "
    "rm $TMPDIR/{wildcards.sample}.fa"

rule intersect_variants:
  input:
    vcf="resources/personal_genomes/vcf/{dataset}/{sample}.vcf.gz",
    index="resources/personal_genomes/vcf/{dataset}/{sample}.vcf.gz.tbi",
    events_bed="resources/personal_genomes/wt_bed/MicEvents_known_novel.feature_intervals.hg38.bed"
  output:
    "results/personal_genomes/bed/MicEvents_known_novel.variant_overlaps.[{dataset}--{sample}].bed"
  conda: "../envs/bcftools.yml"
  shell:
    f"bcftools view {{input.vcf}} -i {VCF_FILTER} | "
     "bedtools intersect -a {input.events_bed} -b stdin -loj > {output}"

rule postprocess_intersect_variants:
  input:
    rules.intersect_variants.output[0]
  output:
    "results/personal_genomes/variant_intersection/MicEvents_known_novel.variants_overlaps.[{dataset}--{sample}].tsv"
  conda: "../envs/pybedtools.yaml"
  script: "../scripts/postprocess_intersect_variants.py"


rule liftover:
  input:
    chain=rules.consensus.output.chain,
    up="resources/personal_genomes/wt_bed/MicEvents_known_novel.hg38.upInt.bed",
    dn="resources/personal_genomes/wt_bed/MicEvents_known_novel.hg38.dnInt.bed"
  output:
    up="results/personal_genomes/bed/MicEvents_known_novel.[{dataset}--{sample}].upInt.bed",
    dn="results/personal_genomes/bed/MicEvents_known_novel.[{dataset}--{sample}].dnInt.bed",
    up_un="results/personal_genomes/bed/MicEvents_known_novel.[{dataset}--{sample}].upInt.unmapped.bed",
    dn_un="results/personal_genomes/bed/MicEvents_known_novel.[{dataset}--{sample}].dnInt.unmapped.bed"
  conda : "../envs/bcftools.yml"
  shell:
    """liftOver {input.up} {input.chain} {output.up} {output.up_un}
       liftOver {input.dn} {input.chain} {output.dn} {output.dn_un}"""


rule bed_to_mic_events:
  input:
    wt_events="resources/inputs/personal_genomes/MicEvents_known_novel.hg38.tab.gz",
    up_bed=rules.liftover.output.up,
    dn_bed=rules.liftover.output.dn,
    up_un_bed=rules.liftover.output.up_un,
    dn_un_bed=rules.liftover.output.dn_un,
    wt_genome="resources/genomes/hg38.2bit",
    mt_genome=rules.consensus.output.twobit
  output:
    mt_events="resources/inputs/personal_genomes/MicEvents_known_novel.[{dataset}--{sample}].tab.gz",
    skipped_events="resources/inputs/personal_genomes/MicEvents_known_novel.skipped_events.[{dataset}--{sample}].csv"
  params:
    intron_window=300,
    exon_window=100
  conda: "../envs/base_env.yml"
  script: "../scripts/bed_to_mic_events.py"
