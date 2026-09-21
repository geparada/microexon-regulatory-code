import pandas as pd

known_mics_hg38 = pd.read_table("resources/inputs/wt/MicEvents.hg38.tab.gz")
known_mics_hg38 = known_mics_hg38.loc[known_mics_hg38.label == "hi"]["event"]

rule saturation_ism_make_variant_files_hg38:
  input:
    "resources/inputs/wt/MicEvents.hg38.tab.gz",
    "resources/genomes/hg38.2bit"
  output:
    [f"resources/inputs/saturation_mutagenesis/{event}_variants.hg38.tab.gz" for event in known_mics_hg38]
  params:
    window=300,
    label="hi"
  conda:
    "../envs/base_env.yml"
  script:
    "../scripts/make_saturation_ism_variant_files.py"

known_mics_mm10 = pd.read_table("resources/inputs/wt/MicEvents.mm10.tab.gz")
known_mics_mm10 = known_mics_mm10.loc[known_mics_mm10.label == "hi"]["event"]

rule saturation_ism_make_variant_files_mm10:
  input:
    "resources/inputs/wt/MicEvents.mm10.tab.gz",
    "resources/genomes/mm10.2bit"
  output:
    [f"resources/inputs/saturation_mutagenesis/{event}_variants.mm10.tab.gz" for event in known_mics_mm10]
  params:
    window=300,
    label="hi"
  conda:
    "../envs/base_env.yml"
  script:
    "../scripts/make_saturation_ism_variant_files.py"

putative_detected_hg38, = glob_wildcards(
  "resources/inputs/putative_ISM_input/{event}_variants.hg38.tab.gz"
)
rule predict_saturation_in_silico_mutagenesis_putative_detected:
  input:
    expand(
      "results/predictions/{model}/putative_ISM_input/predictions.{event}_variants.hg38.csv.gz",
      event=putative_detected_hg38,
      model=["known_microexons", "novel_microexons", "original_model"]
    )

control_variants, = glob_wildcards(
  "resources/inputs/Control_variants/{event}.tab.gz"
)
rule predict_control_variants:
  input:
    expand(
      "results/predictions/{model}/Control_variants/predictions.{event}.csv.gz",
      event=control_variants,
      model=["known_microexons", "novel_microexons", "original_model"]
    )

mssng_variants, = glob_wildcards(
  "resources/inputs/MSSNG/{event}.tab.gz"
)
rule predict_mssng_variants:
  input:
    expand(
      "results/predictions/{model}/MSSNG/predictions.{event}.csv.gz",
      event=mssng_variants,
      model=["known_microexons", "novel_microexons", "original_model"]
    )
