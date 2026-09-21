from math import ceil

import yaml
import pathlib


feature_keys = yaml.load(
  open(pathlib.Path(config["model_path"]) / "feature_set.yaml"),
  Loader=yaml.FullLoader
)['feature_set']

def divup(a, b):
  return int(ceil(float(a) / b))


rule predictions_for_model:
  input:
    model="resources/models/{model}/xgboost_model.json",
    exons=config["input_file_schema"],
    length=rules.features_length.output[0],
    ss_strength=rules.features_ss_strength.output[0],
    gc_content=rules.features_gc_content.output[0],
    bp_scores=rules.features_svm_bp.output[0],
    manual_rbp=rules.features_manual_rbp_motifs.output[0],
    cis_bp=rules.features_cisbp.output[0],
    cossmo=rules.features_cossmo.output[0],
    seqweaver=rules.features_seqweaver.output[0],
    nsr100=rules.features_nsr100.output[0],
    kmers=rules.features_kmers.output[0],
    cisbp_features="resources/cisbp_feature_names.txt"
  output:
    "results/predictions/{model}/{task}/predictions.{prefix}.{species}.csv.gz"
  params:
    feature_keys=feature_keys
  conda:
    "../envs/predictions.yml"
  wildcard_constraints:
    model="known_microexons|novel_microexons|original_model|known_microexons_hg38",
    species=r"(\w{2,6}\d{1,2})|(\[.+\])"
  script:
    "../scripts/predict_microexons.py"
