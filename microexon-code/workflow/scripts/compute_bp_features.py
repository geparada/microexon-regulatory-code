"""Computes binned branch point features
from the output of the SVM-BP algorithm"""

import pandas as pd


def compute_bp_features(svm_bp, output_file):
  """Takes svm-bp scores, and computes the maximum score BP for each sequence"""
  max_svm_bp = svm_bp.loc[svm_bp.groupby(level=0).svm_scr.idxmax()]\
    .reset_index().set_index('seq_id')
  max_svm_bp = max_svm_bp.loc[~max_svm_bp.index.duplicated(keep='first')]
  alt_ss_bp_score = max_svm_bp.loc[
    [i for i in max_svm_bp.index if i.endswith('_alt')]]
  alt_ss_bp_score.index = [i[:i.rindex('_')] for i in alt_ss_bp_score.index]
  dn_ss_bp_score = max_svm_bp.loc[
    [i for i in max_svm_bp.index if i.endswith('_dn')]]
  dn_ss_bp_score.index = [i[:i.rindex('_')] for i in dn_ss_bp_score.index]

  bp_features = pd.DataFrame(dict(
    alt_ss_bp_score=alt_ss_bp_score.svm_scr,
    alt_ss_bp_dist=alt_ss_bp_score.ss_dist,
    dn_ss_bp_score=dn_ss_bp_score.svm_scr,
    dn_ss_bp_dist=alt_ss_bp_score.ss_dist
  ))
  bp_features.to_csv(output_file, index_label="event")


if __name__ == "__main__":
  svm_bp_file = snakemake.input['svm_bp_scores']
  output_file = snakemake.output[0]

  svm_bp = pd.read_csv(
    svm_bp_file,
    sep="\t", index_col=("seq_id", "ss_dist"),
    converters={'seq_id': lambda x: x.split(' ')[0]}
  )
  compute_bp_features(svm_bp, output_file)
