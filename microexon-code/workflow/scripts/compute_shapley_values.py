import xgboost as xg
import pandas as pd
from xgboost_utils import load_data
import pathlib
import pickle
import numpy as np
import shap


def fill_missing_features(features, bst):
  missing_features = [f for f in bst.feature_names if f not in features.columns]
  for feature in missing_features:
    features[feature] = np.nan
  features = features[bst.feature_names]
  return features


def main(input_dict, output_path, feature_keys=None):

  if feature_keys:
    input_dict = {k: input_dict[k] for k in feature_keys + ["exons", "model"]}
  exons, features = load_data(input_dict)
  bst = pickle.load(open(input_dict['model'], 'rb'))

  explainer = shap.TreeExplainer(bst)

  features = fill_missing_features(features, bst)
  shap_values = explainer.shap_values(features)

  shap_values = pd.DataFrame(shap_values, columns=features.columns, index=features.index)
  shap_values.to_csv(output_path, index_label="event")


if __name__ == "__main__":
  main(dict(snakemake.input), snakemake.output[0], **snakemake.params)
