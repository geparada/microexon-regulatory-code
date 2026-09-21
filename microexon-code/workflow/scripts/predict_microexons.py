import xgboost as xg
import pandas as pd
from xgboost_utils import load_data
import pathlib
import pickle
import numpy as np


def fill_missing_features(features, bst):
  missing_features = [f for f in bst.feature_names if f not in features.columns]
  for feature in missing_features:
    features[feature] = np.nan
  features = features[bst.feature_names]
  return features


def main(input_dict, output_path, feature_keys=None):

  if "cisbp_features" in input_dict:
    cisbp_features = open(input_dict["cisbp_features"], "rt").read().split("\n")
  else:
    cisbp_features = None

  if feature_keys:
    input_dict = {k: input_dict[k] for k in feature_keys + ["exons", "model"]}

  exons, features = load_data(input_dict, cisbp_feature_names=cisbp_features)

  if input_dict['model'].endswith('pkl'):
    bst = pickle.load(open(input_dict['model'], 'rb'))
  elif input_dict['model'].endswith('json'):
    bst = xg.XGBClassifier()
    bst.load_model(input_dict['model'])
  else:
    raise ValueError(f"Invalid file extension: {input_dict['model']}")
  
  # features = fill_missing_features(features, bst)
  # features_xg = xg.DMatrix(features)
  predictions = bst.predict_proba(features)[:, 1]
  predictions = pd.DataFrame(dict(
    neural_high_pred=predictions),
    index=exons.index
  )
  if 'label' in exons:
    predictions['label'] = exons.label
  predictions.to_csv(output_path, index_label="event")


if __name__ == "__main__":
  main(dict(snakemake.input), snakemake.output[0], **snakemake.params)
