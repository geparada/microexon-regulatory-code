import json
import os
import pickle

import pandas as pd
import numpy as np
import xgboost as xg

def load_data(config, event_subset=None, cisbp_feature_names=None):
  exons = pd.read_csv(
    config["exons"],
    sep="\t",
    index_col="event",
    converters={
      "upIntStart": lambda pos: int(pos) - 1,
      "dnIntStart": lambda pos: int(pos) - 1,
    },
  )

  if event_subset:
    idx_keep = set(np.where(np.isin(exons.index.values, event_subset))[0] + 1)
    idx_keep.add(0)

    def subset_f(idx):
      return not idx in idx_keep

    exons = exons.loc[exons.index.intersection(event_subset)]
  else:

    def subset_f(idx):
      return False

  features = pd.DataFrame(index=exons.index)

  if "length" in config:
    length_features = pd.read_csv(
      config["length"], index_col="event", skiprows=subset_f
    )
    features = features.join(length_features)

  if "ss_strength" in config:
    ss_strength_features = pd.read_csv(
      config["ss_strength"], index_col="event", skiprows=subset_f
    )
    features = features.join(ss_strength_features)

  if "gc_content" in config:
    gc_content_features = pd.read_csv(
      config["gc_content"], index_col="event", skiprows=subset_f
    )
    features = features.join(gc_content_features)

  if "bp_scores" in config:
    bp_features = pd.read_csv(
      config["bp_scores"], index_col="event", skiprows=subset_f
    )
    features = features.join(bp_features)

  if "manual_rbp" in config:
    manual_rbp_features = pd.read_csv(
      config["manual_rbp"], index_col="event", skiprows=subset_f
    )
    features = features.join(manual_rbp_features)

  if "cis_bp" in config:
    cis_bp_features = pd.read_csv(
      config["cis_bp"], index_col="event", skiprows=subset_f
    )

    if cisbp_feature_names:
      for f in cisbp_feature_names:
        if f not in cis_bp_features.columns:
          cis_bp_features[f] = np.nan
      cis_bp_features = cis_bp_features[cisbp_feature_names]

    features = features.join(cis_bp_features)

  if "cossmo" in config:
    cossmo_features = pd.read_csv(
      config["cossmo"], index_col="event", skiprows=subset_f
    )
    features = features.join(cossmo_features)

  if "seqweaver" in config:
    seqweaver_features = pd.read_csv(
      config["seqweaver"], index_col="event", skiprows=subset_f
    )
    features = features.join(seqweaver_features)

  if "nsr100" in config:
    nsr100_features = pd.read_csv(
      config["nsr100"], index_col="event", skiprows=subset_f
    )
    features = features.join(nsr100_features)

  if "kmers" in config:
    kmer_features = pd.read_csv(
      config["kmers"], index_col="event", skiprows=subset_f
    )
    features = features.join(kmer_features)

  features = features.astype('float')
  return exons, features


def load_data_multispecies(config, species, bg=False, add_species_suffix=False):
  exons = {}
  features = {}

  for s in species:
    exons[s] = pd.read_csv(
      config["exons"][s],
      sep="\t",
      index_col="event",
      converters={
        "upIntStart": lambda pos: int(pos) - 1,
        "dnIntStart": lambda pos: int(pos) - 1,
      },
    )
    exons[s]["species"] = s
    features[s] = pd.DataFrame(index=exons[s].index)

    if "length" in config:
      length_features = pd.read_csv(config["length"][s], index_col="event")
      features[s] = features[s].join(length_features)

    if "ss_strength" in config:
      ss_strength_features = pd.read_csv(
        config["ss_strength"][s], index_col="event"
      )
      features[s] = features[s].join(ss_strength_features)

    if "gc_content" in config:
      gc_content_features = pd.read_csv(
        config["gc_content"][s], index_col="event"
      )
      features[s] = features[s].join(gc_content_features)

    if "bp_scores" in config:
      bp_features = pd.read_csv(config["bp_scores"][s], index_col="event")
      features[s] = features[s].join(bp_features)

    if "manual_rbp" in config:
      manual_rbp_features = pd.read_csv(
        config["manual_rbp"][s], index_col="event"
      )
      features[s] = features[s].join(manual_rbp_features)

    if "rna_compete" in config:
      rna_compete_features = pd.read_csv(
        config["rna_compete"][s], index_col="event"
      )
      features[s] = features[s].join(rna_compete_features)

    if "cossmo" in config:
      cossmo_features = pd.read_csv(config["cossmo"][s], index_col="event")
      features[s] = features[s].join(cossmo_features)

    if "seqweaver" in config:
      seqweaver_features = pd.read_csv(config["seqweaver"][s], index_col="event")
      features[s] = features[s].join(seqweaver_features)

    if "nsr100" in config:
      nsr100_features = pd.read_csv(config["nsr100"][s], index_col="event")
      features[s] = features[s].join(nsr100_features)

    if "kmers" in config:
      kmer_features = pd.read_csv(config["kmers"][s], index_col="event")
      features[s] = features[s].join(kmer_features)

    if bg and add_species_suffix:
      exons[s].index = ["{}_{}".format(ex, s) for ex in exons[s].index]
      features[s].index = ["{}_{}".format(ex, s) for ex in features[s].index]

  exons = pd.concat(exons.values())
  features = pd.concat(features.values())

  features = features.astype('float')
  return exons, features


def prepare_data(files, keys=None, set_labels=None):

  if keys:
    files["foreground"] = {
      k: v for k, v in files["foreground"].items() if k in keys + ["exons"]
    }
    files["background"] = {
      k: v for k, v in files["background"].items() if k in keys + ["exons"]
    }

  foreground_exons, foreground_features = load_data(files["foreground"])
  background_exons, background_features = load_data(files["background"])
  background_exons["label"] = "crypt"

  all_exons = pd.concat([foreground_exons, background_exons], sort=False)
  all_exons.label[all_exons.label.isna()] = "na"

  all_features = pd.concat([foreground_features, background_features], sort=False)

  if not set_labels:
    foreground_targets = pd.Series(1, index=foreground_exons.index)
    background_targets = pd.Series(0, index=background_exons.index)
  else:
    if isinstance(set_labels["foreground"], (list, tuple)):
      foreground_labels = set(set_labels["foreground"])
    else:
      foreground_labels = set([set_labels["foreground"]])

    if isinstance(set_labels["background"], (list, tuple)):
      background_labels = set(set_labels["background"])
    else:
      background_labels = set([set_labels["background"]])

    foreground_exons = all_exons[all_exons.label.isin(foreground_labels)]
    background_exons = all_exons[all_exons.label.isin(background_labels)]
    all_exons = pd.concat([foreground_exons, background_exons])
    all_features = all_features.loc[all_exons.index]

    foreground_targets = pd.Series(1, foreground_exons.index)
    background_targets = pd.Series(0, background_exons.index)

  all_targets = pd.concat([foreground_targets, background_targets])
  return all_exons, all_features, all_targets


def prepare_data_multispecies(
    files, species, keys=None, set_labels=None, add_species_suffix=False, multiclass=False
):

  if keys:
    files["foreground"] = {
      k: v for k, v in files["foreground"].items() if k in keys + ["exons"]
    }
    files["background"] = {
      k: v for k, v in files["background"].items() if k in keys + ["exons"]
    }

  foreground_exons, foreground_features = load_data_multispecies(
    files["foreground"], species
  )
  background_exons, background_features = load_data_multispecies(
    files["background"], species, bg=True, add_species_suffix=add_species_suffix
  )
  background_exons["label"] = "crypt"

  all_exons = pd.concat([foreground_exons, background_exons], sort=False)
  all_exons.label[all_exons.label.isna()] = "na"

  all_features = pd.concat([foreground_features, background_features], sort=False)

  if multiclass:
    hi_targets = pd.Series(1, index=foreground_exons.index[foreground_exons.label == 'hi'])
    not_hi_targets = pd.Series(2, index=foreground_exons.index[foreground_exons.label != 'hi'])
    background_targets = pd.Series(0, index=background_exons.index)
    all_targets = pd.concat([hi_targets, not_hi_targets, background_targets])

    all_exons = all_exons.loc[all_targets.index]
    all_features = all_features.loc[all_targets.index]
  elif set_labels:
    if isinstance(set_labels["foreground"], (list, tuple)):
      foreground_labels = set(set_labels["foreground"])
    else:
      foreground_labels = set([set_labels["foreground"]])

    if isinstance(set_labels["background"], (list, tuple)):
      background_labels = set(set_labels["background"])
    else:
      background_labels = set([set_labels["background"]])

    foreground_exons = all_exons[all_exons.label.isin(foreground_labels)]
    background_exons = all_exons[all_exons.label.isin(background_labels)]
    all_exons = pd.concat([foreground_exons, background_exons])
    all_features = all_features.loc[all_exons.index]

    foreground_targets = pd.Series(1, foreground_exons.index)
    background_targets = pd.Series(0, background_exons.index)
    all_targets = pd.concat([foreground_targets, background_targets])
  else:
    foreground_targets = pd.Series(1, index=foreground_exons.index)
    background_targets = pd.Series(0, index=background_exons.index)
    all_targets = pd.concat([foreground_targets, background_targets])

  return all_exons, all_features, all_targets


def get_data_splits(exons, test_set_event_ids, validation_split=0.1):
  train_events = exons[~exons.index.isin(test_set_event_ids)].index
  n_train_events = len(train_events)
  n_validation_events = int(validation_split * n_train_events)
  validation_event_idx = np.random.choice(n_train_events, n_validation_events)
  train_event_idx = list(set(range(n_train_events)) - set(validation_event_idx))
  validation_events = train_events[validation_event_idx]
  train_events = train_events[train_event_idx]

  validation_events = validation_events[
    exons.loc[validation_events].species == "hg38"
    ]
  test_events = exons.index[exons.index.isin(test_set_event_ids)]

  return train_events, validation_events, test_events


def train_xgboost(
    params,
    train_features,
    validation_features,
    test_features,
    train_targets,
    validation_targets,
    test_targets,
    early_stopping_rounds=20,
    max_epochs=200,
    callbacks=None
):
  xg_train = xg.DMatrix(train_features, label=train_targets)
  xg_validation = xg.DMatrix(validation_features, label=validation_targets)
  xg_test = xg.DMatrix(test_features, label=test_targets)
  evallist = [(xg_train, "train"), (xg_validation, "eval")]
  evals_result = {}
  bst = xg.train(
    params,
    xg_train,
    max_epochs,
    evallist,
    early_stopping_rounds=early_stopping_rounds,
    evals_result=evals_result,
    callbacks=callbacks
  )
  train_performance = bst.attributes()["best_msg"]
  test_predictions = bst.predict(xg_test, ntree_limit=bst.best_ntree_limit)
  return bst, train_performance, test_predictions, evals_result


def save_model(bst, evals_result, test_predictions, test_targets, save_path):
  os.makedirs(save_path, exist_ok=True)
  pickle.dump(bst, open(save_path / "xgboost_model.pkl", "wb"))
  json.dump(evals_result, open(save_path / "training_stats.json", "w"))

  if len(test_predictions.shape) > 1:
    test_predictions = pd.DataFrame(
      dict(
        neural_high_pred_hi=test_predictions[:, 0],
        neural_high_pred_crypt=test_predictions[:, 1],
        neural_high_pred_non_hi=test_predictions[:, 2],
        neural_high_true=test_targets
      ),
      index=test_targets.index
    )
  else:
    test_predictions = pd.DataFrame(
      dict(neural_high_pred=test_predictions, neural_high_true=test_targets),
      index=test_targets.index,
    )
  test_predictions.to_csv(save_path / "test_predictions.csv.gz", index_label="event")


def get_species_from_filename(fname, species):
  for s in species:
    if s in str(fname):
      return s


def reformat_input_dict(input_dict, species):
  files = {
    "foreground": {
      "exons": {
        get_species_from_filename(f, species): f
        for f in input_dict["foreground_exons"]
      },
      "length": {
        get_species_from_filename(f, species): f
        for f in input_dict["foreground_length"]
      },
      "ss_strength": {
        get_species_from_filename(f, species): f
        for f in input_dict["foreground_ss_strength"]
      },
      "gc_content": {
        get_species_from_filename(f, species): f
        for f in input_dict["foreground_gc_content"]
      },
      "bp_scores": {
        get_species_from_filename(f, species): f
        for f in input_dict["foreground_bp_scores"]
      },
      "manual_rbp": {
        get_species_from_filename(f, species): f
        for f in input_dict["foreground_manual_rbp"]
      },
      "rna_compete": {
        get_species_from_filename(f, species): f
        for f in input_dict["foreground_rna_compete"]
      },
      "cossmo": {
        get_species_from_filename(f, species): f
        for f in input_dict["foreground_cossmo"]
      },
      "seqweaver": {
        get_species_from_filename(f, species): f
        for f in input_dict["foreground_seqweaver"]
      },
      "nsr100": {
        get_species_from_filename(f, species): f
        for f in input_dict["foreground_nsr100"]
      },
      "kmers": {
        get_species_from_filename(f, species): f
        for f in input_dict["foreground_kmers"]
      },
    },
    "background": {
      "exons": {
        get_species_from_filename(f, species): f
        for f in input_dict["background_exons"]
      },
      "length": {
        get_species_from_filename(f, species): f
        for f in input_dict["background_length"]
      },
      "ss_strength": {
        get_species_from_filename(f, species): f
        for f in input_dict["background_ss_strength"]
      },
      "gc_content": {
        get_species_from_filename(f, species): f
        for f in input_dict["background_gc_content"]
      },
      "bp_scores": {
        get_species_from_filename(f, species): f
        for f in input_dict["background_bp_scores"]
      },
      "manual_rbp": {
        get_species_from_filename(f, species): f
        for f in input_dict["background_manual_rbp"]
      },
      "rna_compete": {
        get_species_from_filename(f, species): f
        for f in input_dict["background_rna_compete"]
      },
      "cossmo": {
        get_species_from_filename(f, species): f
        for f in input_dict["background_cossmo"]
      },
      "seqweaver": {
        get_species_from_filename(f, species): f
        for f in input_dict["background_seqweaver"]
      },
      "nsr100": {
        get_species_from_filename(f, species): f
        for f in input_dict["background_nsr100"]
      },
      "kmers": {
        get_species_from_filename(f, species): f
        for f in input_dict["background_kmers"]
      },
    },
  }
  return files