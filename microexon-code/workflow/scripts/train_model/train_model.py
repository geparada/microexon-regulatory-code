import json
import re
from collections import defaultdict

import pandas as pd
import xgboost


def train_xgboost(
    params,
    train_data,
    validation_data,
    test_data,
    early_stopping_rounds=20,
    max_epochs=200,
    callbacks=None
):
  evallist = [(train_data, "train"), (validation_data, "eval")]
  evals_result = {}
  bst = xgboost.train(
    params,
    train_data,
    max_epochs,
    evallist,
    early_stopping_rounds=early_stopping_rounds,
    evals_result=evals_result,
    callbacks=callbacks
  )
  # train_performance = bst.attributes()["best_msg"]
  test_predictions = bst.predict(test_data, ntree_limit=bst.best_ntree_limit)
  return bst, test_predictions, evals_result


def main(
    training_data,
    test_data,
    validation_data,
    test_events_fg,
    test_events_bg,
    xgboost_model_out,
    training_stats_out,
    test_predictions_out,
    max_epochs=500,
    early_stopping_rounds=20,
    xgboost_params=None,
    scale_pos_weight=False
):
  params = {
    "booster": "gbtree",
    "max_depth": 2,
    "eta": 0.1,
    "objective": "reg:logistic",
    "lambda": 0.1,
    "alpha": 0.8,
    "nthread": 12,
    "eval_metric": ["error", "auc", "aucpr", "logloss"],
  }
  if xgboost_params is not None:
    params.update(xgboost_params)

  training_data = xgboost.DMatrix(training_data)
  test_data = xgboost.DMatrix(test_data)
  validation_data = xgboost.DMatrix(validation_data)

  labels = training_data.get_label()
  if scale_pos_weight:
    scale_pos_weight = (labels == 0).sum()/labels.sum()
    params.update({"scale_pos_weight": scale_pos_weight})

  bst, test_predictions, evals_result = train_xgboost(
    params,
    training_data,
    validation_data,
    test_data,
    max_epochs=max_epochs,
    early_stopping_rounds=early_stopping_rounds
  )

  bst.save_model(xgboost_model_out)
  json.dump(evals_result, open(training_stats_out, "w"))

  test_events_fg = pd.read_csv(
    test_events_fg,
    sep="\t",
    index_col="event",
  )
  test_events_bg = pd.read_csv(
    test_events_bg,
    sep="\t",
    index_col="event",
  )
  test_events = pd.concat([test_events_fg, test_events_bg])

  test_predictions = pd.DataFrame(
    dict(neural_high_pred=test_predictions, neural_high_true=test_data.get_label()),
    index=test_events.index,
  )
  test_predictions.to_csv(test_predictions_out, index_label="event")

  return bst, test_predictions, evals_result


if __name__ == "__main__":
  main(**snakemake.input, **snakemake.output, **snakemake.params)
