import re
from collections import defaultdict

import pandas as pd
import xgboost

FEATURE_FILE_PATTERN = re.compile(r".+\.(.+)\.csv\.gz")


def load_events(file):
  events = pd.read_csv(
    file,
    sep="\t",
    index_col="event",
  )

  return events


def load_data(file_names, events: pd.Index, feature_names):
  features_all = defaultdict(list)

  for f in file_names:
    m = FEATURE_FILE_PATTERN.match(str(f))
    feature = m.group(1)

    df = pd.read_csv(f, index_col="event")
    df = df.loc[events.intersection(df.index)]

    features_all[feature].append(df)

  features_out = None
  for feature in feature_names:
    f = pd.concat(features_all[feature])
    if features_out is None:
      features_out = f
    else:
      features_out = features_out.join(f)

  return features_out


def main(
    events_fg,
    events_bg,
    features_fg,
    features_bg,
    features_out,
    feature_names
):

  events_fg = load_events(events_fg).index
  events_bg = load_events(events_bg).index
  features_fg = load_data(features_fg, events_fg, feature_names)
  features_bg = load_data(features_bg, events_bg, feature_names)
  labels_fg = pd.Series(1, index=events_fg)
  labels_bg = pd.Series(0, index=events_bg)

  features = pd.concat([features_fg, features_bg])
  labels = pd.concat([labels_fg, labels_bg])

  features_xg = xgboost.DMatrix(features, labels)
  features_xg.save_binary(features_out)


if __name__ == "__main__":
  main(features_out=snakemake.output[0], **snakemake.input, **snakemake.params)