import pandas as pd


def load_events(files, filter_neural_high=True):
  if not isinstance(files, (list, tuple)):
    files = [files]

  events = []

  for file in files:
    df = pd.read_csv(
      file,
      sep="\t",
      index_col="event",
    )

    if "label" in df and filter_neural_high:
      df = df.loc[df.label == 'hi']

    events.append(df)

  events = pd.concat(events)
  return events


def main(known_foreground_events, known_background_events, test_set_event_ids,
         training_set_foreground, training_set_background, test_set_foreground,
         test_set_background, validation_set_foreground, validation_set_background,
         validation_split):
  events_fg = load_events(known_foreground_events)
  events_bg = load_events(known_background_events, filter_neural_high=False)
  test_set_event_ids = pd.read_table(test_set_event_ids, index_col=0, header=None).index

  test_events_fg = events_fg.loc[events_fg.index.isin(test_set_event_ids)]
  test_events_bg = events_bg.loc[events_bg.index.isin(test_set_event_ids)]

  train_events_fg = events_fg.loc[~events_fg.index.isin(test_set_event_ids)]
  train_events_bg = events_bg.loc[~events_bg.index.isin(test_set_event_ids)]

  validation_events_fg = train_events_fg.sample(frac=validation_split)
  validation_events_bg = train_events_bg.sample(frac=validation_split)

  train_events_fg = train_events_fg.loc[~train_events_fg.index.isin(validation_events_fg.index)]
  train_events_bg = train_events_bg.loc[~train_events_bg.index.isin(validation_events_bg.index)]

  train_events_fg.to_csv(training_set_foreground, sep="\t", index_label="event")
  train_events_bg.to_csv(training_set_background, sep="\t", index_label="event")
  test_events_fg.to_csv(test_set_foreground, sep="\t", index_label="event")
  test_events_bg.to_csv(test_set_background, sep="\t", index_label="event")
  validation_events_fg.to_csv(validation_set_foreground, sep="\t", index_label="event")
  validation_events_bg.to_csv(validation_set_background, sep="\t", index_label="event")


if __name__ == "__main__":
  main(
    **snakemake.input,
    **snakemake.output,
    **snakemake.params
  )
