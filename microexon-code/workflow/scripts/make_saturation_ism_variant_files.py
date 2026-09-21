import pandas as pd
import twobitreader
from tqdm import tqdm

def get_variants_for_interval(window, pos_labels):
  variant_events = []
  for pos, pos_label in zip(range(*window), pos_labels):
    pos_1 = pos + 1
    ref = genome[event.chrom][pos].upper()

    for alt in set('ACGT').difference(ref):
      pos_1 = pos + 1
      variant_string = f"{event.chrom}:{pos_1}:{ref}:{alt}"
      event_var = event.copy()
      event_var.name = f"{event_id}_{pos_label}_{ref}:{alt}"
      event_var['variant'] = variant_string
      variant_events.append(event_var)

  return variant_events


if __name__ == "__main__":
  genome = twobitreader.TwoBitFile(snakemake.input[1])
  exons = pd.read_csv(
    snakemake.input[0],
    sep="\t", index_col="event",
    converters={
      'upIntStart': lambda pos: int(pos) - 1,
      'dnIntStart': lambda pos: int(pos) - 1,
    }
  )
  exons = exons.loc[exons.label == snakemake.params['label']]

  window_size = snakemake.params['window']

  for event_id, event in tqdm(exons.iterrows(), total=len(exons)):
    variant_events = []
    c1_pos_labels = [f"c1{i:+d}" for i in range(-window_size, window_size + 1)]

    exon_pos_labels = (
      [f"up{i}" for i in range(-window_size, 0)] +
      [f"ex{i}" for i in range(1, event.lengthDiff + 1)] +
      [f"dn{i:+d}" for i in range(1, window_size + 1)]
    )

    c2_pos_labels = [f"c2{i:+d}" for i in range(-window_size, window_size + 1)]

    if event.strand == '+':
      c1_window = (event.upIntStart - window_size, event.upIntStart + window_size)
      exon_window = (event.upIntEnd - window_size, event.dnIntStart + window_size)
      c2_window = (event.dnIntEnd - window_size, event.dnIntEnd + window_size)
    else:
      c1_window = (event.upIntEnd - window_size, event.upIntEnd + window_size)
      c1_pos_labels = reversed(c1_pos_labels)

      exon_window = (event.dnIntEnd - window_size, event.upIntStart + window_size)
      exon_pos_labels = reversed(exon_pos_labels)

      c2_window = (event.dnIntStart - window_size, event.dnIntStart + window_size)
      c2_pos_labels = reversed(c2_pos_labels)

    variant_events.extend(get_variants_for_interval(c1_window, c1_pos_labels))
    variant_events.extend(get_variants_for_interval(exon_window, exon_pos_labels))
    variant_events.extend(get_variants_for_interval(c2_window, c2_pos_labels))

    variant_events = pd.DataFrame(variant_events)
    variant_events_pos1 = variant_events.copy()
    variant_events_pos1.upIntStart = variant_events_pos1.upIntStart + 1
    variant_events_pos1.dnIntStart = variant_events_pos1.dnIntStart + 1

    output_file = [f for f in snakemake.output if event_id in f][0]
    variant_events_pos1.to_csv(output_file, sep="\t", index_label="event")
