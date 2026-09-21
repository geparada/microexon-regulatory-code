import twobitreader
import pandas as pd
import gzip
import csv

from variant_genome import Interval, get_variant_genome_from_series

NSR100_MOTIFS = ['TGCTGC', 'CTGCTG', 'GTGCC', 'GTGCT', 'TGGACG']


def compute_nsr100_features_for_event(
    event: pd.Series, wt_genome: twobitreader.TwoBitFile, bin_width: int, upstream_window: int) -> dict:

  mt_genome = get_variant_genome_from_series(event, wt_genome)

  if event.strand == '+':
    upstream_int = Interval(
      event.chrom, event.strand, event.upIntEnd - upstream_window, event.upIntEnd, event.upIntEnd, None)
  else:
    upstream_int = Interval(
      event.chrom, event.strand, event.upIntStart, event.upIntStart + upstream_window, event.upIntStart, None)

  upstream_seq = mt_genome.get_seq(upstream_int).upper()

  feature_profile = [upstream_seq[i:i+6] in NSR100_MOTIFS for i in range(0, len(upstream_seq) - 5)]
  features = dict(event=event.name)
  features.update({
    "nsr100_hexamer_{}_{}".format(i - upstream_window, i - upstream_window + bin_width):
      any(feature_profile[i:i+bin_width])
    for i in range(0, len(feature_profile)-bin_width, bin_width)
  })
  return features


def test_compute_nsr100_motifs():
  event = pd.Series(dict(chrom='chr17', strand='+', upIntEnd=76933291), name='HsaEX0039104')
  genome = twobitreader.TwoBitFile('EXTERNAL_DATA/hg38.2bit')
  features = compute_nsr100_features_for_event(event, genome, 3, 60)
  assert features['nsr100_hexamer_-48_-45']


def compute_nsr100_features(event_file, genome_file, output_file, bin_width, upstream_window):
  events = pd.read_csv(
    event_file,
    sep="\t", index_col="event",
    converters={
      'upIntStart': lambda pos: int(pos) - 1,
      'dnIntStart': lambda pos: int(pos) - 1,
    }
  )

  genome = twobitreader.TwoBitFile(genome_file)

  with gzip.open(output_file, 'wt') as output_f:
    csv_writer = None
    for event_id, event in events.iterrows():
      try:
        features = compute_nsr100_features_for_event(event, genome, bin_width, upstream_window)
      except:
        print(f"Error in '{event_id}.")
        raise

      if csv_writer is None:
        csv_writer = csv.DictWriter(output_f, fieldnames=features.keys())
        csv_writer.writeheader()

      csv_writer.writerow(features)


if __name__ == "__main__":
  compute_nsr100_features(*snakemake.input, *snakemake.output, **snakemake.params)
