"""Computes the following length features for each alternative splicing event:
- length_exon: Length of the alternative exon
- length_upInt: Length of the upstream intron
- length_dnInt: Length of the downstream intron
- loglen_{exon|upInt|dnInt}: Log10 lengths
- zlen_{exon|upInt|dnInt}: z-scores
"""

#%% Imports
import numpy as np
import pandas as pd


#%% Length features
def zlen(x, exon_length_stats):
  return ((x - exon_length_stats['mean']) / exon_length_stats['std'])


def compute_length_features(exons, output_file, exon_length_stats):
  info_lengths = pd.DataFrame(dict(
    length_exon=exons.lengthDiff,
    length_upInt=exons.upIntEnd - exons.upIntStart,
    length_dnInt=exons.dnIntEnd - exons.dnIntStart,
    loglen_exon=np.log10(exons.lengthDiff),
    loglen_upInt=np.log10(exons.upIntEnd - exons.upIntStart),
    loglen_dnInt=np.log10(exons.dnIntEnd - exons.dnIntStart),
    zlen_exon=zlen(exons.lengthDiff, exon_length_stats),
    zlen_upInt=zlen(exons.upIntEnd - exons.upIntStart, exon_length_stats),
    zlen_dnInt=zlen(exons.dnIntEnd - exons.dnIntStart, exon_length_stats)
  ))
  info_lengths.to_csv(output_file, index_label="event")


if __name__ == "__main__":
  try:
    snakemake
  except NameError:
    snakemake = None

  if snakemake:
    exon_file = snakemake.input['exons']
    exon_length_stats_file = snakemake.input['exon_length_stats']
    output_file = snakemake.output[0]
  else:
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument(
      "exome_file", help="Tab-delimited file with exon events from Vast-tools")
    parser.add_argument("exon_length_stats")
    parser.add_argument("output_file_path", help="CSV output file")
    args = parser.parse_args()
    exon_file = args.exon_file
    exon_length_stats_file = args.exon_length_stats
    output_file = args.output_file

  exons = pd.read_csv(
    exon_file,
    sep="\t", index_col="event",
    converters={
      'upIntStart': lambda pos: int(pos) - 1,
      'dnIntStart': lambda pos: int(pos) - 1,
    }
  )
  exon_length_stats = np.load(exon_length_stats_file)
  compute_length_features(exons, output_file, exon_length_stats)
