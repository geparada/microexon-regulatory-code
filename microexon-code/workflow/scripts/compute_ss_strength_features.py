"""Compute splice site strength features using MaxEntScan"""

#%%

# %%
import csv
import gzip

import tqdm
import twobitreader
from maxentpy import maxent_fast as maxent
from maxentpy.maxent import load_matrix5, load_matrix3

from variant_genome import Interval, get_variant_genome_from_dict

matrix3 = load_matrix3()
matrix5 = load_matrix5()


#%% Splice site strength features
def maxent3_score(event, ss_name, genome):
  if event["strand"] == '+':
    if ss_name == "up":
      ss_3p = event["upIntEnd"]
    elif ss_name == "dn":
      ss_3p = event["dnIntEnd"]
    else:
      raise ValueError
    # seq = wt_genome[event["chrom"]][ss_3p - 20:ss_3p + 3].upper()
    seq = genome.get_seq(Interval(event["chrom"], event["strand"], ss_3p - 20, ss_3p + 3, ss_3p, None))
  else:
    if ss_name == "up":
      ss_3p = event["upIntStart"]
    elif ss_name == "dn":
      ss_3p = event["dnIntStart"]
    else:
      raise ValueError
    # seq = reverse_complement(wt_genome[event["chrom"]][ss_3p - 3:ss_3p + 20].upper())
    seq = genome.get_seq(Interval(event["chrom"], event["strand"], ss_3p - 3, ss_3p + 20, ss_3p, None))
  try:
    return maxent.score3(seq, matrix3)
  except KeyError:
    return float('nan')


def maxent5_score(event, ss_name, genome):
  if event["strand"] == '+':
    if ss_name == 'dn':
      ss_5p = event["dnIntStart"]
    elif ss_name == 'up':
      ss_5p = event["upIntStart"]
    else:
      raise ValueError
    # seq = wt_genome[event["chrom"]][ss_5p - 3:ss_5p + 6].upper()
    seq = genome.get_seq(Interval(event["chrom"], event["strand"], ss_5p - 3, ss_5p + 6, ss_5p, None))
  else:
    if ss_name == 'dn':
      ss_5p = event["dnIntEnd"]
    elif ss_name == 'up':
      ss_5p = event["upIntEnd"]
    else:
      raise ValueError
    # seq = reverse_complement(wt_genome[event["chrom"]][ss_5p - 6:ss_5p + 3].upper())
    seq = genome.get_seq(Interval(event["chrom"], event["strand"], ss_5p - 6, ss_5p + 3, ss_5p, None))
  try:
    return maxent.score5(seq, matrix5)
  except KeyError:
    return float('nan')


def compute_ss_strength_features(exon, genome):
  ss_features = dict(
      maxent_upInt5=maxent5_score(exon, 'up', genome),
      maxent_upInt3=maxent3_score(exon, 'up', genome),
      maxent_dnInt5=maxent5_score(exon, 'dn', genome),
      maxent_dnInt3=maxent3_score(exon, 'dn', genome)
    )
  ss_features['maxent_diff5'] = \
    ss_features["maxent_dnInt5"] - ss_features["maxent_upInt5"]
  ss_features['maxent_diff3'] = \
    ss_features["maxent_upInt3"] - ss_features["maxent_dnInt3"]
  return ss_features


def compute_ss_strength(exon_file_path, genome_file, output_file_path):
  wt_genome = twobitreader.TwoBitFile(genome_file)
  n_exons = sum(1 for _ in gzip.open(exon_file_path, 'rt'))
  with gzip.open(exon_file_path, 'rt') as exon_file:
    with gzip.open(output_file_path, 'wt') as output_file:
      reader = csv.DictReader(exon_file, delimiter="\t")
      writer = csv.DictWriter(
        output_file,
        ["event", "maxent_upInt5", "maxent_upInt3",
         "maxent_dnInt5", "maxent_dnInt3",
         'maxent_diff5', 'maxent_diff3']
      )
      writer.writeheader()

      for exon in tqdm.tqdm(reader, disable=None, total=n_exons):
        exon["upIntStart"] = int(exon["upIntStart"]) - 1
        exon["upIntEnd"] = int(exon["upIntEnd"])
        exon["dnIntStart"] = int(exon["dnIntStart"]) - 1
        exon["dnIntEnd"] = int(exon["dnIntEnd"])

        mt_genome = get_variant_genome_from_dict(exon, wt_genome)
        ss_features = compute_ss_strength_features(exon, mt_genome)
        ss_features["event"] = exon["event"]
        writer.writerow(ss_features)


if __name__ == "__main__":
  compute_ss_strength(snakemake.input['exons'], snakemake.input['genome'], snakemake.output[0])
  # compute_ss_strength(
  #   "CONTROLLED_ACCESS/inputs/MSSNG/putative_microexon_000035518_variants.hg38.tab.gz",
  #   "EXTERNAL_DATA/hg38.2bit",
  #   "CONTROLLED_ACCESS/out.csv.gz"
  # )
