"""Compute GC content feature"""

import csv
import gzip

import tqdm
import twobitreader

from variant_genome import Interval, get_variant_genome_from_dict


def gc_content_for_region(interval, genome):
  # seq = wt_genome[chrom][start:stop].upper()
  seq = genome.get_seq(interval)
  return sum(s in {'G', 'C'} for s in seq) / len(seq)


def gc_content_up_int(event, genome, intron_window):
  if event["strand"] == '+':
    interval = Interval(
      event["chrom"], event["strand"], event["upIntEnd"] - intron_window, event["upIntEnd"], event["upIntEnd"], None)
    # start = event["upIntEnd"] - intron_window
    # stop = event["upIntEnd"]
  else:
    interval = Interval(
      event["chrom"], event["strand"], event["upIntStart"], event["upIntStart"] + intron_window,
      event["upIntStart"], None)
    # start = event["upIntStart"]
    # stop = event["upIntStart"] + intron_window

  return gc_content_for_region(interval, genome)


def gc_content_dn_int(event, genome, intron_window):
  if event["strand"] == '+':
    interval = Interval(
      event["chrom"], event["strand"], event["dnIntStart"], event["dnIntStart"] + intron_window,
      event["dnIntStart"], None)
    # start = event["dnIntStart"]
    # stop = event["dnIntStart"] + intron_window
  else:
    interval = Interval(
      event["chrom"], event["strand"], event["dnIntEnd"] - intron_window, event["dnIntEnd"],
      event["upIntEnd"], None)
    # start = event["upIntEnd"] - intron_window
    # stop = event["upIntEnd"]

  return gc_content_for_region(interval, genome)


def gc_content_ex(event, genome):
  if event["strand"] == '+':
    interval = Interval(
      event["chrom"], event["strand"], event["upIntEnd"], event["dnIntStart"], event["upIntEnd"], None)
    # start = event["upIntEnd"]
    # stop = event["dnIntStart"]
  else:
    interval = Interval(
      event["chrom"], event["strand"], event["dnIntEnd"], event["upIntStart"], event["upIntStart"], None)
    # start = event["dnIntStart"]
    # stop = event["upIntEnd"]

  return gc_content_for_region(interval, genome)


# exons = pd.read_csv(
#   exon_file_path,
#   sep="\t_up", index_col="event",
#   converters={
#     'upIntStart': lambda pos: int(pos) - 1,
#     'dnIntStart': lambda pos: int(pos) - 1,
#   }
# )

def get_gc_content(exon, genome, intron_window):
  return dict(
    event=exon["event"],
    gc_content_upInt=gc_content_up_int(exon, genome, intron_window),
    gc_content_ex=gc_content_ex(exon, genome),
    gc_content_dnInt=gc_content_dn_int(exon, genome, intron_window)
  )


def gc_content(exon_file_path, genome_file, output_file_path, intron_window):
  n_exons = sum(1 for _ in gzip.open(exon_file_path, 'rt'))
  wt_genome = twobitreader.TwoBitFile(genome_file)

  with gzip.open(exon_file_path, 'rt') as exon_file:
    with gzip.open(output_file_path, 'wt') as output_file:
      reader = csv.DictReader(exon_file, delimiter="\t")
      writer = csv.DictWriter(
        output_file,
        fieldnames=[
          "event", "gc_content_upInt", "gc_content_ex", "gc_content_dnInt"],
      )
      writer.writeheader()

      for i, exon in enumerate(tqdm.tqdm(reader, total=n_exons, disable=None)):
        exon["upIntStart"] = int(exon["upIntStart"]) - 1
        exon["upIntEnd"] = int(exon["upIntEnd"])
        exon["dnIntStart"] = int(exon["dnIntStart"]) - 1
        exon["dnIntEnd"] = int(exon["dnIntEnd"])

        mt_genome = get_variant_genome_from_dict(exon, wt_genome)

        f = get_gc_content(exon, mt_genome, intron_window)
        writer.writerow(f)


if __name__ == "__main__":
  gc_content(snakemake.input['exons'], snakemake.input['genome'], snakemake.output[0], **snakemake.params)
