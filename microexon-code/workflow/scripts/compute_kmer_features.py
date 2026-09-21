import csv
import gzip
from collections import Counter

import pandas as pd
import twobitreader
from scipy.stats import fisher_exact
from tqdm import tqdm

from utils import reverse_complement
from variant_genome import Interval, get_variant_genome_from_series

SPECIES = ("hg38", "mm10", "bosTau6", "monDom5", "rheMac2", "panTro4", "rn6")


def get_kmers_from_events(events, upstream_window, k, genome_by_species):
  kmer_counts = Counter()
  for _, event in events.iterrows():
    if event.strand == "+":
      seq = genome_by_species[event.species][event.chrom][
            event.upIntEnd - upstream_window : event.upIntEnd
            ].upper()
    else:
      seq = reverse_complement(
        genome_by_species[event.species][event.chrom][
        event.upIntStart : event.upIntStart + upstream_window
        ].upper()
      )

    for i in range(upstream_window - k):
      kmer = seq[i : i + k]
      if "N" in kmer:
        continue
      kmer_counts[kmer] += 1

  return kmer_counts


def compute_kmer_features_for_event(
    event: pd.Series,
    kmers: pd.DataFrame,
    wt_genome: twobitreader.TwoBitFile,
    bin_width: int,
    upstream_window: int,
):
  mt_genome = get_variant_genome_from_series(event, wt_genome)

  if event.strand == "+":
    upstream_int = Interval(
      event.chrom,
      event.strand,
      event.upIntEnd - upstream_window,
      event.upIntEnd,
      event.upIntEnd,
      None,
      )
  else:
    upstream_int = Interval(
      event.chrom,
      event.strand,
      event.upIntStart,
      event.upIntStart + upstream_window,
      event.upIntStart,
      None,
      )

  upstream_seq = mt_genome.get_seq(upstream_int).upper()
  features = dict(event=event.name)

  for kmer, kmer_details in kmers.iterrows():
    feature_profile = [
      upstream_seq[i : i + len(kmer)] == kmer
      for i in range(0, upstream_window - len(kmer) + 1)
    ]
    kmer_feature = {
      f"kmer_{kmer_details.set}_{kmer}_{i-upstream_window}_{i-upstream_window+bin_width}": any(
        feature_profile[i : i + bin_width]
      )
      for i in range(0, upstream_window, bin_width)
    }
    features.update(kmer_feature)
  return features


def compute_enriched_kmers(
    foreground,
    background,
    upstream_window,
    k,
    genome_by_species,
    n_kmers=None,
    p_value_cutoff=None,
):

  assert n_kmers or p_value_cutoff

  foreground_kmers = pd.Series(
    get_kmers_from_events(foreground, upstream_window, k, genome_by_species),
    name="fg_count",
  )
  background_kmers = pd.Series(
    get_kmers_from_events(background, upstream_window, k, genome_by_species),
    name="bg_count",
  )

  kmers = pd.concat([foreground_kmers, background_kmers], axis=1).fillna(0)

  fg_total = kmers.fg_count.sum()
  bg_total = kmers.bg_count.sum()

  def fisher_exact_test(motif, alternative="greater"):
    cont_table = [
      [int(motif["fg_count"]), int(motif["bg_count"])],
      [int(fg_total - motif["fg_count"]), int(bg_total - motif["bg_count"])],
    ]
    return fisher_exact(cont_table, alternative)[1]

  p_value_high = kmers.apply(
    fisher_exact_test, axis=1, alternative="greater"
  ).sort_values()
  p_value_low = kmers.apply(
    fisher_exact_test, axis=1, alternative="less"
  ).sort_values()

  p_value_high.name = "p_val_min"
  p_value_low.name = "p_val_min"

  p_value_high = pd.concat(p_value_high.align(kmers, join="left"), axis=1)
  p_value_low = pd.concat(p_value_low.align(kmers, join="left"), axis=1)

  p_value_high = p_value_high.astype(dict(fg_count="int", bg_count="int"))
  p_value_low = p_value_low.astype(dict(fg_count="int", bg_count="int"))

  if n_kmers:
    return p_value_high.iloc[:n_kmers], p_value_low.iloc[:n_kmers]
  else:
    return (
      p_value_high.loc[p_value_high < p_value_cutoff],
      p_value_low.loc[p_value_low < p_value_cutoff],
    )


def compute_enriched_kmer_set(
    fg_event_files,
    bg_event_files,
    genome_files,
    output_file,
    upstream_window,
    k,
    n_kmer=None,
    p_value_cutoff=None,
    species_list=SPECIES,
    test_set_id_file=None,
    fg_label=None,
):
  foreground_events = []
  background_events = []

  def get_species_from_filename(fname):
    for s in species_list:
      if s in str(fname):
        return s

  fg_event_files = {get_species_from_filename(f): f for f in fg_event_files}
  bg_event_files = {get_species_from_filename(f): f for f in bg_event_files}
  genome_files = {get_species_from_filename(f): f for f in genome_files}

  if test_set_id_file:
    test_set_ids = pd.Index(
      [id.strip() for id in open(test_set_id_file) if id.strip()]
    )
  else:
    test_set_ids = pd.Index([])

  for species in species_list:
    fg_events = pd.read_csv(
      fg_event_files[species],
      sep="\t",
      index_col="event",
      converters={
        "upIntStart": lambda pos: int(pos) - 1,
        "dnIntStart": lambda pos: int(pos) - 1,
      },
    )
    fg_events["species"] = species
    if fg_label:
      fg_events = fg_events.loc[fg_events.label == fg_label]

    fg_events = fg_events.loc[fg_events.index.difference(test_set_ids)]

    foreground_events.append(fg_events)

    bg_events = pd.read_csv(
      bg_event_files[species],
      sep="\t",
      index_col="event",
      converters={
        "upIntStart": lambda pos: int(pos) - 1,
        "dnIntStart": lambda pos: int(pos) - 1,
      },
    )
    bg_events["species"] = species

    bg_events = bg_events.loc[bg_events.index.difference(test_set_ids)]

    background_events.append(bg_events)

  foreground_events = pd.concat(foreground_events)
  background_events = pd.concat(background_events)

  genomes = {
    species: twobitreader.TwoBitFile(genome_files[species])
    for species in species_list
  }

  kmers_fg, kmers_bg = compute_enriched_kmers(
    foreground_events,
    background_events,
    upstream_window,
    k,
    genomes,
    n_kmer,
    p_value_cutoff,
  )

  kmers_fg["set"] = "fg"
  kmers_bg["set"] = "bg"

  kmers = pd.concat([kmers_fg, kmers_bg], axis=0)

  kmers.to_csv(output_file, index_label="kmer")


def compute_kmer_features(
    event_file, genome_file, kmer_file, output_file, upstream_window, bin_width
):
  events = pd.read_csv(
    event_file,
    sep="\t",
    index_col="event",
    converters={
      "upIntStart": lambda pos: int(pos) - 1,
      "dnIntStart": lambda pos: int(pos) - 1,
    },
  )

  wt_genome = twobitreader.TwoBitFile(genome_file)
  kmers = pd.read_csv(kmer_file, index_col="kmer")

  with gzip.open(output_file, "wt") as output_f:

    csv_writer = None

    for event_name, event in tqdm(events.iterrows(), total=len(events)):
      try:
        kmer_features = compute_kmer_features_for_event(
          event, kmers, wt_genome, bin_width, upstream_window
        )
        if csv_writer is None:
          csv_writer = csv.DictWriter(output_f, fieldnames=kmer_features.keys())
          csv_writer.writeheader()
        csv_writer.writerow(kmer_features)
      except:
        print(f"Error in event '{event_name}':")
        raise


if __name__ == "__main__":
  compute_kmer_features(*snakemake.input, *snakemake.output, **snakemake.params)
