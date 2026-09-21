#%%
import os
import pathlib
import re
from itertools import combinations, product

import pandas as pd
import twobitreader
from tqdm import tqdm

from utils import get_sequence_for_event


def complement(nt):
  return 'TGCA'['ACGT'.index(nt)]


def main(
    motifs,
    event_file,
    genome_file,
    output_file,
    event_id,
    motif_len=6
):

  events = pd.read_csv(
    event_file,
    sep="\t",
    converters=dict(
      upIntStart=lambda pos: int(pos) - 1,
      dnIntStart=lambda pos: int(pos) - 1
    ),
    index_col="event",
  )
  events = events.loc[events.label == 'hi']
  event = events.loc[event_id]

  genome = twobitreader.TwoBitFile(genome_file)

  postfix_dict = dict(variants=0)

  seq_up = get_sequence_for_event(event, genome)[0].upper()
  variant_events = []
  j = 0
  for motif_a, motif_b, in combinations(motifs, 2):
    motif_a_matches = [m.start() for m in re.finditer(motif_a, seq_up)]
    motif_b_matches = [m.start() for m in re.finditer(motif_b, seq_up)]

    for motif_a_rel_pos, motif_b_rel_pos in product(motif_a_matches, motif_b_matches):
      if abs(motif_a_rel_pos - motif_b_rel_pos) < motif_len:
        continue

      for motif_a_internal_pos, motif_b_internal_pos in product(range(motif_len), range(motif_len)):

        motif_a_ref_allele = motif_a[motif_a_internal_pos]
        motif_b_ref_allele = motif_b[motif_b_internal_pos]

        motif_a_alt_alleles = set('ACGT').difference(motif_a_ref_allele)
        motif_b_alt_alleles = set('ACGT').difference(motif_b_ref_allele)

        if event.strand == '-':
          motif_a_ref_allele = complement(motif_a_ref_allele)
          motif_b_ref_allele = complement(motif_b_ref_allele)

        for motif_a_alt_allele, motif_b_alt_allele in product(motif_a_alt_alleles, motif_b_alt_alleles):
          event_var = event.copy()

          if event.strand == '+':
            pos_a = event.upIntEnd - 150 + motif_a_rel_pos + motif_a_internal_pos
            pos_b = event.upIntEnd - 150 + motif_b_rel_pos + motif_b_internal_pos

            assert genome[event.chrom][pos_a:pos_a+1].upper() == motif_a_ref_allele
            assert genome[event.chrom][pos_b:pos_b+1].upper() == motif_b_ref_allele

          else:

            motif_a_alt_allele = complement(motif_a_alt_allele)
            motif_b_alt_allele = complement(motif_b_alt_allele)
            pos_a = event.upIntStart + 150 - motif_a_rel_pos - motif_a_internal_pos - 1
            pos_b = event.upIntStart + 150 - motif_b_rel_pos - motif_b_internal_pos - 1

            assert genome[event.chrom][pos_a:pos_a+1].upper() == motif_a_ref_allele
            assert genome[event.chrom][pos_b:pos_b+1].upper() == motif_b_ref_allele

          variant_a = f"{event.chrom}:{pos_a+1}:{motif_a_ref_allele}:{motif_a_alt_allele}"
          variant_b = f"{event.chrom}:{pos_b+1}:{motif_b_ref_allele}:{motif_b_alt_allele}"

          event_var["variant"] = [variant_a, variant_b]
          event_var["motif_a"] = motif_a
          event_var["motif_b"] = motif_b
          event_var["motif_a_rel_pos"] = motif_a_rel_pos - 150
          event_var["motif_b_rel_pos"] = motif_b_rel_pos - 150
          event_var["motif_a_internal_pos"] = motif_a_internal_pos
          event_var["motif_b_internal_pos"] = motif_b_internal_pos
          event_var["motif_a_alt"] = motif_a_alt_allele
          event_var["motif_b_alt"] = motif_b_alt_allele
          event_var.upIntStart += 1
          event_var.dnIntStart += 1

          event_var.name = f"{event.name}_{j}"
          variant_events.append(event_var)

          j += 1
          postfix_dict['variants'] += 1
  variant_events = pd.DataFrame(variant_events)

  if not variant_events.empty:
    variant_events.to_csv(output_file, sep="\t", index_label="event")


if __name__ == '__main__':
  motifs = snakemake.params['motifs']
  event_file = snakemake.input['event_file']
  genome_file = snakemake.input['genome_file']
  output_file = snakemake.output[0]
  event_id = snakemake.wildcards['event_id']
  output_file_pattern = "MicEvents.neural.double_mutants.{event_id}.hg38.tab.gz"

  os.makedirs(output_file, exist_ok=True)
  main(motifs, event_file, genome_file, output_file, event_id)
