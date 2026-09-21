import twobitreader
import pandas as pd
import re
from utils import reverse_complement


def get_motif_matches_for_re(events, genome, motif_re, window_size=150):
  motifs = []

  for event_id, event in events.iterrows():
    if event.strand == '+':
      upstream_seq = genome[event.chrom][event.upIntEnd - window_size:event.upIntEnd].upper()
      downstream_seq = genome[event.chrom][event.dnIntStart:event.dnIntStart + window_size].upper()
    else:
      upstream_seq = reverse_complement(
        genome[event.chrom][event.upIntStart:event.upIntStart + window_size].upper())
      downstream_seq = reverse_complement(
        genome[event.chrom][event.dnIntEnd - window_size:event.dnIntEnd].upper())

    for match in motif_re.finditer(upstream_seq):
      m_start, m_end = match.span()

      e = event.copy()
      e["motif"] = match.group()
      e["affected_region"] = "up"
      e["m_start"] = m_start - window_size
      e["m_end"] = m_end - window_size

      assert upstream_seq[e.m_start:e.m_end] == e.motif

      motifs.append(e)

    for match in motif_re.finditer(downstream_seq):
      m_start, m_end = match.span()

      e = event.copy()
      e["motif"] = match.group()
      e["affected_region"] = "dn"
      e["m_start"] = m_start
      e["m_end"] = m_end

      assert downstream_seq[e.m_start:e.m_end] == e.motif

      motifs.append(e)

  motifs = pd.DataFrame(motifs)
  return motifs


def get_variants_for_motif(motifs, genome, output_file):
  variants = []

  for event_id, motif in motifs.iterrows():
    for rel_pos, ref_allele in enumerate(motif.motif):

      for alt_allele in {nt for nt in 'ACGT' if nt != ref_allele}:
        variant = motif.copy()
        if motif.strand == '+':
          if motif.affected_region == 'up':
            pos = motif.upIntEnd + motif.m_start + rel_pos
          elif motif.affected_region == 'dn':
            pos = motif.dnIntStart + motif.m_start + rel_pos
          else:
            raise ValueError(motif.affected_region)
          assert genome[motif.chrom][pos].upper() == ref_allele
          variant["variant"] = f"{variant.chrom}:{pos + 1}:{ref_allele}:{alt_allele}"
        else:
          if motif.affected_region == 'up':
            pos = motif.upIntStart - motif.m_start - rel_pos - 1
          elif motif.affected_region == 'dn':
            pos = motif.dnIntEnd - motif.m_start - rel_pos - 1
          else:
            raise ValueError(motif.affected_region)
          assert genome[motif.chrom][pos].upper() == reverse_complement(ref_allele)
          variant["variant"] = \
            f"{variant.chrom}:{pos + 1}:{reverse_complement(ref_allele)}:{reverse_complement(alt_allele)}"

        assert ref_allele != alt_allele

        variant.name = f"{variant.name}_{variant.variant}"
        variants.append(variant)

  variants = pd.DataFrame(variants)
  variants.to_csv(output_file, sep="\t", index_label="event")


def iupac_motif_to_re(motif):
  iupac_dict = dict(
    A="A",
    C="C",
    G="G",
    T="T",
    W="[AT]",
    S="[CG]",
    M="[AC]",
    K="[GT]",
    R="[AG]",
    Y="[CT]",
    B="[CGT]",
    D="[AGT]",
    H="[ACT]",
    V="[ACG]",
    N="[ACGT]"
  )

  return "".join(iupac_dict[nt] for nt in motif)


def main(events, genome, motif, window_size, output_file, event_label_filter=None):
  mic_events = pd.read_csv(
    events,
    sep="\t", index_col="event",
    converters={
      'upIntStart': lambda pos: int(pos) - 1,
      'dnIntStart': lambda pos: int(pos) - 1,
    }
  )

  if event_label_filter:
    mic_events = mic_events.loc[mic_events.label == event_label_filter]

  genome = twobitreader.TwoBitFile(genome)
  motif_re = re.compile(iupac_motif_to_re(motif))

  motifs = get_motif_matches_for_re(mic_events, genome, motif_re, window_size=window_size)
  get_variants_for_motif(motifs, genome, output_file)


if __name__ == "__main__":
  main(
    snakemake.input[0],
    snakemake.input[1],
    snakemake.wildcards["motif"],
    snakemake.params.get("window_size", 150),
    snakemake.output[0],
    snakemake.params.get("event_label_filter")
  )