"""Extract sequences to generate an input fasta-file for FIMO"""
import pandas as pd
import twobitreader
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

from variant_genome import Interval, get_variant_genome_from_series

INTRON_SEQ_LEN = 300


def get_seqs_for_event(event, genome):
  if event.strand == '+':
    acc_pos = event.upIntEnd
    don_pos = event.dnIntStart

    up300_seq = Seq(
      genome.get_seq(Interval(event.chrom, event.strand, acc_pos - INTRON_SEQ_LEN, acc_pos, acc_pos, None))
    )
    exon_seq = Seq(
      genome.get_seq(Interval(event.chrom, event.strand, acc_pos, don_pos, acc_pos, None))
    )
    dn300_seq = Seq(
      genome.get_seq(Interval(event.chrom, event.strand, don_pos, don_pos + INTRON_SEQ_LEN, don_pos, None))
    )
  else:
    acc_pos = event.upIntStart
    don_pos = event.dnIntEnd

    up300_seq = Seq(
      genome.get_seq(Interval(event.chrom, event.strand, acc_pos, acc_pos + INTRON_SEQ_LEN, acc_pos, None))
    )
    exon_seq = Seq(
      genome.get_seq(Interval(event.chrom, event.strand, don_pos, acc_pos, None, None))
    )
    dn300_seq = Seq(
      genome.get_seq(Interval(event.chrom, event.strand, don_pos - INTRON_SEQ_LEN, don_pos, don_pos, None))
    )

  return SeqRecord(up300_seq, id="{}_{}".format(event.name, 'up300')), \
         SeqRecord(exon_seq, id="{}_{}".format(event.name, 'exon')), \
         SeqRecord(dn300_seq, id="{}_{}".format(event.name, 'dn300'))


def generate_seqs_mt(events, wt_genome):
  for _, event in tqdm(events.iterrows(), total=events.shape[0], disable=None):
    mt_genome = get_variant_genome_from_series(event, wt_genome)
    for record in get_seqs_for_event(event, mt_genome):
      yield record


def make_seqs(event_file, genome_2bit_file, output_fasta_file):
  exons = pd.read_csv(
    event_file,
    sep="\t", index_col="event",
    converters={
      'upIntStart': lambda pos: int(pos) - 1,
      'dnIntStart': lambda pos: int(pos) - 1,
    }
  )
  wt_genome = twobitreader.TwoBitFile(genome_2bit_file)
  sequences = generate_seqs_mt(exons, wt_genome)
  with open(output_fasta_file, "wt") as f:
    SeqIO.write(sequences, f, "fasta")


if __name__ == "__main__":
  make_seqs(*snakemake.input, *snakemake.output)