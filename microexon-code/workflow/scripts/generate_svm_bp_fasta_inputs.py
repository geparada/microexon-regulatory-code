"""Extract sequences to generate an input fasta-file for SVM-BP"""


import pandas as pd
import twobitreader
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from variant_genome import Interval, get_variant_genome_from_series


def get_seq_for_acceptor(event, acceptor_name, genome, seq_len=100):
  if event.strand == '+':
    if acceptor_name == 'alt':
      ss_pos = event.upIntEnd
    elif acceptor_name == 'dn':
      ss_pos = event.dnIntEnd

    seq = Seq(genome.get_seq(Interval(event.chrom, event.strand, ss_pos - seq_len, ss_pos, ss_pos, None)))
  else:
    if acceptor_name == 'alt':
      ss_pos = event.upIntStart
    elif acceptor_name == 'dn':
      ss_pos = event.dnIntStart

    seq = Seq(genome.get_seq(Interval(event.chrom, event.strand, ss_pos, ss_pos + seq_len, ss_pos, None)))

  return SeqRecord(seq, id="{}_{}".format(event.name, acceptor_name))


def generate_seqs_mt(events, wt_genome, seq_len):
  for _, event in events.iterrows():
    mt_genome = get_variant_genome_from_series(event, wt_genome)
    yield get_seq_for_acceptor(event, 'alt', mt_genome, seq_len)
    yield get_seq_for_acceptor(event, 'dn', mt_genome, seq_len)


def make_seqs(event_file, genome_2bit_file, output_fasta_file, seq_len):
  exons = pd.read_csv(
    event_file,
    sep="\t", index_col="event",
    converters={
      'upIntStart': lambda pos: int(pos) - 1,
      'dnIntStart': lambda pos: int(pos) - 1,
    }
  )
  wt_genome = twobitreader.TwoBitFile(genome_2bit_file)
  sequences = generate_seqs_mt(exons, wt_genome, seq_len)
  with open(output_fasta_file, "wt") as f:
    SeqIO.write(sequences, f, "fasta")


if __name__ == "__main__":
  make_seqs(*snakemake.input, snakemake.output[0], snakemake.params['seq_len'])