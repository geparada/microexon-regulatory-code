"""Extract splice site features from the COSSMO model."""
import csv
import gzip

import numpy as np
import pandas as pd
import tensorflow as tf
import twobitreader
from tqdm import tqdm

from variant_genome import Interval, get_variant_genome_from_series

FEATURE_NAMES = [
  'event', 'c1_donor_psi', 'a_acceptor_psi', 'a_donor_psi', 'c2_acceptor_psi',
  'c1_donor_logit', 'a_acceptor_logit', 'a_donor_logit','c2_acceptor_logit'
] + [
  "{}_cossmo_feat_{:03d}".format(a, b)
  for a in ('a_acceptor', 'c2_acceptor', 'a_donor', 'c1_donor')
  for b in range(200)
]


def load_cossmo_model(model_path):
  sm = tf.saved_model.load(model_path)

  placeholders = {n.name[n.name.rfind('/') + 1:-2]: n
                  for n in sm.graph.get_collection('inputs')}
  psi_prediction, logits = sm.graph.get_collection('outputs')
  lstm_output = sm.graph.get_tensor_by_name(
    'cossmo_model/scoring_network/output_LSTM/Reshape_1:0')
  cossmo_f = sm.prune(
    [
      placeholders['alt_dna_seq'],
      placeholders['const_dna_seq'],
      placeholders['rna_seq'],
      placeholders['n_alt_ss'],
      placeholders['alt_ss_position'],
      placeholders['const_site_position'],
    ],
    [psi_prediction, logits, lstm_output]
  )

  return cossmo_f


def get_seqs_for_event(exon, window_size, genome):
  """Prepare the input sequences for the COSSMO model from a Vast-Tools
  alternative splicing event."""

  if exon.strand == '+':
    c1_donor = exon.upIntStart
    a_acceptor = exon.upIntEnd
    a_donor = exon.dnIntStart
    c2_acceptor = exon.dnIntEnd
  elif exon.strand == '-':
    c1_donor = exon.upIntEnd
    a_acceptor = exon.upIntStart
    a_donor = exon.dnIntEnd
    c2_acceptor = exon.dnIntStart
  else:
    raise ValueError

  c1_donor_seq = genome.get_seq(
    Interval(exon.chrom, exon.strand, c1_donor-window_size, c1_donor + window_size, c1_donor, None))
  a_acceptor_seq = genome.get_seq(
    Interval(exon.chrom, exon.strand, a_acceptor - window_size, a_acceptor + window_size, a_acceptor, None))
  a_donor_seq = genome.get_seq(
    Interval(exon.chrom, exon.strand, a_donor-window_size, a_donor + window_size, a_donor, None))
  c2_acceptor_seq = genome.get_seq(
    Interval(exon.chrom, exon.strand, c2_acceptor-window_size, c2_acceptor + window_size, c2_acceptor, None))

  # if exon.strand == '-':
  #   c1_donor_seq = reverse_complement(c1_donor_seq)
  #   a_acceptor_seq = reverse_complement(a_acceptor_seq)
  #   a_donor_seq = reverse_complement(a_donor_seq)
  #   c2_acceptor_seq = reverse_complement(c2_acceptor_seq)

  acceptor_inputs = {
    'alt_dna_seq': np.array([[a_acceptor_seq, c2_acceptor_seq]]),
    'const_dna_seq': np.array([c1_donor_seq]),
    'rna_seq': np.array([[
      c1_donor_seq[:window_size] + a_acceptor_seq[window_size:],
      c1_donor_seq[:window_size] + c2_acceptor_seq[window_size:]
    ]]),
    'n_alt_ss': np.array([2], np.int32),
    'alt_ss_position': np.array([[a_acceptor, c2_acceptor]], np.int32),
    'const_site_position': np.array([c1_donor], np.int32)
  }
  acceptor_inputs = {k: tf.constant(v) for k, v in acceptor_inputs.items()}

  donor_inputs = {
    'alt_dna_seq': np.array([[a_donor_seq, c1_donor_seq]]),
    'const_dna_seq': np.array([c2_acceptor_seq]),
    'rna_seq': np.array([[
      a_donor_seq[:window_size] + c2_acceptor_seq[window_size:],
      c1_donor_seq[:window_size] + c2_acceptor_seq[window_size:]
    ]]),
    'n_alt_ss': np.array([2], np.int32),
    'alt_ss_position': np.array([[a_donor, c1_donor]], np.int32),
    'const_site_position': np.array([c2_acceptor], np.int32)
  }
  donor_inputs = {k: tf.constant(v) for k, v in donor_inputs.items()}

  return acceptor_inputs, donor_inputs


def cossmo_features_for_event(exon, acceptor_cossmo_f, donor_cossmo_f, genome):
  """Run COSSMO to get splice site features."""

  acceptor_inputs, donor_inputs = get_seqs_for_event(exon, 100, genome)
  acceptor_psi, acceptor_logits, acceptor_feats = \
    acceptor_cossmo_f(*acceptor_inputs.values())
  donor_psi, donor_logits, donor_feats = \
    donor_cossmo_f(*donor_inputs.values())

  features = np.concatenate(
    [acceptor_feats, donor_feats]).reshape((800,))
  features = np.concatenate(
    [np.stack((donor_psi[0, 0, 1], acceptor_psi[0, 0, 0],
               donor_psi[0, 0, 0], acceptor_psi[0, 0, 1],
               donor_logits[0, 0, 1], acceptor_logits[0, 0, 0],
               donor_logits[0, 0, 0], acceptor_logits[0, 0, 1])),
     features])
  cossmo_features = dict(zip(FEATURE_NAMES[1:], features))
  cossmo_features['event'] = exon.name
  return cossmo_features


def cossmo_predict(
    event_file, acceptor_model_path, donor_model_path, genome_2bit_file, output_file):

  exons = pd.read_csv(
    event_file,
    sep="\t", index_col="event",
    converters={
      'upIntStart': lambda pos: int(pos) - 1,
      'dnIntStart': lambda pos: int(pos) - 1,
    }
  )
  acceptor_cossmo_f = load_cossmo_model(acceptor_model_path)
  donor_cossmo_f = load_cossmo_model(donor_model_path)
  wt_genome = twobitreader.TwoBitFile(genome_2bit_file)
  with gzip.open(output_file, 'wt') as f:
    csvwriter = csv.DictWriter(f, FEATURE_NAMES)
    csvwriter.writeheader()
    for event_id, event in tqdm(exons.iterrows(), total=exons.shape[0]):
      try:
        mt_genome = get_variant_genome_from_series(event, wt_genome)
        cossmo_features_for_exon = cossmo_features_for_event(
          event, acceptor_cossmo_f, donor_cossmo_f, mt_genome)
        csvwriter.writerow(cossmo_features_for_exon)
      except:
        print("Exception in event: {}".format(event_id))
        raise


if __name__ == "__main__":
  cossmo_predict(**snakemake.input, **snakemake.output)