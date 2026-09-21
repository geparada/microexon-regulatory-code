#%%
import csv
import gzip
import sys
from collections import namedtuple

import numpy as np
import pandas as pd
import torch
import twobitreader
from torch.utils.serialization import load_lua
from tqdm import tqdm

from utils import encode_seq, one_hot
from variant_genome import Interval, get_variant_genome_from_series

#%%
Regions = namedtuple(
  "Regions", ("c1_donor_exon", "c1_donor_intron", "alt_intron_up", "alt_exon",
              "alt_intron_dn", "c2_acceptor_intron", "c2_acceptor_exon"))


def regions_for_exon(exon, window_intron, window_exon):
  if exon.strand == '+':
    return Regions(
      c1_donor_exon=[exon.upIntStart - window_exon, exon.upIntStart, exon.upIntStart],  # start, end, anchor
      c1_donor_intron=[exon.upIntStart, exon.upIntStart + window_intron, exon.upIntStart],
      alt_intron_up=[exon.upIntEnd - window_intron, exon.upIntEnd, exon.upIntEnd],
      alt_exon=[exon.upIntEnd, exon.dnIntStart, None],
      alt_intron_dn=[exon.dnIntStart, exon.dnIntStart + window_intron, exon.dnIntStart],
      c2_acceptor_intron=[exon.dnIntEnd - window_intron, exon.dnIntEnd, exon.dnIntEnd],
      c2_acceptor_exon=[exon.dnIntEnd, exon.dnIntEnd + window_exon, exon.dnIntEnd]
    )
  else:
    return Regions(
      c1_donor_exon=[exon.upIntEnd, exon.upIntEnd + window_exon, exon.upIntEnd],
      c1_donor_intron=[exon.upIntEnd - window_intron, exon.upIntEnd, exon.upIntEnd],
      alt_intron_up=[exon.upIntStart, exon.upIntStart + window_intron, exon.upIntStart],
      alt_exon=[exon.dnIntEnd, exon.upIntStart, None],
      alt_intron_dn=[exon.dnIntEnd - window_intron, exon.dnIntEnd, exon.dnIntEnd],
      c2_acceptor_intron=[exon.dnIntStart, exon.dnIntStart + window_intron, exon.dnIntStart],
      c2_acceptor_exon=[exon.dnIntStart - window_exon, exon.dnIntStart, exon.dnIntStart]
    )


def predict_exon(batch_size, deepsea_model, exon, exon_id, feature_names, genome, input_half, input_size, use_cuda):
  regions = regions_for_exon(exon, 300, 100)
  feature_dict = {"event": exon_id}
  seq_torch = []
  for region, (start, end, anchor) in regions._asdict().items():
    # seq = wt_genome[exon.chrom][start - input_half:end + input_half].upper()
    seq = genome.get_seq(Interval(exon.chrom, exon.strand, start-input_half, end+input_half, anchor, None))
    # if exon.strand == '-':
    #   seq = reverse_complement(seq)

    seq_one_hot = one_hot(encode_seq(seq), np.float32, seqweaver=True)

    seq_batch = [seq_one_hot[i:i + input_size] for i in
                 range(seq_one_hot.shape[0] - input_size)]

    if seq_batch:
      seq_torch.append(torch.Tensor(
        np.expand_dims(np.swapaxes(np.stack(seq_batch), 1, 2), -1)))
  seq_torch = torch.cat(seq_torch, 0)
  torch_preds = []
  for i in range(0, seq_torch.shape[0], batch_size):
    seq_torch_batch = seq_torch[i:i + batch_size]

    if use_cuda:
      seq_torch_batch = seq_torch_batch.cuda()
    torch_preds.append(
      deepsea_model.forward(seq_torch_batch).cpu().numpy())
  torch_preds = np.concatenate(torch_preds, 0)
  i = 0
  for region, (start, end, anchor) in regions._asdict().items():
    j = i + (end - start)

    max_features = torch_preds[i:j].max(0, initial=0)
    feature_dict.update(
      {"{}_{}".format(f, region): v
       for f, v in zip(feature_names, max_features)}
    )
    i = j
  return feature_dict


def compute_pytorch_features(
    exons, genome_file, deepsea_model, input_size,
    feature_names, output_file, batch_size=512, use_cuda=False):

  input_half = input_size // 2

  with gzip.open(output_file, 'wt') as f:
    fieldnames = ["event"] + [
      "{}_{}".format(f, r)  for r in Regions._fields for f in feature_names]

    csvwriter = csv.DictWriter(f, fieldnames)
    csvwriter.writeheader()
    wt_genome = twobitreader.TwoBitFile(genome_file)

    for exon_id, exon in tqdm(exons.iterrows(), total=exons.shape[0]):

      mt_genome = get_variant_genome_from_series(exon, wt_genome)

      try:
        feature_dict = predict_exon(batch_size, deepsea_model, exon, exon_id, feature_names, mt_genome, input_half,
                                    input_size, use_cuda)
        csvwriter.writerow(feature_dict)

      except:
        print("Exception in '{}': {}.".format(exon_id, sys.exc_info()[0]))
        raise


def run_pytorch(
    event_file, genome_file, pytorch_model_file, feature_name_file, output_file,
    use_cuda, batch_size
):

  model = load_lua(pytorch_model_file)
  if use_cuda:
    model.cuda()
  exons = pd.read_csv(
    event_file,
    sep="\t", index_col="event",
    converters={
      'upIntStart': lambda pos: int(pos) - 1,
      'dnIntStart': lambda pos: int(pos) - 1,
    }
  )

  feature_names = [line.rstrip().replace('|', '_')
                   for line in open(feature_name_file)]
  compute_pytorch_features(
    exons, genome_file, model, 1000,
    feature_names, output_file, batch_size, use_cuda
  )


if __name__ == "__main__":
  run_pytorch(*snakemake.input, *snakemake.output, **snakemake.params)
