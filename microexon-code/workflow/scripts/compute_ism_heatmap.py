#%%
import pickle

import pandas as pd
import twobitreader
import pathlib
import os
from utils import reverse_complement
import re
from itertools import product
import tqdm
import matplotlib.pyplot as plt
import seaborn as sns

#%%
events = pd.read_csv(
  snakemake.input["events"],
  sep="\t",
  index_col="event",
  converters=dict(
    upIntStart=lambda pos: int(pos) - 1,
    dnIntStart=lambda pos: int(pos) - 1
  )
)

#%%
dpsi = pd.read_csv(
  snakemake.input["ism_score_diff"],
  index_col=["event", "region", "rel_pos", "ref", "alt"],
  ).sort_index(level=[0, 1, 2])

#%%
k = 6
motifs = [''.join(motif) for motif in product('ACGT', repeat=k)]

#%%
genome = twobitreader.TwoBitFile(snakemake.input["genome"])


def get_sequence_for_event(event, window=150):
  if event.strand == '+':
    up_seq = genome[event.chrom][event.upIntEnd-window:event.upIntEnd].upper()
    ex_seq = genome[event.chrom][event.upIntEnd:event.dnIntStart].upper()
    dn_seq = genome[event.chrom][event.dnIntStart:event.dnIntStart+window].upper()
  else:
    up_seq = reverse_complement(genome[event.chrom][event.upIntStart:event.upIntStart+window].upper())
    ex_seq = reverse_complement(genome[event.chrom][event.dnIntEnd:event.upIntStart].upper())
    dn_seq = reverse_complement(genome[event.chrom][event.dnIntEnd-window:event.dnIntEnd].upper())

  return up_seq, ex_seq, dn_seq


def complement_alleles(df):
  compl = dict(A="T", C="G", G="C", T="A")

  df.index = pd.MultiIndex.from_tuples(
    (region, pos, compl[ref], compl[alt]) for region, pos, ref, alt in df.index
  )
  return df


def get_motif_location_for_event(event, dpsi, motif, window=150):
  up_seq, ex_seq, dn_seq = get_sequence_for_event(event, window)

  motif_len = len(motif)

  up_motif_pos = [i - window for i in range(len(up_seq)) if up_seq[i:i+motif_len] == motif]
  ex_motif_pos = [i for i in range(len(ex_seq)) if ex_seq[i:i+motif_len] == motif]
  dn_motif_pos = [i for i in range(len(dn_seq)) if dn_seq[i:i+motif_len] == motif]

  dpsi_event = dpsi.loc[event.name]

  up_dpsi_motif = []
  for pos in up_motif_pos:
    dpsi_pos = dpsi_event.xs(('up', slice(pos, pos+motif_len-1)), level=("region", "rel_pos"), drop_level=False)
    dpsi_pos = dpsi_pos[["wt_psi", "mt_psi", "dpsi"]].rename(columns=dict(
      wt_psi="neural_high_pred_wt",
      mt_psi="neural_high_pred_mt",
      dpsi="neural_high_pred_diff"
    ))
    dpsi_pos["motif_pos"] = pos
    dpsi_pos["in_motif_pos"] = dpsi_pos.index.get_level_values(1) - pos

    if event.strand == '-':
      dpsi_pos = complement_alleles(dpsi_pos)
    up_dpsi_motif.append(dpsi_pos)
  if up_dpsi_motif:
    up_dpsi_motif = pd.concat(up_dpsi_motif)
    up_dpsi_motif.index = pd.MultiIndex.from_tuples((
      (event.name, region, pos, ref, alt) for region, pos, ref, alt in up_dpsi_motif.index),
      names=("event", "region", "position", "ref", "alt"))
  else:
    up_dpsi_motif = None

  ex_dpsi_motif = []
  for pos in ex_motif_pos:
    dpsi_pos = dpsi_event.xs(('ex', slice(pos+1, pos+motif_len)), level=("region", "rel_pos"), drop_level=False)
    dpsi_pos = dpsi_pos[["wt_psi", "mt_psi", "dpsi"]].rename(columns=dict(
      wt_psi="neural_high_pred_wt",
      mt_psi="neural_high_pred_mt",
      dpsi="neural_high_pred_diff"
    ))
    dpsi_pos["motif_pos"] = pos
    dpsi_pos["in_motif_pos"] = dpsi_pos.index.get_level_values(1) - 1 - pos
    if event.strand == '-':
      dpsi_pos = complement_alleles(dpsi_pos)
    ex_dpsi_motif.append(dpsi_pos)
  if ex_dpsi_motif:
    ex_dpsi_motif = pd.concat(ex_dpsi_motif)
    ex_dpsi_motif.index = pd.MultiIndex.from_tuples((
      (event.name, region, pos, ref, alt) for region, pos, ref, alt in ex_dpsi_motif.index),
      names=("event", "region", "position", "ref", "alt"))
  else:
    ex_dpsi_motif = None

  dn_dpsi_motif = []
  for pos in dn_motif_pos:
    dpsi_pos = dpsi_event.xs(('dn', slice(pos+1, pos+motif_len)), level=("region", "rel_pos"), drop_level=False)
    dpsi_pos = dpsi_pos[["wt_psi", "mt_psi", "dpsi"]].rename(columns=dict(
      wt_psi="neural_high_pred_wt",
      mt_psi="neural_high_pred_mt",
      dpsi="neural_high_pred_diff"
    ))
    dpsi_pos["motif_pos"] = pos
    dpsi_pos["in_motif_pos"] = dpsi_pos.index.get_level_values(1) - 1 - pos
    if event.strand == '-':
      dpsi_pos = complement_alleles(dpsi_pos)
    dn_dpsi_motif.append(dpsi_pos)
  if dn_dpsi_motif:
    dn_dpsi_motif = pd.concat(dn_dpsi_motif)
    dn_dpsi_motif.index = pd.MultiIndex.from_tuples((
      (event.name, region, pos, ref, alt) for region, pos, ref, alt in dn_dpsi_motif.index),
      names=("event", "region", "position", "ref", "alt"))
  else:
    dn_dpsi_motif = None

  return up_dpsi_motif, ex_dpsi_motif, dn_dpsi_motif


#%%
up_dpsi, ex_dpsi, dn_dpsi = {}, {},  {}

for motif in tqdm.tqdm(motifs):
  motif_up_dpsi, motif_ex_dpsi, motif_dn_dpsi = [], [], []
  for event_name in dpsi.index.unique(0):
    up_dpsi_motif, ex_dpsi_motif, dn_dpsi_motif = \
      get_motif_location_for_event(events.loc[event_name], dpsi, motif)

    if up_dpsi_motif is not None:
      motif_up_dpsi.append(up_dpsi_motif)
    if ex_dpsi_motif is not None:
      motif_ex_dpsi.append(ex_dpsi_motif)
    if dn_dpsi_motif is not None:
      motif_dn_dpsi.append(dn_dpsi_motif)

  if motif_up_dpsi:
    motif_up_dpsi = pd.concat(motif_up_dpsi)
    motif_by_pos_up = (
      motif_up_dpsi
        .groupby(["motif_pos", "in_motif_pos"])
        .neural_high_pred_diff
        .mean()
        .groupby(level=0)
        .mean()
    )
    up_dpsi[motif] = motif_by_pos_up

  if motif_ex_dpsi:
    motif_ex_dpsi = pd.concat(motif_ex_dpsi)
    motif_by_pos_ex = (
      motif_ex_dpsi
        .groupby(["motif_pos", "in_motif_pos"])
        .neural_high_pred_diff
        .mean()
        .groupby(level=0)
        .mean()
    )
    ex_dpsi[motif] = motif_by_pos_ex

  if motif_dn_dpsi:
    motif_dn_dpsi = pd.concat(motif_dn_dpsi)
    motif_by_pos_dn = (
      motif_dn_dpsi
        .groupby(["motif_pos", "in_motif_pos"])
        .neural_high_pred_diff
        .mean()
        .groupby(level=0)
        .mean()
    )
    dn_dpsi[motif] = motif_by_pos_dn

#%%
up_dpsi_all = pd.DataFrame.from_dict(up_dpsi).T
ex_dpsi_all = pd.DataFrame.from_dict(ex_dpsi).T
dn_dpsi_all = pd.DataFrame.from_dict(dn_dpsi).T

# up_dpsi_all.to_csv(snakemake.output["ism_by_pos_up"])
# ex_dpsi_all.to_csv(snakemake.output["ism_by_pos_ex"])
# dn_dpsi_all.to_csv(snakemake.output["ism_by_pos_dn"])

#%%
up_dpsi_all_long = up_dpsi_all\
  .reset_index()\
  .melt(id_vars='index', var_name="position", value_name="neural_high_pred_diff")
up_dpsi_all_long["position_interval"] = pd.cut(
  up_dpsi_all_long.position,
  [-160, -140, -120, -100, -80, -60, -40, -20, 0],
)
up_dpsi_heatmap = up_dpsi_all_long.groupby(by=["index", "position_interval"]).neural_high_pred_diff.mean().unstack(1)
up_dpsi_heatmap.columns = up_dpsi_heatmap.columns.astype(str)

ex_dpsi_all_long = ex_dpsi_all\
  .reset_index()\
  .melt(id_vars='index', var_name="position", value_name="neural_high_pred_diff")
ex_dpsi_all_long["position_interval"] = pd.cut(
  ex_dpsi_all_long.position + 1, [0, 27], labels=["ex"]).astype(str)
ex_dpsi_heatmap = ex_dpsi_all_long.groupby(by=["index", "position_interval"]).neural_high_pred_diff.mean().unstack(1)
ex_dpsi_heatmap.columns = ex_dpsi_heatmap.columns.astype(str)

dn_dpsi_all_long = dn_dpsi_all\
  .reset_index()\
  .melt(id_vars='index', var_name="position", value_name="neural_high_pred_diff")
dn_dpsi_all_long["position_interval"] = pd.cut(
  dn_dpsi_all_long.position + 1,
  [0, 20, 40, 60, 80, 100, 120, 140, 160],
)
dn_dpsi_heatmap = dn_dpsi_all_long.groupby(by=["index", "position_interval"]).neural_high_pred_diff.mean().unstack(1)
dn_dpsi_heatmap.columns = dn_dpsi_heatmap.columns.astype(str)

dpsi_heatmap = up_dpsi_heatmap.join(ex_dpsi_heatmap).join(dn_dpsi_heatmap)
dpsi_heatmap = dpsi_heatmap.loc[dpsi_heatmap.mean(axis=1).sort_values().index]


#%%
pvals_up, pvals_dn = pickle.load(open(snakemake.input["hexamer_pvals"], "rb"))
pvals_up_bins = pd.cut(
  pvals_up.index,
  [-150, -140, -120, -100, -80, -60, -40, -20, 0]
)
pvals_up_binned = pvals_up.groupby(pvals_up_bins).mean()
pvals_up_min = pvals_up_binned.T.min(axis=1)

pvals_dn = pvals_dn.iloc[1:]
pvals_dn_bins = pd.cut(
  pvals_dn.index + 1,
  [0, 20, 40, 60, 80, 100, 120, 140, 150],
)
pvals_dn_binned = pvals_dn.groupby(pvals_dn_bins).mean()
pvals_dn_min = pvals_dn_binned.T.min(axis=1)

pvals_up_binned.index = pvals_up_binned.index.astype(str)
pvals_dn_binned.index = pvals_dn_binned.index.astype(str)
pvals = pvals_up_binned.T.join(pvals_dn_binned.T)
pvals_sorted = pvals.loc[pvals.min(axis=1).sort_values().index]

#%%
pvals_sorted.to_csv(snakemake.output["motif_pvals"], index_label="motif")
dpsi_heatmap.to_csv(snakemake.output["neural_high_pred_diff_heatmap"], index_label="motif")
dpsi_heatmap.loc[pvals_sorted.index.intersection(dpsi_heatmap.index)].to_csv(
  snakemake.output["neural_high_pred_diff_heatmap_by_pval"], index_label="motif")
