#%%
import pybedtools
from collections import defaultdict


bed = pybedtools.BedTool(snakemake.input[0])

variant_overlaps = defaultdict(list)
genotypes = defaultdict(list)

for interval in bed:
  if interval[7] != "-1":
    variant = f"{interval[6]}:{interval[7]}:{interval[9]}:{interval[10]}"
    variant_overlaps[interval.name].append(variant)

    gt_idx = interval[14].split(":").index("GT")
    genotype = interval[15].split(":")[gt_idx]
    genotypes[interval.name].append(genotype)

with open(snakemake.output[0], "wt") as out_f:
  out_f.write("interval\tvariants\tgenotypes\n")
  for interval, variants in variant_overlaps.items():
    out_f.write(f"{interval}\t{'|'.join(variants)}\t{'|'.join(genotypes[interval])}\n")
