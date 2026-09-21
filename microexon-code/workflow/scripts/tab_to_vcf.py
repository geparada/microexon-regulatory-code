import csv
import pathlib
from collections import defaultdict
import gzip
import os

VCF_samples = defaultdict(list)

with gzip.open(snakemake.input[0], "rt") as file:
  reader = csv.DictReader(file, delimiter="\t")

  for row in reader:
    if row["FILTER"] == "PASS":
      sample = row["#Sample"]

      vcf_row = (row["CHROM"], row["POS"], row["ID"], row["REF"], row["ALT"], row["QUAL"], row["FILTER"], row["INFO"])

      end = int(row["end"])
      start = int(row["start"])

      VCF_samples[sample].append(vcf_row)

output_path = pathlib.Path(snakemake.output[0])
os.makedirs(output_path, exist_ok=True)

for sample in VCF_samples:

  VCF = VCF_samples[sample]
  VCF_sorted = sorted(VCF, key=lambda x: (x[0], int(x[1])))

  with open(output_path / f"{sample}.A.vcf", 'wt', newline='') as out_file:

    tsv_writer = csv.writer(out_file, delimiter='\t')
    tsv_writer.writerow(["##fileformat=VCFv4.1"])
    tsv_writer.writerow(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample])

    for row in VCF_sorted:
      out = list(row) + ["GT", "1/1"]
      tsv_writer.writerow(out)
