import pybedtools
import pandas as pd
import gzip


def sort_variants(variants, mic_events, metadata, output_file):
  event_ids = pd.Index([event.fields[3] for event in pybedtools.BedTool(mic_events)])
  metadata = pd.read_csv(metadata, sep="\t", index_col="Sample ID")

  with gzip.open(output_file, "wt", newline='') as output_f:
    write_header = True
    for variant_file in variants:
      variants_for_sample = pd.read_csv(variant_file, sep="\t", index_col="#Sample")
      variants_for_sample = variants_for_sample.join(metadata.Affection)
      variants_for_sample.loc[event_ids.intersection(variants_for_sample.index)]

      variants_for_sample.to_csv(output_f, sep="\t", header=write_header, index_label="#Sample")
      write_header = False


if __name__ == "__main__":
  sort_variants(**snakemake.input, output_file=snakemake.output[0])
