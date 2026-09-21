import pybedtools
import pandas as pd
import gzip
import re


def sort_variants(variants, mic_events, novel_mic_list, metadata, output_file):
  event_ids = [event.fields[3] for event in pybedtools.BedTool(mic_events)]
  metadata = pd.read_csv(metadata, sep="\t", index_col="Sample ID")
  novel_mic_list = pd.read_table(novel_mic_list, index_col="event")
  match_idx = re.compile("(putative_microexon_\d{9})")
  event_ids = pd.Index([idx for idx in event_ids if match_idx.search(idx).group(1) in novel_mic_list.index])

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
