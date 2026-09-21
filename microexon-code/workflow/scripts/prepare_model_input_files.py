import pandas as pd
import csv
import gzip
from tqdm import tqdm


def get_variant_distance(event):
  chrom, pos, ref, alt = event.variant.split(':')
  pos = int(pos) - 1

  if event.strand == '+':
    start = event.upIntEnd
    end = event.dnIntStart
  else:
    start = event.dnIntEnd
    end = event.upIntStart

  var_end = pos + len(ref)

  if var_end <= start:
    return var_end - start
  elif pos <= start < var_end <= end:
    return 0
  elif pos <= start < end <= var_end:
    return 0
  elif start <= pos < var_end <= end:
    if (pos - start) < (end - var_end):
      return pos - start
    else:
      return var_end - end
  elif start <= pos <= end < var_end:
    return 0
  elif end <= pos:
    return pos - end
  else:
    raise ValueError(f"How did we get here? pos: {pos}, var_end: {var_end}, start: {start}, end: {end}, event: {event.name}")


def main(variant_files, mic_events, metadata, novel_mic_list, output_file):
  mic_events = pd.read_csv(mic_events, sep="\t", index_col="event")
  mic_events = mic_events.loc[mic_events.label == 'hi']

  novel_mic_list = pd.read_table(novel_mic_list, index_col="event")

  metadata = pd.read_csv(metadata, sep="\t", index_col="Sample ID")

  mic_variant_events = []
  for variant_file in tqdm(variant_files):
    with gzip.open(variant_file, "rt") as variant_f:
      reader = csv.DictReader(variant_f, dialect=csv.excel_tab)

      for v in reader:
        for interval_id in v["microexon_interval_identifier"].split(','):
          mic_id, *region_id = interval_id.split('_')

          # Handle putative microexons
          if mic_id == "putative":
            mic_id = f"putative_microexon_{region_id[1]}"
            region_id = region_id[2:]

          if mic_id in mic_events.index:
            mic_variant = mic_events.loc[mic_id].copy()

          elif mic_id in novel_mic_list.index:
            mic_variant = novel_mic_list.loc[mic_id].copy()

          else:
            continue

          var_str = f"{v['CHROM']}:{v['POS']}:{v['REF']}:{v['ALT']}"
          mic_variant["variant"] = var_str

          try:
            if not region_id:
              mic_variant["affected_exon"] = "mic"
              mic_variant["affected_region"] = "ex"
            elif region_id[0] in ("c1", "c2"):
              mic_variant["affected_exon"] = region_id[0]
              try:
                mic_variant["affected_region"] = region_id[1]
              except IndexError:
                mic_variant["affected_region"] = "ex"
            else:
              mic_variant["affected_exon"] = "mic"
              try:
                mic_variant["affected_region"] = region_id[0]
              except IndexError:
                mic_variant["affected_region"] = "ex"
          except:
            print(f"MIC_ID: {mic_id}, Region ID: {region_id}, Interval ID: {interval_id}")
            raise

          mic_variant["sample_id"] = v["#Sample"]
          mic_variant["affection"] = ["non-ASD", "control", "ASD"][int(metadata.loc[v["#Sample"], "Affection"])]

          mic_var_distance = get_variant_distance(mic_variant)
          if mic_variant.strand == '-':
            mic_var_distance *= -1
          mic_variant["mic_var_distance"] = mic_var_distance

          mic_variant.name = f"{interval_id}_{var_str}"
          mic_variant_events.append(mic_variant)


  mic_variant_events = pd.DataFrame(mic_variant_events)
  # mic_variant_events = mic_variant_events[~mic_variant_events.index.duplicated(keep='first')]

  print("De-duplicating variants:")
  idx_counts = mic_variant_events.index.value_counts()
  idx_counts = idx_counts[idx_counts > 1]
  for duplicated_idx, count in tqdm(idx_counts.iteritems(), total=len(idx_counts)):
    events = mic_variant_events.loc[duplicated_idx]
    mic_variant = events.iloc[0].copy()
    mic_variant["sample_id"] = ','.join(events.sample_id.unique())
    mic_variant["affection"] = ','.join(events.affection.unique())
    mic_variant_events.drop(duplicated_idx, inplace=True)
    mic_variant_events = mic_variant_events.append(mic_variant)

  assert mic_variant_events.index.is_unique

  mic_variant_events.to_csv(output_file, sep="\t", index_label="event")


if __name__ == "__main__":
  print(snakemake.input)
  print(snakemake.output)
  main(snakemake.input["variant_files"], snakemake.input["mic_events"],
       snakemake.input["metadata"], snakemake.input["novel_mic_list"],
       snakemake.output[0])
