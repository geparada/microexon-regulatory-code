from collections import Counter
import csv
import gzip
import pathlib
from concurrent.futures import ThreadPoolExecutor as Executor
from itertools import repeat, product
import pandas as pd
from tqdm import tqdm
from collections import defaultdict


dict_counter = lambda: defaultdict(Counter)


def get_stats_for_file(file_path, regions):
  sample_id = pathlib.Path(file_path).name
  sample_id = sample_id[:sample_id.index(".tsv.gz")]

  with gzip.open(file_path, "rt") as f:
    reader = csv.DictReader(f, dialect=csv.excel_tab)

    counts = dict_counter()

    for variant in reader:
      variant_region = variant["microexon_interval_identifier"]

      for region_set, region_ids in regions.items():
        if variant_region in region_ids:
          counts['variants'][region_set] += 1

          if len(variant["REF"]) == 1 and len(variant["ALT"]) == 1:
            counts['SNVs'][region_set] += 1

          if variant["high_confidence_denovo"] == 'true':
            counts['de_novo'][region_set] += 1

  return sample_id, counts


def read_regions(regions_file):
  name = pathlib.Path(regions_file).name
  name = name[:name.index(".bed.gz")]

  regions = set()
  with gzip.open(regions_file, "rt") as f:
    reader = csv.reader(f, dialect=csv.excel_tab)
    for region in reader:
      regions.add(region[3])
  return name, regions


def main(input_files, region_definitions, metadata, html_output_file, csv_output_file):
  metadata = pd.read_csv(metadata, sep="\t", index_col="Sample ID")

  with Executor() as executor:
    job = executor.map(read_regions, region_definitions)
    regions = dict(job)

    counts = dict(
      affected=dict_counter(),
      unaffected=dict_counter()
    )
    for sample_id, file_counts in \
        tqdm(executor.map(get_stats_for_file, input_files, repeat(regions)), total=len(input_files)):
      if metadata.loc[sample_id, "Affection"] == 1:
        affection = 'unaffected'
      elif metadata.loc[sample_id, "Affection"] == 2:
        affection = 'affected'
      else:
        continue
      counts[affection]['variants'] += file_counts['variants']
      counts[affection]['SNVs'] += file_counts['SNVs']
      counts[affection]['de_novo'] += file_counts['de_novo']

  counts_df = pd.concat([pd.DataFrame(counts['affected']), pd.DataFrame(counts['unaffected'])], axis=1)
  counts_df.columns = pd.MultiIndex.from_product([('affected', 'unaffected'), counts['affected'].keys()])
  counts_df.to_html(html_output_file)
  counts_df.to_csv(csv_output_file)


if __name__ == "__main__":
  # input_files = list(pathlib.Path("CONTROLLED_ACCESS/tcag/filtered_variants/MSSNG_CG").glob("*.tsv.gz"))[:100]
  # region_definitions = list(pathlib.Path("CONTROLLED_ACCESS/MSSNG.variants.2020-07-24").glob("*.bed.gz"))
  # metadata = "CONTROLLED_ACCESS/tcag/metadata/MSSNG_metadata.tsv"
  # output_file = "CONTROLLED_ACCESS/counts.html"
  input_files = snakemake.input["variants"]
  region_definitions = snakemake.input["regions"]
  metadata = snakemake.input["metadata"]
  html_output_file = snakemake.output["html_output"]
  csv_output_file = snakemake.output["csv_output"]
  main(input_files, region_definitions, metadata, html_output_file, csv_output_file)
