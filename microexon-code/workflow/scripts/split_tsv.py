import gzip
import tqdm


def int_ceil(a, b):
  return (a - 1) // b + 1

n_exons = sum(1 for _ in gzip.open(snakemake.input[0], 'rt'))

output_files = {
  k: gzip.open(v, "wt") for k, v in
  enumerate(snakemake.output)
}
n_files = len(output_files)

n_exons_per_file = int_ceil(n_exons, n_files)

with gzip.open(snakemake.input[0], 'rt') as input_file:
  header = input_file.__next__()

  for f in output_files.values():
    f.write(header)

  for i, line in tqdm.tqdm(enumerate(input_file), total=n_exons, disable=None):
    # start_pos = int(line.split('\t_up')[5])
    file_idx = i // n_exons_per_file
    output_files[file_idx].write(line)

for f in output_files.values():
  f.close()
