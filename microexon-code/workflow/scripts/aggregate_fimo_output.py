import csv
import gzip
from collections import defaultdict
from subprocess import Popen, PIPE
import sys
import tqdm


output_file = sys.argv[1]

def count_lines_in_gzip_file(file):
  p1 = Popen(["gzip", "-dc", file], stdout=PIPE)
  p2 = Popen(["wc", "-l"], stdin=p1.stdout, stdout=PIPE)
  p1.stdout.close()
  return int(p2.communicate()[0].strip())

# nlines = count_lines_in_gzip_file(snakemake.input[0]) - 1

# df = pd.read_csv(
#   snakemake.input[0],
#   sep="\t_up",
#   dtype={'start': 'Int64', 'stop': 'Int64'},
#   comment='#'
# )
# df['event'] = df.sequence_name.apply(lambda x: x[:x.rindex('_')])
# df['sequence_name'] = df.sequence_name.apply(lambda x: x[x.rindex('_')+1:])
# df.set_index(['event', 'sequence_name', 'motif_alt_id'], inplace=True)

# dfg = df.groupby(level=[0, 1, 2])
# fieldnames = ["event"] + [
#   "{}_{}".format(rbp, region) for region in df.index.levels[1]
#   for rbp in df.index.levels[2]]

# with gzip.open(snakemake.input[0], 'rt') as input_f:


with sys.stdin as input_f:
  csv_reader = csv.DictReader(input_f, delimiter="\t")

  lag_motif = None
  lag_event = None
  lag_region = None
  output_dict = defaultdict(dict)
  max_score = float('-inf')
  fieldnames = set()

  for record in tqdm.tqdm(csv_reader, disable=None):
    event = record['sequence_name'][:record['sequence_name'].rindex('_')]
    region = record['sequence_name'][record['sequence_name'].rindex('_') + 1:]
    score = float(record['score'])

    if region != lag_region:
      if lag_motif is not None:
        # Doesn't_up run in first iteration
        field = "{}_{}".format(lag_motif, lag_region)
        fieldnames.add(field)
        output_dict[lag_event][field] = max_score

      lag_event = event
      lag_region = region
      lag_motif = record['motif_alt_id'] if record['motif_alt_id'] else record['motif_id']
      max_score = float('-inf')

    if score > max_score:
      max_score = float(record['score'])

fieldnames = ["event"] + sorted(list(fieldnames))
with gzip.open(output_file, 'wt') as output_f:
  csv_writer = csv.DictWriter(output_f, fieldnames)
  csv_writer.writeheader()
  for event, features in output_dict.items():
    features['event'] = event
    csv_writer.writerow(features)


# max_score_match = dfg.score.max()
#
# #%%
# out = max_score_match.unstack(2).unstack(1)
# out.columns = ['_'.join(c) for c in out.columns]
# out[out.isna()] = -np.inf
# out.to_csv(snakemake.output[0], index_label='event')
