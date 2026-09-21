import gzip


def main(input_file, output_file, filter_field_name, filter_value):
  with gzip.open(input_file, "rt") as input_file:
      with gzip.open(output_file, "wt") as output_file:
        header = next(input_file)
        header_fields = header.split('\t')
        hq_field_idx = header_fields.index(filter_field_name)
        output_file.write(header)

        for line in input_file:
          line_fields = line.split('\t')
          if line_fields[hq_field_idx] == filter_value:
            output_file.write(line)


if __name__ == "__main__":
  import sys

  if len(sys.argv) == 1:
    main(snakemake.input[0], snakemake.output[0],
         snakemake.params['filter_field_name'], snakemake.params['filter_value'])
  else:
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
