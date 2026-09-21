storage:
    provider="http",

rule download_genome_2bit:
  """Download the genomes"""
  input:
    lambda wildcards: storage.http(
      f"http://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.species}/bigZips/{wildcards.species}.{wildcards.extension}")
  output:
    "resources/genomes/{species}.{extension}"
  wildcard_constraints:
    species=r"\w{2,6}\d{1,2}",
    extension=r"(2bit)|(fa.gz)"
  shell:
    "mv {input} {output}"
