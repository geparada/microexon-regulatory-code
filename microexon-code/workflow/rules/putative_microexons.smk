rule predict_putative_microexons:
  input:
    expand("results/putative_microexons/shapley_values.putative_microexons.{split:03d}.{species}.csv.gz",
      split=range(500), species=["hg38", "mm10"]),
    expand("results/putative_microexons/predictions.putative_microexons.{split:03d}.{species}.csv.gz",
      split=range(500), species=["hg38", "mm10"]),
