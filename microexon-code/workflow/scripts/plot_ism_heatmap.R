# install.packages(c("heatmaply", "tidyverse", "cetcolor"), repos='http://cran.rstudio.com/')

library("heatmaply")
library("tidyverse")
library("cetcolor")

scores_heatmap <- read.csv(snakemake@input[[1]], row.names = "motif")
# names(scores_heatmap) <- c("(-160, -140", "(-140, -120)", "(-120, -100)", "(-100, -80)", "(-80, -60)",
#                            "(-60, -40)", "(-40, -20)", "(-20, 0)", "exon", "(0, 20)", "(20, 40)", "(40, 60)",
#                            "(60, 80)", "(80, 100)", "(100, 120)", "(120, 140)", "(140, 160)")

names(scores_heatmap) <- c(-150, -130, -110, -90, -70, -50, -30, -10, "exon", 10, 30, 50, 70, 90, 110, 130, 150)

heatmaply(
  scores_heatmap[1:51,],
  colors = rev(cet_pal(n=100, name = "d1")),
  na.value = "grey50",
  limits = c(-.3, .3),
  #colors = viridis(n = 256, alpha = 1, begin = 0, end = 1, option = "viridis"),
  dendrogram = "row",
  file = snakemake@output[[1]],
  plot_method = "plotly",
  column_text_angle=315,
  xlab="Mutation position",
  ylab="Hexamer",
  key.title="Avg. score change\nby mutation",
  colorbar_len=.3,
  fontsize_row=8,
  fontsize_col=8,
)
