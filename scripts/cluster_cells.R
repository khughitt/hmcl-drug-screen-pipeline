#
# Clusters cells in pairwise cell similarity matrix
#
library(tidyverse)

set.seed(1)

# number of clusters to detect?
k <- snakemake@config[["num_cell_clusters"]]

infile <- snakemake@input[[1]]

cor_mat <- read_tsv(infile, show_col_types=FALSE) %>%
  column_to_rownames("cell")

clusters <- kmeans(cor_mat, k)$cluster

res <- data.frame(cell=names(clusters),
                  cluster=factor(as.vector(clusters)))

res %>%
  write_tsv(snakemake@output[[1]])
