#
# Clusters drugs in pairwise drug similarity matrix
#
library(tidyverse)

set.seed(1)

# number of clusters to detect?
k <- 9

infile <- snakemake@input[[1]]

cor_mat <- read_tsv(infile, show_col_types=FALSE) %>%
  column_to_rownames("drug_id")

clusters <- kmeans(cor_mat, k)$cluster

res <- data.frame(drug_id=names(clusters),
                  cluster=factor(as.vector(clusters)))

res %>%
  write_tsv(snakemake@output[[1]])
