#
# consruct pairwise cell similarity matrix
#
library(tidyverse)


df <- read_tsv(snakemake@input[[1]], show_col_types=FALSE) %>%
  column_to_rownames("cell")

cor_mat <- cor(t(df), method="spearman")

cor_mat %>%
  as.data.frame() %>%
  rownames_to_column("cell") %>%
  write_tsv(snakemake@output[[1]])
