#
# consruct pairwise drug similarity matrix
#
library(tidyverse)

df <- read_tsv(snakemake@input[[1]], show_col_types=FALSE) %>%
  column_to_rownames("drug_id")

cor_mat <- cor(t(df), method="pearson")

cor_mat %>%
  as.data.frame() %>%
  rownames_to_column("drug_id") %>%
  write_tsv(snakemake@output[[1]])
