#
# Creates a filtered version of the drug curve table with two outlier cell lines (KMS21BM_JCRB and
# Karpas417_ECACC) excluded.
#
library(tidyverse)

drug_curves <- read_tsv(snakemake@input[[1]], show_col_types=FALSE)

outliers <- c("KMS21BM_JCRB", "Karpas417_ECACC")

drug_curves %>%
  filter(!cell_line %in% outliers) %>%
  write_tsv(snakemake@output[[1]])
