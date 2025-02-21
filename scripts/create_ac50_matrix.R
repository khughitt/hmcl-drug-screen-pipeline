#
# Creates drug x cell ac-50 matrix
#
library(tidyverse)

snek <- snakemake

# load drug curve data
drug_curves <- read_tsv(snek@input[[1]], show_col_types=FALSE)

out_dir <- dirname(snek@output[[1]])

# iterate over response variables and generate cell x drug matrices for each
drug_mat <- drug_curves %>%
  select(cell_line, drug_id, all_of("ac50")) %>%
  pivot_wider(names_from = drug_id, values_from = all_of("ac50"))

write_tsv(drug_mat, snek@output[[1]])
