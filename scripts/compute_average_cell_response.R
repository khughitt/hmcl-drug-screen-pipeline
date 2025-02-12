#
# Computes cell line viability averaged across all drugs
#
library(tidyverse)

drug_curves <- read_tsv(snakemake@input[[1]], show_col_types=FALSE)

data_cols <- sprintf("dose_%d", 0:10)

cell_factors <- drug_curves %>%
  select(`Cell Line`=cell_line, all_of(data_cols)) %>%
  pivot_longer(-`Cell Line`, names_to="Dose", values_to="viability") %>%
  group_by(`Cell Line`, Dose) %>%
  summarize(Viability=mean(viability))

cell_factors$Dose <- factor(cell_factors$Dose, levels=paste0("dose_", 0:10), labels=1:11)

cell_factors %>%
  write_tsv(snakemake@output[[1]])
