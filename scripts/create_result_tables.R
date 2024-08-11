#
# create .xlsx spreadsheets for some of the main summary tables
#
library(tidyverse)
library(openxlsx)

# load input datasets
drug_ac50 <- read_tsv(snakemake@input[[1]], show_col_types = FALSE) %>%
  select(cell_line, drug_id, ac50) %>%
  arrange(cell_line, drug_id)

cell_clusters <- read_tsv(snakemake@input[[2]], show_col_types=FALSE) %>%
  rename(cell_line = cell)

drug_clusters <- read_tsv(snakemake@input[[3]], show_col_types=FALSE)

mdat <- read_tsv(snakemake@input[[4]], show_col_types = FALSE) %>%
  select(drug_id=sample_id, drug_name=sample_name)

# 1) create ac-50 table
drug_ac50 <- drug_ac50 %>%
  inner_join(mdat, by="drug_id") %>%
  select(cell_line, drug_name, drug_id, ac50)

openxlsx::write.xlsx(list(ac50=drug_ac50), file=snakemake@output[[1]])

# 2) create cell cluster table
cell_tbls <- list()

for (cell_cluster in sort(unique(cell_clusters$cluster))) {
  cell_tbls[[sprintf("Cell cluster %d", cell_cluster)]] <- cell_clusters %>%
    filter(cluster == cell_cluster) %>%
    select(-cluster) %>%
    arrange(cell_line)
}

openxlsx::write.xlsx(cell_tbls, file=snakemake@output[[2]])

# 3) create drug cluster table
drug_tbls <- list()

drug_clusters <- drug_clusters %>%
  inner_join(mdat, by="drug_id") %>%
  select(drug_name, drug_id, cluster)

for (drug_cluster in sort(unique(drug_clusters$cluster))) {
  drug_tbls[[sprintf("Drug cluster %d", drug_cluster)]] <- drug_clusters %>%
    filter(cluster == drug_cluster) %>%
    select(-cluster) %>%
    arrange(drug_name)
}

openxlsx::write.xlsx(drug_tbls, file=snakemake@output[[3]])
