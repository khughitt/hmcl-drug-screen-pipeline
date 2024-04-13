#
# Combines viability measurements for all doses by cell and drug and stores the resulting matrices.
#
library(tidyverse)

# load viability score matrices
infiles <- unlist(snakemake@input)
infiles <- infiles[grepl("dose", infiles)]

drug_mats <- list()

for (infile in infiles) {
  dose <- sub(".tsv", "", basename(infile))
  drug_mats[[dose]] <- read_tsv(infile, show_col_types = FALSE)
}

# store cell/drug ids
cell_lines <- drug_mats[[1]]$cell_line
drug_ids <- colnames(drug_mats[[1]])[-1]

for (dose in names(drug_mats)) {
  drug_mats[[dose]] <- drug_mats[[dose]] %>%
    select(-cell_line) %>%
    as.matrix()
}

# 1. cell-indexed response matrice;
dat <- as.data.frame(do.call(cbind, drug_mats))
rownames(dat) <- cell_lines
colnames(dat) <- sprintf("X%d", seq_len(ncol(dat)))

dat %>%
  rownames_to_column("cell") %>%
  write_tsv(snakemake@output[[1]])

# 2. drug-indexed response matrice
dat <- as.data.frame(t(do.call(rbind, drug_mats)))
rownames(dat) <- drug_ids

dat %>%
  rownames_to_column("drug_id") %>%
  write_tsv(snakemake@output[[2]])
