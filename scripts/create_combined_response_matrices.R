#
# Creates reduced-dimension versions of the adjusted drug screen data with respect to drugs/cells
#
# in order to make use of all of the available data, instead of projecting the drug x cell matrices
# for a single dose or the inferred AC-50 scores, viability measurements for all doses are
# concatenated and used as input for the dimension reduction step.
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

cell_lines <- drug_mats[[1]]$cell_line

for (dose in names(drug_mats)) {
  drug_mats[[dose]] <- drug_mats[[dose]] %>%
    select(-cell_line) %>%
    as.matrix()
}

# 1. cell-indexed response matrice
dat <- as.data.frame(do.call(cbind, drug_mats))
write_tsv(dat, snakemake@output[[1]])

# 2. drug-indexed response matrice
dat <- as.data.frame(t(do.call(rbind, drug_mats)))
write_tsv(dat, snakemake@output[[2]])

# pca <- prcomp(t(dat), scale = TRUE)
#
# pca_dat <- pca$x[, 1:3]
# colnames(pca_dat) <- c("PC1", "PC2", "PC3")
#
# compute variance explained
# var_explained <- summary(pca)$importance[2, 1:3] * 100


