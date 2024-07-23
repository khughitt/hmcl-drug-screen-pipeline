#
# Creates reduced-dimension versions of the adjusted drug screen data with respect to drugs/cells
#
# in order to make use of all of the available data, instead of projecting the drug x cell matrices
# for a single dose or the inferred AC-50 scores, viability measurements for all doses are
# concatenated and used as input for the dimension reduction step.
#
library(tidyverse)
library(uwot)

set.seed(1)

# load concatenated data
cells <- read.delim(snakemake@input[[1]], sep="\t", row.names=1)
drugs <- read.delim(snakemake@input[[2]], sep="\t", row.names=1)

# PCA (cells)
pca <- prcomp(cells, scale=TRUE)
pca_dat <- pca$x[, 1:2]
colnames(pca_dat) <- c("PC1", "PC2")

pca_dat %>%
  as.data.frame() %>%
  rownames_to_column("cell") %>%
  write_tsv(snakemake@output[[1]])

var_explained <- summary(pca)$importance[2, 1:2] * 100

fp <- file(snakemake@output[[2]], "w")
writeLines(as.character(var_explained), fp)
close(fp)

# UMAP (cells)
umap_cells <- umap(cells, n_neighbors=5, n_components=2)
colnames(umap_cells) <- c("UMAP1", "UMAP2")

umap_cells %>%
  as.data.frame() %>%
  rownames_to_column("cell") %>%
  write_tsv(snakemake@output[[3]])

# PCA (drugs)
pca <- prcomp(drugs, scale=TRUE)
pca_dat <- pca$x[, 1:2]
colnames(pca_dat) <- c("PC1", "PC2")

pca_dat %>%
  as.data.frame() %>%
  rownames_to_column("drug_id") %>%
  write_tsv(snakemake@output[[4]])

var_explained <- summary(pca)$importance[2, 1:2] * 100

fp <- file(snakemake@output[[5]], "w")
writeLines(as.character(var_explained), fp)
close(fp)

# UMAP (drugs)
umap_drugs <- umap(drugs, n_neighbors=15, n_components=2)
colnames(umap_drugs) <- c("UMAP1", "UMAP2")

umap_drugs %>%
  as.data.frame() %>%
  rownames_to_column("drug_id") %>%
  write_tsv(snakemake@output[[6]])
