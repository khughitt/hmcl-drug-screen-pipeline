#
# Creates a filtered version of the drug plates table with outlier cell lines removed
#
library(tidyverse)

plate_mat <- read_tsv(snakemake@input[[1]], show_col_types=FALSE) %>%
  as.matrix()

plate_mdata <- read_tsv(snakemake@input[[2]], show_col_types=FALSE)

# get outlier plate ids
outlier_ids <- plate_mdata %>%
  filter(cell_line %in% snakemake@config[["outlier_cells"]]) %>%
  pull(plate)

# filter and store result
mask <- !(colnames(plate_mat) %in% outlier_ids)
plate_mat <- plate_mat[, mask]

write_tsv(as.data.frame(plate_mat), snakemake@output[[1]])
