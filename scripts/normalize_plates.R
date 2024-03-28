#
# normalize_data.R
#
library(tidyverse)

plate_mat <- read_tsv(snakemake@input[[1]], show_col_types=FALSE) %>%
  as.matrix()

# indices of controls in column vector plate representations
POS_CONTROLS <- 33:96
NEG_CONTROLS <- 97:128

# quantile to clip upper values at
clip_quantile <- 0.995

# for each plate:
#  1. compute mean positive control
#  2. compute mean negitive control
#  3. apply norm eqn to all cells
for (i in seq_len(ncol(plate_mat))) {
  # clip negative values (only affects a very small number of measurements)
  plate_mat[, i] <- pmax(plate_mat[, i], 0)

  # clip at 99.5% quantile to reduce impact of outliers
  clip_upper <- quantile(plate_mat[, i], clip_quantile)
  plate_mat[, i] <- pmin(plate_mat[, i], clip_upper)

  # median of positive controls (cells + bort, cols 2 & 3)
  median_pos <- median(plate_mat[POS_CONTROLS, i])

  # median of negative controls (cells + dmso, col 4)
  median_neg <- median(plate_mat[NEG_CONTROLS, i])

  # normalized viability (%)=100 * ( well - pos control ) / ( neg control - pos control)
  plate_mat[, i] <- 100 * (plate_mat[, i] - median_pos) / (median_neg - median_pos)
}

# store normalized version of combined matrix
write_tsv(as.data.frame(plate_mat), snakemake@output[[1]])
