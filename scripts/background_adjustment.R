#
# background_adjustment.R
#
library(tidyverse)

plate_mat <- read_tsv(snakemake@input[[1]], show_col_types=FALSE) %>%
  as.matrix()

# constants
PLATE_NUM_ROWS <- 32
NUM_CONTROL_WELLS <- 4

# plate column names
cnames <- c("dmso", "pos1", "pos2", "neg",
            paste0(c("drug4_dose", "drug3_dose"), rep(10:0, each=2)),
            paste0(c("drug2_dose", "drug1_dose"), rep(10:0, each=2)))

# average wells across all plates; each row corresponds to a single <x, y> plate position,
# so here we are just finding the average values for each such position across all plates
# and then reshaping those values into a rectangular matrix with the same shape as the actual plates
background_im <- matrix(apply(plate_mat, 1, mean), nrow=PLATE_NUM_ROWS)
colnames(background_im) <- cnames

# cols 1-4 controls
# subtract column-wise average values from each of the first four columns (controls)
col_averages <- apply(background_im[, 1:NUM_CONTROL_WELLS], 2, mean)
background_im[, 1:NUM_CONTROL_WELLS] <- sweep(background_im[, 1:NUM_CONTROL_WELLS], 2, col_averages, "-")

BASE_COL <- NUM_CONTROL_WELLS + 1

# cols  5 - 48
# for the remaining columns (treated wells), iterate over each individual drug concentration,
# compute a contration-specific background, and substract that from the average
for (conc in 0:10) {
  conc_offset <- 2 * conc

  # column indices for a single concentration;
  # there are 11 concentrations in total, going from Dose10 -> Dose0, from left to
  # right in the form:  D4C10, D3C10, D4C9, D3C9 ... | D2C10, D1C10, D2C9, D1C9 ...
  # where "D4C10" stands for "Drug 4, concentration 10"
  indices <- c(BASE_COL + conc_offset, BASE_COL + 1 + conc_offset,
               BASE_COL + 22 + conc_offset, BASE_COL + 23 + conc_offset)

  # determine average viability for all wells of the specified concentration
  conc_average <- mean(background_im[, indices])

  # background_im[, indices] is a <32, 4> rectangular matrix representing the
  # average plate intensities for all wells at a fixed concentration.

  # compute concentration-specific background and update main background image

  # for each position, determine the difference between the global average
  # viability for the specified concentration, and the average observed at that
  # specific well; this difference represents the background to be adjusted for.
  background_im[, indices] <- background_im[, indices] - conc_average
}

# subtract background from normalized plates
plate_mat <- sweep(plate_mat, 1, as.vector(background_im), "-")

# scale to range [0, 100]
for (i in seq_len(ncol(plate_mat))) {

  plate_min <- min(plate_mat[, i])
  plate_max <- max(plate_mat[, i])

  plate_mat[, i] <- ((plate_mat[, i] - plate_min) / (plate_max - plate_min)) * 100
}

write_tsv(as.data.frame(background_im), snakemake@output[[1]])
write_tsv(as.data.frame(plate_mat), snakemake@output[[2]])
