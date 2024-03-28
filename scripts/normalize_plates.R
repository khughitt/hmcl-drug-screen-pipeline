#
# normalize_data.R
#
library(tidyverse)

#plate_mat <- read_tsv("../output/drug_plates/raw.tsv.gz", show_col_types=FALSE) %>%
plate_mat <- read_tsv(snakemake@input[[1]], show_col_types=FALSE) %>%
  as.matrix()

# plate_mdata <- read_tsv("../output/metadata/plate-metadata.tsv", show_col_types=FALSE)
plate_mdata <- read_tsv(snakemake@input[[2]], show_col_types=FALSE)

out_dir <- dirname(dirname(snakemake@output[[1]]))

# determine list of output plate ids; allows subset of plates to be used during development
plate_ids <- unique(basename(dirname(unlist(snakemake@output))))

#
# constants
#
PLATE_NUM_ROWS <- 32
PLATE_NUM_COLS <- 48
POS_CONTROLS <- 2:3
NEG_CONTROLS <- 4
TRMT_WELLS <- 5:48

# plate column names
cnames <- c("DMSO", "Pos1", "Pos2", "Neg",
            paste0(c("Drug4_dose", "Drug3_dose"), rep(10:0, each=2)),
            paste0(c("Drug2_dose", "Drug1_dose"), rep(10:0, each=2)))

# quantile to clip upper values at
clip_quantile <- 0.995

# for each plate:
#  1. compute mean positive control
#  2. compute mean negitive control
#  3. apply norm eqn to all cells
for (plate_id in plate_ids) {
  out_dir_plate <- file.path(out_dir, plate_id)

  cell_line <- plate_mdata %>%
    filter(plate == plate_id) %>%
    pull(cell_line)

  plate_raw <- matrix(plate_mat[, plate_id], nrow=PLATE_NUM_ROWS, ncol=PLATE_NUM_COLS)

  colnames(plate_raw) <- cnames

  # clip negative values (only affects a very small number of measurements)
  plate_normed <- pmax(plate_raw, 0)

  # clip at 99.5% quantile to reduce impact of outliers
  clip_upper <- quantile(plate_raw, clip_quantile)
  plate_normed <- pmin(plate_raw, clip_upper)

  # median of positive controls (cells + bort, cols 2 & 3)
  median_pos <- median(plate_normed[, POS_CONTROLS])

  # median of negative controls (cells + dmso, col 4)
  median_neg <- median(plate_normed[, NEG_CONTROLS])

  # normalized viability (%)=100 * ( well - pos control ) / ( neg control - pos control)
  plate_normed <- 100 * (plate_normed - median_pos) / (median_neg - median_pos)

  # save raw and normalized plate matrices
  write_tsv(as.data.frame(plate_raw), file.path(out_dir_plate, "01-raw.tsv"))
  write_tsv(as.data.frame(plate_normed), file.path(out_dir_plate, "02-normed.tsv"))
}
