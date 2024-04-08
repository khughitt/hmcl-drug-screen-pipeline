#
# create_plate_matrices.R
#
# Creates matrices containing 1d representations of the concentrations and measurements for
# each drug plate, as well as some basic metadata relating to the plates and sets of drugs used
# together for the same plates.
#
library(tidyverse)

#
# constants
#
NUM_PLATES_PER_CELL_LINE <- 15
NUM_CONTROLS <- 4
PLATE_CONTROL_INDICES <- 1:4
PLATE_NUM_DRUGS_PER_ROW <- 4
PLATE_NUM_ROWS <- 32
PLATE_NUM_COLS <- 48
PLATE_DRUG_BASE_COLS   <- c(25, 26, 47, 48)  # columns containing highest concentrations for each drug
PLATE_DRUG_REL_OFFSETS <- seq(0, -20, by=-2) # positions of well indices relative to base column

# maximum number of drugs tested per cell line
max_num_drugs_tested <- NUM_PLATES_PER_CELL_LINE * PLATE_NUM_ROWS * PLATE_NUM_DRUGS_PER_ROW

# column types
ctypes <- cols(
  cell_line = col_factor(),
  date = col_date(),
  experiment = col_factor(),
  plate = col_factor(),
  layer_name = col_factor(),
  row = col_integer(),
  col = col_integer(),
  well_value = col_double(),
  sample_id = col_character(),
  concentration = col_double()
)

# load raw plate data
raw_plate_dat <- read_tsv(snakemake@input[[1]], col_types = ctypes)

# create lists to store well & concentration values for each plate
plate_vals <- list()
plate_conc <- list()

# order by plate_id, row, and column
raw_plate_dat <- raw_plate_dat %>%
  arrange(plate, row, col)

for (plate_id in unique(raw_plate_dat$plate)) {
  well_vals <- raw_plate_dat %>%
    filter(plate == plate_id) %>%
    pull(well_value)
  plate_vals[[plate_id]] <- matrix(well_vals, PLATE_NUM_ROWS, PLATE_NUM_COLS, byrow=TRUE)

  conc_vals <- raw_plate_dat %>%
    filter(plate == plate_id) %>%
    pull(concentration)
  plate_conc[[plate_id]] <- matrix(conc_vals, PLATE_NUM_ROWS, PLATE_NUM_COLS, byrow=TRUE)
}

# create a simple table showing which drugs & cell lines were used for each plate
plate_mdata <- raw_plate_dat %>%
  group_by(cell_line, sample_id) %>%
  arrange(desc(col)) %>%
  slice(1) %>%
  select(plate, cell_line, date, drug_id=sample_id) %>%
  ungroup

# for each plate, 128 different drugs were applied; while the 128 drugs were chosen at random, the
# same sets of 128 drugs were then always applied together on a given plate.
# as a result, the drugs can be grouped together into a "drug group", and each plate can be
# associated with a specific group, which can be useful when checking for differences across plates
# due to the specific set of drugs applied.
# note that because the number of drugs tested is not a multiple of 128, there is also one drug
# group / set of plates with fewer drugs applied.
drug_groups <- plate_mdata %>%
  group_by(plate) %>%
  arrange(drug_id) %>%
  summarise(drug_group=paste(drug_id, collapse=", "))

num_groups <- length(unique(drug_groups$drug_group))

drug_groups$drug_group <- factor(drug_groups$drug_group, labels=sprintf("Drug Group %d",
                                                                        seq_len(num_groups)))

# create a table mapping from "drug group" -> drug
drug_group_members <- list()

for (group_id in sort(unique(drug_groups$drug_group))) {
  plate_id <- drug_groups %>%
    filter(drug_group == group_id) %>%
    head(1) %>%
    pull(plate)

  drug_group_members[[group_id]] <- plate_mdata %>%
    filter(plate == plate_id) %>%
    pull(drug_id) %>%
    unique() %>%
    sort()
}

drug_group_mapping <- stack(drug_group_members) %>%
  select(drug_group=ind, ncgc_id=values) %>%
  filter(ncgc_id != "DMSO")

# next, each plate is scored based on how often viability measurements match an idealized sigmoidal
# drug response curve with increasing dose.
# this can be useful to help flag potentially problematic plates with behavior that deviates
# significantly from this.
plate_scores <- list()

# model using logistic function
x <- -5:5
idealized_viability_scores <- rev(1 / (1 + exp(-x)))

for (plate_id in names(plate_vals)) {
  # create an empty vector to store drug behavior scores for each plate
  plate_scores[[plate_id]] <- c()

  # iterate over rows on plate
  for (row_num in 1:PLATE_NUM_ROWS) {
    # iterate over drugs on row
    for (col_offset in PLATE_DRUG_BASE_COLS) {
      # viability scores for a single drug (from lowest -> highest concentration);
      # each drug has 11 doses, and is zebra-striped with another drug
      col_indices <- col_offset + PLATE_DRUG_REL_OFFSETS
      viability_scores <- plate_vals[[plate_id]][row_num, col_indices]

      # compute correlation between actual and idealized viability scores
      drug_score <- cor(viability_scores, idealized_viability_scores)

      # store result
      plate_scores[[plate_id]] <- c(plate_scores[[plate_id]], drug_score)
    }
  }
}

# median quality score for 128 drugs on each plate
median_plate_scores <- sapply(plate_scores, median) %>%
  enframe(name="plate", value="quality")

# add drug gruops and plate scores to plate metadata
plate_mdata <- plate_mdata %>%
  select(plate, date, cell_line) %>%
  distinct() %>%
  inner_join(drug_groups, by="plate") %>%
  inner_join(median_plate_scores, by="plate")

# compute average viability of control wells
control_averages_vec <- c()

for (plate_id in names(plate_vals)) {
  control_averages_vec <- c(control_averages_vec,
                            colMeans(plate_vals[[plate_id]][, PLATE_CONTROL_INDICES]))
}

control_averages <- matrix(control_averages_vec, nrow=length(plate_vals), byrow=TRUE) %>%
  as.data.frame()

colnames(control_averages) <- sprintf("control%d_mean", 1:NUM_CONTROLS)

# convert to z-scores
control_zscores <- as.data.frame(scale(control_averages))

num_stds <- 3
num_outliers <- rowSums(abs(control_zscores) > num_stds)
names(num_outliers) <- names(plate_vals)

# add column indicating the number of control columns with average well values
# three standard deviations away from the column averages across all plates
plate_mdata$num_control_outliers <- num_outliers[match(plate_mdata$plate,
                                                       names(num_outliers))]

# based on the above, plates from one experimental date and two cell lines were found to
# have unexpected control behavior, possibly due to DMSO resistance or a problem distributing the
# drug (cell lines: MM1S_ATCC, U266_ATCC)
plate_mdata$incomplete_controls <- plate_mdata$date == "2013-09-18"

# based downstream investigation of plate images, KMS21BM and Karpas417 were found to be excessively
# sensitive to many drugs, and are flagged as outliers below to bring attention to them.
#
# together, KMS21BM & Karpas417 have the lowest median pairwise correlations with other
# cell lines (0.331 and 0.440 for kms21bm and karpas417, respectively; compared with
# 0.525 or higher for all other cell lines).
#
# plates from both cell lines were found to have significant edge effects.
#
# KMS21BM has also been flagged as a possibly misidentified cell line;
# https://web.expasy.org/cellosaurus/CVCL_2991
#
plate_mdata$outlier_cell_line <- plate_mdata$cell_line %in% c("KMS21BM_JCRB", "Karpas417_ECACC")

# next, "combined" plate measurement and concentration matrices are created where
# each column corresponds to a 1d vector representation of well values or concentrations for
# a single plate.

# initialize empty matrices to well & concentration data tables
well_mat <- matrix(0, nrow=PLATE_NUM_ROWS * PLATE_NUM_COLS, ncol=length(plate_vals))
conc_mat <- matrix(0, nrow=PLATE_NUM_ROWS * PLATE_NUM_COLS, ncol=length(plate_conc))

colnames(well_mat) <- names(plate_vals)
colnames(conc_mat) <- names(plate_conc)

for (plate_id in names(plate_vals)) {
  well_mat[, plate_id] <- as.vector(plate_vals[[plate_id]])
  conc_mat[, plate_id] <- as.vector(plate_conc[[plate_id]])
}

# create a table mapping from drug -> plate + row/column indices of the lowest dose
drug_pos <- raw_plate_dat %>%
  select(plate, cell_line, drug_id=sample_id, concentration, row, col) %>%
  filter(drug_id != "DMSO") %>%
  group_by(cell_line, drug_id) %>%
  slice(which.min(concentration))

# create a table containing the indices and concetrations of each drug dose
drug_inds <- raw_plate_dat %>%
  select(plate, cell_line, drug_id=sample_id, concentration, row, col) %>%
  filter(drug_id != "DMSO") %>%
  group_by(cell_line, drug_id) %>%
  mutate(dose = dense_rank(concentration)) %>%
  arrange(cell_line, drug_id, dose)

# exclude plates outside of requested ones; useful during development
well_mat <- well_mat[, snakemake@params[["plate_ids"]]]
conc_mat <- conc_mat[, snakemake@params[["plate_ids"]]]

drug_inds <- drug_inds %>%
  filter(plate %in% snakemake@params[["plate_ids"]])

plate_mdata <- plate_mdata %>%
  filter(plate %in% snakemake@params[["plate_ids"]])

# write outputs
write_tsv(as.data.frame(well_mat), snakemake@output[[1]])
write_tsv(as.data.frame(conc_mat), snakemake@output[[2]])
write_tsv(plate_mdata, snakemake@output[[3]])
write_tsv(drug_inds, snakemake@output[[4]])
write_tsv(drug_group_mapping, snakemake@output[[5]])
