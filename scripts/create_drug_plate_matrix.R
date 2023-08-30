#
# create_drug_plate_matrix.R
#
# Creates a matrix containing 1d representations each drug plate.
# This is useful as an intermediate format, and makes it easy to compare plates with each other.
#
library(tidyverse)

#
# constants
#
num_drugs <- 1912
num_plates_per_cell_line <- 15
plate_num_rows <- 32
plate_num_cols <- 48
num_drugs_per_row <- 4
num_control_wells <- 4
num_doses <- 11

plate_drug_base_cols <- c(25, 26, 47, 48)      # columns containing highest concentrations for each drug
plate_drug_rel_offsets <- seq(0, -20, by = -2) # positions of well indices relative to base column

# maximum number of drugs tested per cell line
max_num_drugs_tested <- num_plates_per_cell_line * plate_num_rows * num_drugs_per_row

# load raw plate data
raw_plate_dat <- read_tsv(snakemake@input[[1]], show_col_types = FALSE)

# create a list with date for each plate stored as a separate entry
plates <- list()

# order by plate_id, row, and column
raw_plate_dat <- raw_plate_dat %>%
  arrange(plate, row, col)

for (plate_id in unique(raw_plate_dat$plate)) {
  dat <- raw_plate_dat %>%
    filter(plate == plate_id) %>%
    pull(well_value)
  plates[[plate_id]] <- matrix(dat, plate_num_rows, plate_num_cols, byrow = TRUE)
}

# create a mapping from plate -> cell & drugs
plate_mapping <- raw_plate_dat %>%
  dplyr::select(plate, cell_line, date, experiment, layer_name) %>%
  group_by(plate) %>%
  slice(1)

plate_mapping$date <- factor(plate_mapping$date)

# plates generates on the two earliest dates are missing some control columns
plate_mapping$incomplete_controls <- as.character(plate_mapping$date) %in% c("2013-09-09", "2013-09-18")

# create a simple table showing which drugs & cell lines were used for each plate
plate_mapping <- raw_plate_dat %>%
  group_by(cell_line, sample_id) %>%
  arrange(desc(col)) %>%
  slice(1) %>%
  select(plate, cell_line, drug_id = sample_id) %>%
  ungroup

# for each plate, 128 different drugs were applied; while the 128 drugs were chosen at random, the
# same sets of 128 drugs were then always applied together on a given plate.
# as a result, the drugs can be grouped together into a "drug group", and each plate can be
# associated with a specific group, which can be useful when checking for differences across plates
# due to the specific set of drugs applied.
# note that because the number of drugs tested is not a multiple of 128, there is also one drug
# group / set of plates with fewer drugs applied.
drug_groups <- plate_mapping %>%
  group_by(plate) %>%
  arrange(drug_id) %>%
  summarise(drug_group = paste(drug_id, collapse = ", "))

num_groups <- length(unique(drug_groups$drug_group))

drug_groups$drug_group <- factor(drug_groups$drug_group, labels = sprintf("Drug Group %d",
                                                                          seq_len(num_groups)))

# create a table mapping from "drug group" -> drug
drug_group_members <- list()

for (group_id in sort(unique(drug_groups$drug_group))) {
  plate_id <- drug_groups %>%
    filter(drug_group == group_id) %>%
    head(1) %>%
    pull(plate)

  drug_group_members[[group_id]] <- plate_mapping %>%
    filter(plate == plate_id) %>%
    pull(drug_id) %>%
    unique() %>%
    sort()
}

drug_group_mapping <- stack(drug_group_members) %>%
  select(drug_group = ind, ncgc_id = values) %>%
  filter(ncgc_id != "DMSO")

# next, each plate is scored based on how often viability measurements match an idealized sigmoidal
# drug response curve with increasing dose.
# this can be useful to help flag potentially problematic plates with behavior that deviates
# significantly from this.
plate_scores <- list()

# model using logistic function
x <- -5:5
idealized_viability_scores <- rev(1 / (1 + exp(-x)))

for (plate_id in names(plates)) {
  # create an empty vector to store drug behavior scores for each plate
  plate_scores[[plate_id]] <- c()

  # iterate over rows on plate
  for (row_num in 1:plate_num_rows) {
    # iterate over drugs on row
    for (col_offset in plate_drug_base_cols) {
      # viability scores for a single drug (from lowest -> highest concentration);
      # each drug has 11 doses, and is zebra-striped with another drug
      col_indices <- col_offset + plate_drug_rel_offsets
      viability_scores <- plates[[plate_id]][row_num, col_indices]

      # compute correlation between actual and idealized viability scores
      drug_score <- cor(viability_scores, idealized_viability_scores)

      # store result
      plate_scores[[plate_id]] <- c(plate_scores[[plate_id]], drug_score)
    }
  }
}

# median quality score for 128 drugs on each plate
median_plate_scores <- sapply(plate_scores, median) %>%
  enframe(name = "plate", value = "quality")

# create a plate metadata table
plate_mdata <- plate_mapping %>%
  select(plate, cell_line) %>%
  distinct() %>%
  inner_join(drug_groups, by = "plate") %>%
  inner_join(median_plate_scores, by = "plate")

# next, a combined "all plate" matrix is generated with all well values for a single plate stored
# as a 1d column vector;
# this is a useful format for making comparisons across plates

# initialize an empty matrix to final result
plate_mat <- matrix(0, nrow = plate_num_rows * plate_num_cols, ncol = length(plates))

colnames(plate_mat) <- names(plates)

for (plate_id in names(plates)) {
  plate_mat[, plate_id] <- as.vector(plates[[plate_id]])
}

# write outputs
write_tsv(as.data.frame(plate_mat), snakemake@output[[1]])
write_tsv(plate_mdata, snakemake@output[[2]])
write_tsv(drug_group_mapping, snakemake@output[[3]])
