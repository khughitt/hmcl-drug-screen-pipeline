#
# Creates drug x cell ac-50, lac-50 & fixed-dose viability matrices
#
# For each of the drug concentrations tested, a matrix is constructed with the adjusted vibility
# scores for those cell/drug combinations.
#
# Similar matrices are constructed based on the AC-50 and log AC-50 values from the fitted drc dose
# response models.
#
library(tidyverse)

# load drug curve data
drug_curves <- read_tsv(snakemake@input[[1]], show_col_types = FALSE)

# fields to generate matrices for
target_fields <- c("ac50", "lac50", paste0("dose_", 0:10))

out_dir <- dirname(snakemake@output[[1]])

# iterate over response variables and generate cell x drug matrices for each
for (response_var in target_fields) {
  drug_mat <- drug_curves %>%
    select(cell_line, drug_id, all_of(response_var)) %>%
    pivot_wider(names_from = drug_id,
                values_from = all_of(response_var))

  outfile <- file.path(out_dir, sprintf("%s.tsv", response_var))
  write_tsv(drug_mat, outfile)
}
