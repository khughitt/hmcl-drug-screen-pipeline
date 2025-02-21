#
# compute_ac50_summary.R
#
# Quantifies the number of drugs with AC-50 values for 80% or more of the cell lines passing QC
# below a specified threshold.
#
# Specifically, the following bins are calculated:
#
# i.   < 100 nM for >= 80% cell lines
# ii.  < 500 nM for >= 80% cell lines, but > 100nM for one or more cell lines
# iii. < 1000 nM for >= 80% cell lines, but > 500nM for one or more cell lines
# iv.  > 1000 nM for >= 80% cell lines
#
library(tidyverse)

snek <- snakemake

# Percent of cell lines that must meet the specified threshold for a drug to be considered in
# a particular bin?
REQUIRED_CELL_RATIO <- 0.8

ac50 <- read_tsv(snek@input[[1]], show_col_types = FALSE) %>%
  column_to_rownames("cell_line")

# helper functions to count number of drugs falling in each of the chosen bins
bin1 <- function(x) {
  sum(x < 1e-7, na.rm=TRUE)
}
bin2 <- function(x) {
  sum(x < 5e-7, na.rm=TRUE)
}

bin3 <- function(x) {
  sum(x < 1e-6, na.rm=TRUE)
}

bin4 <- function(x) {
  sum(x > 1e-6, na.rm=TRUE)
}

# determine the number of non-missing viability measurements associated with each drug
num_non_missing <- apply(ac50, 2, function(x) {
  sum(!is.na(x))
})

mask1 <- (apply(ac50, 2, bin1) / num_non_missing) >= REQUIRED_CELL_RATIO
mask2 <- (apply(ac50, 2, bin2) / num_non_missing) >= REQUIRED_CELL_RATIO
mask3 <- (apply(ac50, 2, bin3) / num_non_missing) >= REQUIRED_CELL_RATIO
mask4 <- (apply(ac50, 2, bin4) / num_non_missing) >= REQUIRED_CELL_RATIO

mask2 <- mask2 & !mask1
mask3 <- mask3 & !(mask1 | mask2)

# create ac-50 summary table
drug_ids1 <- colnames(ac50)[mask1]
drug_ids2 <- colnames(ac50)[mask2]
drug_ids3 <- colnames(ac50)[mask3]
drug_ids4 <- colnames(ac50)[mask4]

# number of drugs which don't fit nicely into any of the above bins?
num_remaining <- ncol(ac50) - (length(drug_ids1) + length(drug_ids2) +
                               length(drug_ids3) + length(drug_ids4))

# assemble sheets for output spreadsheet
tbl1 <- tibble(
  Criteria=c("i. < 100 nM for >= 80% cell lines",
             "ii. < 500 nM for >= 80% cell lines, but > 100nM for one or more cell lines",
             "iii. < 1000 nM for >= 80% cell lines, but > 500nM for one or more cell lines",
             "iv. > 1000 nM for >= 80% cell lines",
             "v. Remaining drugs"),
  `Number of drugs`=c(length(drug_ids1), length(drug_ids2),
                      length(drug_ids3), length(drug_ids4),
                      num_remaining)
)

write_tsv(tbl1, snek@output[[1]])

