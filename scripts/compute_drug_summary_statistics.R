#
# compute_drug_summary_statistics.R
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
# The drugs falling into each of the bins are reported in separate sheets in a generated excel
# sheet.
#
# Three additional sheets are created with lists of the drugs with AC-50 values < 100nM for _all_
# cells in each of the three cell line clusters.
#
library(tidyverse)

snek <- snakemake

ac50 <- read_tsv(snek@input[[1]], show_col_types = FALSE) %>%
  column_to_rownames("cell_line")

drug_names <- read_tsv(snek@input[[2]]) %>%
  select(drug_id=sample_id, drug_name=sample_name)

drug_clusters <- read_tsv(snek@input[[3]])

cell_clusters <- read_tsv(snek@input[[4]])

# Percent of cell lines that must meet the specified threshold for a drug to be considered in
# a particular bin?
REQUIRED_CELL_RATIO <- 0.8

# colorblind friend palettes
# https://github.com/JLSteenwyk/ggpubfigs
pal_cells <- c("#648FFF", "#FE6100", "#785EF0")

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

# create drug summary table
drug_ids1 <- colnames(ac50)[mask1]
drug_ids2 <- colnames(ac50)[mask2]
drug_ids3 <- colnames(ac50)[mask3]
drug_ids4 <- colnames(ac50)[mask4]

# number of drugs which don't fit nicely into any of the above bins?
num_remaining <- ncol(ac50) - (length(drug_ids1) + length(drug_ids2) +
                                 length(drug_ids3) + length(drug_ids4))

# for each of the cell line clusters, which drugs have AC-50 scores < 100 nM for all cluster
# members?
c1_cells <- cell_clusters %>%
  filter(cluster == 1) %>%
  pull(cell)

c2_cells <- cell_clusters %>%
  filter(cluster == 2) %>%
  pull(cell)

c3_cells <- cell_clusters %>%
  filter(cluster == 3) %>%
  pull(cell)

c1_counts <- apply(ac50[c1_cells, ], 2, function(x) {
  sum(x < 1e-7, na.rm=TRUE)
})
c2_counts <- apply(ac50[c2_cells, ], 2, function(x) {
  sum(x < 1e-7, na.rm=TRUE)
})
c3_counts <- apply(ac50[c3_cells, ], 2, function(x) {
  sum(x < 1e-7, na.rm=TRUE)
})

num_non_missing_c1 <- apply(ac50[c1_cells, ], 2, function(x) {
  sum(!is.na(x))
})
num_non_missing_c2 <- apply(ac50[c2_cells, ], 2, function(x) {
  sum(!is.na(x))
})
num_non_missing_c3 <- apply(ac50[c3_cells, ], 2, function(x) {
  sum(!is.na(x))
})

c1_drugs <- colnames(ac50)[c1_counts == num_non_missing_c1]
c2_drugs <- colnames(ac50)[c2_counts == num_non_missing_c2]
c3_drugs <- colnames(ac50)[c3_counts == num_non_missing_c3]

# assemble sheets for output spreadsheet
sheet1 <- tibble(
  Criteria=c("i. < 100 nM for >= 80% cell lines",
             "ii. < 500 nM for >= 80% cell lines, but > 100nM for one or more cell lines",
             "iii. < 1000 nM for >= 80% cell lines, but > 500nM for one or more cell lines",
             "iv. > 1000 nM for >= 80% cell lines",
             "v. Remaining drugs"),
  `Number of drugs`=c(length(drug_ids1), length(drug_ids2),
                      length(drug_ids3), length(drug_ids4),
                      num_remaining)
)

# compute min, max, mean, and median ac-50 values for each drug
drug_ac50_stats <- tibble(
  drug_id=colnames(ac50),
  min_ac50=apply(ac50, 2, min, na.rm=TRUE),
  max_ac50=apply(ac50, 2, max, na.rm=TRUE),
  mean_ac50=apply(ac50, 2, mean, na.rm=TRUE),
  median_ac50=apply(ac50, 2, median, na.rm=TRUE)
)

# add human-readable drug names and drug cluster assignments
drug_ac50_stats <- drug_ac50_stats %>%
  left_join(drug_names) %>%
  left_join(drug_clusters) %>%
  select(drug_id, drug_name, cluster, everything()) %>%
  arrange(median_ac50)

sheet2 <- drug_ac50_stats %>%
  filter(drug_id %in% drug_ids1)

sheet3 <- drug_ac50_stats %>%
  filter(drug_id %in% drug_ids2)

sheet4 <- drug_ac50_stats %>%
  filter(drug_id %in% drug_ids3)

sheet5 <- drug_ac50_stats %>%
  filter(drug_id %in% drug_ids4)

sheet6 <- drug_ac50_stats %>%
  filter(!(drug_id %in% c(drug_ids1, drug_ids2, drug_ids3, drug_ids4)))

sheet7 <- drug_ac50_stats %>%
  filter(drug_id %in% c1_drugs)

sheet8 <- drug_ac50_stats %>%
  filter(drug_id %in% c2_drugs)

sheet9 <- drug_ac50_stats %>%
  filter(drug_id %in% c3_drugs)

sheets <- list(
  Summary=sheet1,
  `i. < 100 nM`=sheet2,
  `ii. < 500 nM`=sheet3,
  `iii. < 1000 nM`=sheet4,
  `iv. > 1000 nM`=sheet5,
  `v. Remaining drugs`=sheet6,
  `Cell cluster 1 < 100 nM`=sheet7,
  `Cell cluster 2 < 100 nM`=sheet8,
  `Cell cluster 3 < 100 nM`=sheet9
)

openxlsx::write.xlsx(sheets, file=snek@output[[1]])
