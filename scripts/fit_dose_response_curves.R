#
# fit_dose_response_curves.R
#
library(drc)
library(tidyverse)

set.seed(1)

# constants
PLATE_NUM_ROWS <- 32

# use a four-parameter log-logistic model to fit curves;
curve_params <- c("curve_slope", "curve_lower_lim", "curve_upper_lim", "ac50")
fit_fxn <- LL.4(names=curve_params)

# load data
all_plates_resp <- read_tsv(snakemake@input[[1]], show_col_types=FALSE) %>%
  as.matrix()

drug_inds <- read_tsv(snakemake@input[[2]], show_col_types=FALSE) %>%
  filter(!cell_line %in% snakemake@config$outlier_cells)

# create a list of plate matrices
plate_lst <- list()

# iterate over plates
for (i in seq_len(ncol(all_plates_resp))) {
  plate_id <- colnames(all_plates_resp)[i]

  plate_lst[[plate_id]] <- matrix(all_plates_resp[, i], nrow=PLATE_NUM_ROWS)
}

# next, iterate over cell lines and drugs
cell_drugs <- drug_inds %>%
  select(plate, cell_line, drug_id) %>%
  distinct()

# store output rows in a list of lists
res_lst <- list()

for (i in seq_len(nrow(cell_drugs))) {
  # get well positions for the drug/cell combination
  plate_id <- cell_drugs[i, ]$plate
  cell <- cell_drugs[i, ]$cell_line
  drug <- cell_drugs[i, ]$drug_id

  df <- drug_inds %>%
    filter(plate == plate_id & cell_line == cell & drug_id == drug) %>%
    arrange(dose)

  # get row containing desired well values
  row_ind <- df[1, ]$row
  plate_row <- plate_lst[[plate_id]][row_ind, ]

  # get adjusted viability scores
  col_ind <- df$col
  viability <- plate_row[col_ind]

  # get drug concentrations
  conc <- df$concentration

  # fit dose response model
  res <- tryCatch({
    fit <- drm(viability ~ conc, fct=fit_fxn)
    as.numeric(coef(fit))
  }, error = function(e) {
    rep(NA, 4)
  })

  res_lst[[i]] <- c(cell, drug, plate_id, res, viability, conc)
}

# convert nested list to a tibble
cnames <- c("cell_line", "drug_id", "plate_id",
            "slope", "lower_limit", "upper_limit", "ac50",
            paste0("dose_", 0:10),
            paste0("conc_", 0:10))

curve_df <- res_lst %>%
  map_dfr(set_names, cnames) %>%
  as_tibble() %>%
  mutate(across(4:length(cnames), as.double))

# add LAC-50
curve_df$lac50 <- log10(curve_df$ac50 / 1e6)

curve_df <- curve_df %>%
  dplyr::select(cell_line, drug_id, plate_id, slope, lower_limit, upper_limit, ac50, lac50, everything())

write_tsv(curve_df, snakemake@output[[1]])
