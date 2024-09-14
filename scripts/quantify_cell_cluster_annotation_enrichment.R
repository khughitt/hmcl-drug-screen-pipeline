#
# Uses the hypergeometric test to check for over-representation of cell metadata annotations
#
library(tidyverse)

set.seed(1)

# minimum number of times a specific cell metadata variable value (i.e. factor "level") has to
# appear to be considered for enrichment?
MIN_FACTOR_LEVELS <- 5

# load cell similarity matrix
cor_mat <- read_tsv(snakemake@input[[1]], show_col_types=FALSE) %>%
  column_to_rownames("cell")

# load cell cluster assignments
cell_clusters <- read_tsv(snakemake@input[[2]], show_col_types=FALSE)
cell_clusters$cluster <- factor(cell_clusters$cluster)

# load cell metadata
mdat <- read_tsv(snakemake@input[[3]], show_col_types=FALSE) %>%
  select(-cell_name)

# vectors to store result table elements
cluster_nums <- c()
fields <- c()
lvls <- c()
num_in_clusters <- c()
num_with_lvls <- c()
num_without_lvls <- c()
cluster_sizes <- c()
pvals <- c()

# iterate over cell clusters
for (cluster_num in sort(unique(cell_clusters$cluster))) {
  cluster_cells <- cell_clusters %>%
    filter(cluster == cluster_num) %>%
    pull(cell)

  # iterate over cell metadata variables
  for (field in colnames(mdat)[-1]) {
    # get variable values present >= 20 times?
    field_levels <- mdat %>%
      group_by(!!sym(field)) %>%
      summarize(n=n()) %>%
      filter(n >= MIN_FACTOR_LEVELS) %>%
      pull(!!sym(field))

    # iterate over factor levels
    for (lvl in field_levels) {

      # cell annotated with the specified level?
      lvl_cells <- mdat %>%
        filter(!!sym(field) == lvl) %>%
        pull(cell)

      num_in_cluster <- sum(lvl_cells %in% cluster_cells)
      num_with_lvl <- length(lvl_cells)
      num_without_lvl <- nrow(mdat) - num_with_lvl
      cluster_size <- length(cluster_cells)

      pval <- phyper(num_in_cluster - 1, num_with_lvl, num_without_lvl, cluster_size, lower.tail = FALSE)

      # add to result vectors
      cluster_nums <- c(cluster_nums, cluster_num)
      fields <- c(fields, field)
      lvls <- c(lvls, lvl)
      num_in_clusters <- c(num_in_clusters, num_in_cluster)
      num_with_lvls <- c(num_with_lvls, num_with_lvl)
      num_without_lvls <- c(num_without_lvls, num_without_lvl)
      cluster_sizes <- c(cluster_sizes, cluster_size)
      pvals <- c(pvals, pval)
    }
  }
}

# create result table
dat <- data.frame(cluster_num=cluster_nums, field=fields, lvl=lvls, num_in_cluster=num_in_clusters,
                  num_with_lvl=num_with_lvls, num_without_lvl=num_without_lvls,
                  cluster_size=cluster_sizes, pval=pvals)

dat$padj <- p.adjust(dat$pval)

dat <- dat %>%
  arrange(padj)

write_tsv(dat, snakemake@output[[1]])
