#
# Uses the hypergeometric test to check for over-representation of drug metadata annotations
#
library(tidyverse)

set.seed(1)

# minimum number of times a specific drug metadata variable value (i.e. factor "level") has to
# appear to be considered for enrichment?
MIN_FACTOR_LEVELS <- 20

# load drug similarity matrix
cor_mat <- read_tsv(snakemake@input[[1]], show_col_types=FALSE) %>%
  column_to_rownames("drug_id")

# load drug cluster assignments
drug_clusters <- read_tsv(snakemake@input[[2]], show_col_types=FALSE)
drug_clusters$cluster <- factor(drug_clusters$cluster)

# load drug metadata and extract columns of interest
mdat <- read_tsv(snakemake@input[[3]], col_types=cols(.default="c")) %>%
  select(drug_id=sample_id, mechanistic_class, mechanisms_cellular,
         mechanisms_cellular_apoptosis, mechanisms_cellular_angiogenesis,
         mechanisms_cellular_signal_transduction, target_class, gene_symbol_1, moa1,
         binary_pharmacology, general_indication, class)

mdat$class <- factor(mdat$class)

# vectors to store result table elements
cluster_nums <- c()
fields <- c()
lvls <- c()
num_in_clusters <- c()
num_with_lvls <- c()
num_without_lvls <- c()
cluster_sizes <- c()
pvals <- c()

# iterate over drug clusters
for (cluster_num in sort(unique(drug_clusters$cluster))) {
  cluster_drugs <- drug_clusters %>%
    filter(cluster == cluster_num) %>%
    pull(drug_id)

  # iterate over drug metadata variables
  for (field in colnames(mdat)[-1]) {
    # get variable values present >= 20 times?
    field_levels <- mdat %>%
      group_by(!!sym(field)) %>%
      summarize(n=n()) %>%
      filter(n >= MIN_FACTOR_LEVELS) %>%
      pull(!!sym(field))

    # iterate over factor levels
    for (lvl in field_levels) {

      # drug annotated with the specified level?
      lvl_drugs <- mdat %>%
        filter(!!sym(field) == lvl) %>%
        pull(drug_id)

      num_in_cluster <- sum(lvl_drugs %in% cluster_drugs)
      num_with_lvl <- length(lvl_drugs)
      num_without_lvl <- nrow(mdat) - num_with_lvl
      cluster_size <- length(cluster_drugs)

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
