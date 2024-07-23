#
# visualizes drug response curves and curve clusters in a few different ways.
#
# 1. drug similarity umap plot (9 clusters)
# 2. drug similarity umap plot (2 clusters)
# 3. drug AC-50 by cluster/cell line
# 4. drug curves by cluster/all
# 5+. drug curves by cluster/cell line
#
library(tidyverse)
library(uwot)

set.seed(1)

save.image()

# load drug curves
drug_curves <- read_tsv(snakemake@input[[1]])

# load drug similarity matrix UMAP projection
drug_umap <- read_tsv(snakemake@input[[2]], show_col_types=FALSE)

# load drug cluster assignments
drug_clusters <- read_tsv(snakemake@input[[3]], show_col_types=FALSE)
drug_clusters$cluster <- factor(drug_clusters$cluster)

# add column corresponding to membership in the "left" or "right" lobe of the UMAP plot
drug_clusters <- drug_clusters %>%
  mutate(super_cluster=factor(ifelse(cluster %in% c(1, 4, 5, 6, 8), "I", "II")))

# load drug metadata
mdat <- read_tsv(snakemake@input[[4]], col_types=cols(.default="c")) %>%
  select(drug_id=sample_id,
         binary_pharmacology,
         mechanistic_class,
         development_phase,
         general_indication,
         target_class,
         gene_target=gene_symbol_1,
         apoptosis=mechanisms_cellular_apoptosis,
         angiogenesis=mechanisms_cellular_angiogenesis,
         signal_transduction=mechanisms_cellular_signal_transduction)

# create separate columns for some of the interesting drug annotations
mdat$mechanistic_class <- factor(mdat$mechanistic_class)
mdat$apoptosis <- factor(mdat$apoptosis, labels=c("Other", "Apoptosis Inducer"))
mdat$angiogenesis <- factor(mdat$angiogenesis, labels=c("Other", "Angiogenesis Inhibitor"))
mdat$signal_transduction <- factor(mdat$signal_transduction, labels=c("Other", "Signal Transduction Modulator"))

mdat <- mdat %>%
  mutate(cell_cycle = mechanistic_class == "Cell cycle",
         kinase = target_class == "Kinase",
         gpcr = target_class == "GPCR signaling",
         tubb = gene_target == "TUBB",
         pik3ca = gene_target == "PIK3CA",
         egfr = gene_target == "EGFR")

df <- drug_umap %>%
  inner_join(mdat, by="drug_id") %>%
  inner_join(drug_clusters, by="drug_id")

# 1) UMAP plot (color = cluster)
ggplot(df, aes(x=UMAP1, y=UMAP2, color=cluster)) +
  geom_point() +
  ggtitle("HMCL Drug Similarity (UMAP)") +
  guides(fill=guide_legend(title="Cluster")) +
  theme_bw()

ggsave(snakemake@output[[1]], width=1920, height=1080, units="px", dpi=128)

# 2) UMAP plot (color = super cluster)
ggplot(df, aes(x=UMAP1, y=UMAP2, color=super_cluster)) +
  geom_point() +
  ggtitle("HMCL Drug Similarity (UMAP)") +
  guides(fill=guide_legend(title="Cluster")) +
  theme_bw()

ggsave(snakemake@output[[2]], width=1920, height=1080, units="px", dpi=128)

# 3) Mean cell line AC-50 by cluster
ac50_mat <- drug_curves %>%
  select(-plate_id, -slope, -lower_limit, -upper_limit, -lac50, -starts_with("conc"), -starts_with("dose")) %>%
  inner_join(drug_clusters, by="drug_id")

# clip large values to improve resolution
max_ac50 <- quantile(ac50_mat$ac50, 0.99, na.rm=TRUE)
ac50_mat$ac50 <- pmin(ac50_mat$ac50, max_ac50)

average_ac50 <- ac50_mat %>%
  group_by(cluster, cell_line) %>%
  summarize(ac50 = mean(ac50, na.rm=TRUE)) %>%
  select(cell_line, ac50, cluster)

ggplot(average_ac50, aes(x=cell_line, y=ac50)) +
  geom_bar(stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  ggtitle("Mean AC-50 by cluster") +
  facet_wrap(~cluster, ncol=3)

ggsave(snakemake@output[[3]], width=1080, height=1080, units="px", dpi=92)

# 4) drug curves visualized for each cluster, averaged across all cell lines
drug_mat_all <- drug_curves %>%
  select(-cell_line, -plate_id, -slope, -lower_limit, -upper_limit, -ac50, -lac50, -starts_with("conc")) %>%
  group_by(drug_id) %>%
  summarize(across(everything(), mean)) %>%
  inner_join(drug_clusters, by="drug_id") %>%
  select(-super_cluster) %>%
  pivot_longer(!c(drug_id, cluster), names_to="dose", values_to="viability")

# shift dose numbering up one so that the lowest concetration corresponds to "1" and drop the
# "dose_" prefix
drug_mat_all$dose <- factor(as.numeric(sub("dose_", "", drug_mat_all$dose)) + 1)

# compute cluster averages
cluster_means <- drug_mat_all %>%
  group_by(cluster, dose) %>%
  summarize(viability = mean(viability))

cluster_means <- cbind(drug_id="average", cluster_means)

drug_mat_all <- rbind(drug_mat_all, cluster_means)

ggplot(drug_mat_all, aes(x=dose, y=viability, group=drug_id)) +
  geom_line(color="#aaa") +
  geom_line(data=filter(drug_mat_all, drug_id == "average"), aes(x=dose, y=viability), colour="red", linewidth=1) +
  theme_bw() +
  ggtitle("Average dose response curves by cluster (all cells)") +
  facet_wrap(~cluster, ncol=3)

ggsave(snakemake@output[[4]], width=1920, height=1600, units="px", dpi=192)

# 5+) drug curves visualized for each cluster and for each cell line
drug_mat <- drug_curves %>%
  select(-plate_id, -slope, -lower_limit, -upper_limit, -ac50, -lac50, -starts_with("conc")) %>%
  inner_join(drug_clusters, by="drug_id") %>%
  select(-super_cluster)

for (cell in unique(drug_mat$cell_line)) {
  df <- drug_mat %>%
    filter(cell_line == cell) %>%
    select(-cell_line) %>%
    pivot_longer(!c(drug_id, cluster), names_to="dose", values_to="viability")

  # shift dose numbering up one so that the lowest concetration corresponds to "1" and drop the
  # "dose_" prefix
  df$dose <- factor(as.numeric(sub("dose_", "", df$dose)) + 1)

  # compute cluster averages
  cluster_means <- df %>%
    group_by(cluster, dose) %>%
    summarize(viability = mean(viability))

  cluster_means <- cbind(drug_id="average", cluster_means)

  df <- rbind(df, cluster_means)

  ggplot(df, aes(x=dose, y=viability, group=drug_id)) +
    geom_line(color="#aaa") +
    geom_line(data=filter(df, drug_id == "average"), aes(x=dose, y=viability), colour="red", linewidth=1) +
    theme_bw() +
    ggtitle(sprintf("Average drug response by cluster (%s)", cell)) +
    facet_wrap(~cluster, ncol=3)

  outfile <- sub("all_cells", cell, snakemake@output[[4]])
  ggsave(outfile, width=1080, height=1080, units="px", dpi=92)
}

