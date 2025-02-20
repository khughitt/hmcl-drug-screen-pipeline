#
# visualizes drug response curves and curve clusters in a few different ways.
#
# 1. drug similarity pca plot
# 2. drug similarity umap plot
# 3. drug AC-50 by cluster/cell line
# 4. drug curves by cluster/all
# 5. median drug AC-50 values, by cell cluster
# 6+. drug curves by cluster/cell line
#
library(tidyverse)
library(uwot)
library(ggrepel)
library(ggh4x)

set.seed(1)

# plot colors
mean_curve_color <- "#777777"
highlight_color1 <- "#f23639"
highlight_color2 <- "#1B49CD"

# colorblind friend palettes
# https://github.com/JLSteenwyk/ggpubfigs
pal_cells <- c("#648FFF", "#FE6100", "#785EF0")
pal_drugs = c("#332288", "#117733", "#CC6677", "#88CCEE", "#999933", "#882255", "#44AA99", "#DDCC77", "#AA4499")

# mm drugs
mm_drugs <- snakemake@config$mm_drugs$ids
mm_drugs <- setNames(mm_drugs, snakemake@config$mm_drugs$names)

# separate mm drugs into "standard of care" drugs, and others
mm_soc_names <- c("Bortezomib", "Carfilzomib", "Pomalidomide", "Lenalidomide", "Dexamethasone")
mm_soc <- mm_drugs[names(mm_drugs) %in% mm_soc_names]

# show some of the non-soc drugs on the right side of plots to improve legibility
mm_non_soc_right_names <- c("Cytarabine", "Dinaciclib", "Doxorubicin", "Ganetespib", "Oprozomib", "Methotrexate")
mm_non_soc_right <- mm_drugs[names(mm_drugs) %in% mm_non_soc_right_names]
mm_non_soc_left <- mm_drugs[!mm_drugs %in% c(mm_soc, mm_non_soc_right)]

mm_non_soc <- c(mm_non_soc_left, mm_non_soc_right)

# load drug curves
drug_curves <- read_tsv(snakemake@input[[1]])

# load cell x drug AC-50 matrix
ac50_mat <- read_tsv(snakemake@input[[2]]) %>%
  column_to_rownames("cell_line")

# load drug similarity matrix PCA projection
drug_pca <- read_tsv(snakemake@input[[3]], show_col_types=FALSE)
pca_var <- as.numeric(readLines(snakemake@input[[4]]))

# load drug similarity matrix UMAP projection
drug_umap <- read_tsv(snakemake@input[[5]], show_col_types=FALSE)

# load drug cluster assignments
drug_clusters <- read_tsv(snakemake@input[[6]], show_col_types=FALSE)
drug_clusters$cluster <- factor(drug_clusters$cluster)

# load cell cluster assignments
cell_clusters <- read_tsv(snakemake@input[[7]], show_col_types=FALSE) %>%
  rename(cell_line=cell, cell_cluster=cluster) %>%
  mutate(cell_cluster=factor(cell_cluster))

# load drug metadata
mdat <- read_tsv(snakemake@input[[8]], col_types=cols(.default="c")) %>%
  select(drug_id=sample_id,
         drug_name=sample_name,
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

# verify drugs listed in config
if (!all(mdat$drug_id[match(names(mm_drugs), mdat$drug_name)] == mm_drugs)) {
  stop("drug id/name mismatch!")
}

# 1) PCA plot (color = cluster)
df_pca <- drug_pca %>%
  inner_join(mdat, by="drug_id") %>%
  inner_join(drug_clusters, by="drug_id")

ggplot(df_pca, aes(x=PC1, y=PC2, color=cluster)) +
  geom_point() +
  ggtitle("HMCL Drug Similarity (PCA)") +
  guides(color=guide_legend(title="Cluster")) +
  xlab(sprintf("PC1 (%0.2f %% variance)", pca_var[1])) +
  ylab(sprintf("PC2 (%0.2f %% variance)", pca_var[2])) +
  scale_color_manual(values=pal_drugs) +
  theme_bw() +
  theme(axis.title=element_text(size=rel(1.45)),
        axis.text=element_text(size=rel(1.2)),
        legend.text=element_text(size=rel(1.2)),
        legend.title=element_text(size=rel(1.2)))

ggsave(snakemake@output[[1]], width=3375, height=1900, units="px", dpi=300, device="tiff")

# 2) UMAP plot (color = cluster)
df_umap <- drug_umap %>%
  inner_join(mdat, by="drug_id") %>%
  inner_join(drug_clusters, by="drug_id")

ggplot(df_umap, aes(x=UMAP1, y=UMAP2, color=cluster)) +
  geom_point() +
  ggtitle("HMCL Drug Similarity (UMAP)") +
  guides(color=guide_legend(title="Cluster")) +
  scale_color_manual(values=pal_drugs) +
  theme_bw() +
  theme(axis.title=element_text(size=rel(1.45)),
        axis.text=element_text(size=rel(1.2)),
        legend.text=element_text(size=rel(1.2)),
        legend.title=element_text(size=rel(1.2)))

ggsave(snakemake@output[[2]], width=3375, height=1900, units="px", dpi=300, device="tiff")

# 3) Mean cell line AC-50 by cluster
drug_curves_ac50 <- drug_curves %>%
  select(-plate_id, -slope, -lower_limit, -upper_limit, -lac50, -ac50_pval, -starts_with("conc"), -starts_with("dose")) %>%
  inner_join(drug_clusters, by="drug_id")

# clip large values to improve resolution
max_ac50 <- quantile(drug_curves_ac50$ac50, 0.99, na.rm=TRUE)
drug_curves_ac50$ac50 <- pmin(drug_curves_ac50$ac50, max_ac50)

average_ac50 <- drug_curves_ac50 %>%
  group_by(cluster, cell_line) %>%
  summarize(ac50=mean(ac50, na.rm=TRUE)) %>%
  select(cell_line, ac50, cluster)

# add cell line clusters and group cells by cluster
average_ac50 <- average_ac50 %>%
  inner_join(cell_clusters, by="cell_line")

ordered_cells <- cell_clusters %>%
  arrange(cell_cluster, cell_line) %>%
  pull(cell_line)

average_ac50$cell_line <- factor(average_ac50$cell_line, levels=ordered_cells)

# show cluster sizes in sub-headings
num_clusts <- snakemake@config$num_drug_clusters
clust_counts <- table(drug_clusters$cluster)
labels <- sprintf("%s (n=%d)", names(clust_counts), clust_counts)
average_ac50$cluster <- factor(average_ac50$cluster, levels=1:num_clusts, labels=labels)

# color facet strips to match cluster colors
#text_x = elem_list_text(colour=c("#000000"))
strip <- strip_themed(
  background_x=elem_list_rect(fill=pal_drugs),
  text_x = list(element_text(colour="#eeeeee"), NULL, NULL, NULL, NULL,
                element_text(colour="#eeeeee"), NULL, NULL, NULL)
)

ggplot(average_ac50, aes(x=cell_line, y=ac50, fill=cell_cluster)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=pal_cells) +
  theme_bw() +
  theme(axis.text=element_text(face="bold"),
        axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=rel(0.75)),
        axis.title=element_text(size=rel(1.45)),
        legend.text=element_text(size=rel(1.2)),
        strip.text=element_text(face="bold", size=rel(1.0))) +
  guides(fill=guide_legend(title="Cell Line\nCluster")) +
  ggtitle("Mean AC-50 by Drug Cluster") +
  xlab("Human Myeloma Cell Line") +
  ylab("AC-50 (nM)") +
  facet_wrap2(~cluster, strip=strip)
ggsave(snakemake@output[[3]], width=3375, height=1900, units="px", dpi=300, device="tiff")

# 4) drug curves visualized for each cluster, averaged across all cell lines
drug_mat_all <- drug_curves %>%
  select(-cell_line, -plate_id, -slope, -lower_limit, -upper_limit, -ac50, -ac50_pval, -lac50, -starts_with("conc")) %>%
  group_by(drug_id) %>%
  summarize(across(everything(), mean)) %>%
  inner_join(drug_clusters, by="drug_id") %>%
  pivot_longer(!c(drug_id, cluster), names_to="dose", values_to="viability")

# shift dose numbering up one so that the lowest concetration corresponds to "1" and drop the
# "dose_" prefix
drug_mat_all$dose <- factor(as.numeric(sub("dose_", "", drug_mat_all$dose)) + 1)

# compute cluster averages
cluster_means <- drug_mat_all %>%
  group_by(cluster, dose) %>%
  summarize(viability=mean(viability))

cluster_means <- cbind(drug_id="average", cluster_means)

drug_mat_all <- rbind(drug_mat_all, cluster_means)

drug_mat_all$cluster <- factor(drug_mat_all$cluster, levels=1:num_clusts, labels=labels)

# add labels for drugs to hightlight
drug_mat_all$label <- names(mm_drugs)[match(drug_mat_all$drug_id, mm_drugs)]

# data subsets to use for labeling
dat_non_soc_left <- drug_mat_all %>%
  filter(drug_id %in% mm_non_soc_left & dose == 1)

dat_soc_left <- drug_mat_all %>%
  filter(drug_id %in% mm_soc & dose == 1)

dat_non_soc_right <- drug_mat_all %>%
  filter(drug_id %in% mm_non_soc_right & dose == 11)

# visualize drug curves by cluster
ggplot(drug_mat_all, aes(x=dose, y=viability, group=drug_id)) +
  geom_line(color="#ccc", linewidth=0.35) +
  geom_line(data=filter(drug_mat_all, drug_id == "average"), aes(x=dose, y=viability),
            colour=mean_curve_color, linewidth=1.65, alpha=0.90) +
  geom_line(data=filter(drug_mat_all, drug_id %in% mm_non_soc), aes(x=dose, y=viability),
            colour=highlight_color2, linewidth=0.7, linetype="solid", alpha=0.9) +
  geom_line(data=filter(drug_mat_all, drug_id %in% mm_soc), aes(x=dose, y=viability),
            colour=highlight_color1, linewidth=0.7, linetype="solid", alpha=0.9) +
  geom_label_repel(data=dat_non_soc_left, color=highlight_color2, fontface="bold",
                   aes(label=label), hjust=0.8, nudge_x=-0.4, na.rm=TRUE, size=1.8, label.padding=0.15) +
  geom_label_repel(data=dat_non_soc_right, color=highlight_color2, fontface="bold",
                   aes(label=label), hjust=0, nudge_x=-0.4, na.rm=TRUE, size=1.8, label.padding=0.15) +
  geom_label_repel(data=dat_soc_left, color=highlight_color1, fontface="bold",
                   aes(label=label), hjust=1, nudge_x=0.4, na.rm=TRUE, size=1.8, label.padding=0.15) +
  ggtitle("Average dose response curves by cluster (all cells)") +
  ylab("Viability (%)") +
  xlab("Dose (nM)") +
  theme_bw() +
  theme(axis.text=element_text(face="bold"),
        axis.text.x=element_text(angle=0, hjust=1, vjust=0.5, size=rel(1.0)),
        axis.title=element_text(size=rel(1.2)),
        legend.text=element_text(size=rel(1.2)),
        strip.text=element_text(face="bold", size=rel(0.75))) +
  facet_wrap2(~cluster, strip=strip)
ggsave(snakemake@output[[4]], width=1920, height=1620, units="px", dpi=300, device="tiff")

# 5) average drug ac-50 (by cell line cluster)
clust1_cells <- cell_clusters %>%
  filter(cell_cluster == 1) %>%
  pull(cell_line)

clust2_cells <- cell_clusters %>%
  filter(cell_cluster == 2) %>%
  pull(cell_line)

clust3_cells <- cell_clusters %>%
  filter(cell_cluster == 3) %>%
  pull(cell_line)

cluster_averages <- apply(ac50_mat[clust1_cells, ], 2, median, na.rm=TRUE) %>%
  enframe("drug_id", "cluster1")

cluster_averages$cluster2 <- apply(ac50_mat[clust2_cells, ], 2, median, na.rm=TRUE)
cluster_averages$cluster3 <- apply(ac50_mat[clust3_cells, ], 2, median, na.rm=TRUE)

cluster_df <- cluster_averages %>%
  pivot_longer(-drug_id, names_to="cluster", values_to="median_ac50")

cluster_df$Cluster <- factor(cluster_df$cluster)

save.image()

ggplot(cluster_df, aes(x=median_ac50, fill=Cluster, group=Cluster)) +
  geom_histogram(alpha=0.65, position="dodge") +
  scale_fill_manual(values=pal_cells) +
  xlab("AC-50 (median)") +
  ylab("Count") +
  ggtitle("Median AC-50 (By cell line cluster)")

ggsave(snakemake@output[[5]], width=3000, height=1700, units="px", dpi=300, device="tiff")

# 6+) drug curves visualized for each cluster and for each cell line
drug_mat <- drug_curves %>%
  select(-plate_id, -slope, -lower_limit, -upper_limit, -ac50, -ac50_pval, -lac50, -starts_with("conc")) %>%
  inner_join(drug_clusters, by="drug_id")

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
    summarize(viability=mean(viability))

  cluster_means <- cbind(drug_id="average", cluster_means)

  df <- rbind(df, cluster_means)

  ggplot(df, aes(x=dose, y=viability, group=drug_id)) +
    geom_line(color="#aaa") +
    geom_line(data=filter(df, drug_id == "average"), aes(x=dose, y=viability), colour="red", linewidth=1) +
    ylab("Viability (%)") +
    xlab("Dose (nM)") +
    theme_bw() +
    ggtitle(sprintf("Average drug response by cluster (%s)", cell)) +
    facet_wrap(~cluster, ncol=3)

  outfile <- sub("tiff", "png", sub("all_cells", cell, snakemake@output[[4]]))
  ggsave(outfile, width=1080, height=1080, units="px", dpi=92)
}

