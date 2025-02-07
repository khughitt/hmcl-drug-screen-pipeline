#
# cell visualizations.
#
# 1. cell similarity pca plot
# 2. drug similarity umap plot
# 3. average cell viabilities (color = cell line)
# 4. average cell viabilities (color = cluster)
#
library(tidyverse)
library(ggrepel)

set.seed(1)

# colorblind friend palettes
# https://github.com/JLSteenwyk/ggpubfigs
pal_cells <- c("#648FFF", "#FE6100", "#785EF0")

# load cell similarity matrix PCA and UMAP projections
cell_pca <- read_tsv(snakemake@input[[1]], show_col_types=FALSE)
pca_var <- as.numeric(readLines(snakemake@input[[2]]))

cell_umap <- read_tsv(snakemake@input[[3]], show_col_types=FALSE)

# load cell cluster assignments
cell_clusters <- read_tsv(snakemake@input[[4]], show_col_types=FALSE) %>%
  mutate(cluster = factor(cluster))

# most common concentrations for each dose (nM)
drug_conc <- c(0.779, 2.34, 7.01, 21.0, 63.1, 189.0, 568, 1700, 5116, 15300, 46000)

# load cell average viabilities
cell_factors <- read_tsv(snakemake@input[[5]], show_col_types=FALSE)
cell_factors$Dose <- factor(cell_factors$Dose, levels=1:11, labels=as.character(drug_conc))

# load cell metadata
cell_mdat <- read_tsv(snakemake@input[[6]], show_col_types=FALSE)

cell_names <- cell_mdat %>%
  select(cell, cell_name)

# 1) PCA plot
df_pca <- cell_pca %>%
  inner_join(cell_clusters, by="cell") %>%
  inner_join(cell_names, by="cell")

ggplot(df_pca, aes(x=PC1, y=PC2, color=cluster)) +
  geom_point(size=4.0) +
  geom_label_repel(aes(label=cell_name),
                   fontface="bold", hjust=0.8, nudge_x=-0.4, na.rm=TRUE, size=2.2, label.padding=0.15) +
  ggtitle("HMCL Cell Similarity (PCA)") +
  guides(color=guide_legend(title="Cluster")) +
  xlab(sprintf("PC1 (%0.2f %% variance)", pca_var[1])) +
  ylab(sprintf("PC2 (%0.2f %% variance)", pca_var[2])) +
  scale_color_manual(values=pal_cells) +
  theme_bw() +
  theme(axis.text=element_text(size=rel(1.45)),
        axis.title=element_text(size=rel(1.2)),
        legend.text=element_text(size=rel(1.2)),
        legend.title=element_text(size=rel(1.2)))

ggsave(snakemake@output[[1]], width=1440, height=810, units="px", dpi=128)

# 2) UMAP plot
df_umap <- cell_umap %>%
  inner_join(cell_clusters, by="cell") %>%
  inner_join(cell_names, by="cell")

ggplot(df_umap, aes(x=UMAP1, y=UMAP2, color=cluster)) +
  geom_point(size=4.0) +
  geom_label_repel(aes(label=cell_name),
                   fontface="bold", hjust=0.8, nudge_x=-0.4, na.rm=TRUE, size=2.2, label.padding=0.15) +
  ggtitle("HMCL Cell Similarity (UMAP)") +
  guides(color=guide_legend(title="Cluster")) +
  scale_color_manual(values=pal_cells) +
  theme_bw() +
  theme(axis.text=element_text(size=rel(1.3)),
        axis.title=element_text(size=rel(1.3)),
        legend.text=element_text(size=rel(1.3)),
        legend.title=element_text(size=rel(1.2)))

ggsave(snakemake@output[[2]], width=1440, height=810, units="px", dpi=128)

# 3) average cell viability plot

# assign colors based on average viability
cell_order <- cell_factors %>%
  group_by(`Cell Line`) %>%
  summarize(mean_viability=mean(Viability)) %>%
  arrange(-mean_viability) %>%
  pull(`Cell Line`)

cell_factors$`Cell Line` <- factor(cell_factors$`Cell Line`, levels = cell_order)

ggplot(cell_factors, aes(x=Dose, y=Viability, group=`Cell Line`, color=`Cell Line`)) +
  geom_line() +
  theme_bw() +
  theme(axis.text=element_text(face="bold"),
        axis.text.x=element_text(hjust=0.5),
        legend.text=element_text(size=8, face="bold")) +
  ggtitle("Average Cell Line Viability") +
  xlab("Concentration (nM)") +
  ylab("Mean cell viability (scaled)")

ggsave(snakemake@output[[3]], width=1920, height=1080, units="px", dpi=192)

# 4) average cell viability plot (by cell line cluster)
cell_factors$Cluster <- cell_clusters$cluster[match(cell_factors$`Cell Line`, cell_clusters$cell)]

ggplot(cell_factors, aes(x=Dose, y=Viability, group=`Cell Line`, color=`Cluster`)) +
  geom_line() +
  theme_bw() +
  theme(axis.text=element_text(face="bold"),
        axis.text.x=element_text(hjust=0.5),
        legend.text=element_text(size=8, face="bold")) +
  scale_color_manual(values=pal_cells) +
  ggtitle("Average Cell Line Viability") +
  xlab("Concentration (nM)") +
  ylab("Mean cell viability (scaled)")

ggsave(snakemake@output[[4]], width=1920, height=1080, units="px", dpi=192)
