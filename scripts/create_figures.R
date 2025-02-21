#
# Creates manuscript figures
#
library(tidyverse)
library(ggrepel)

set.seed(1)

snek <- snakemake

# colorblind friend palettes
# https://github.com/JLSteenwyk/ggpubfigs
pal_cells <- c("#648FFF", "#FE6100", "#785EF0")

# 3a) average cell viability curves
drug_curves <- read_tsv(snek@input[[1]], show_col_types=FALSE)

data_cols <- sprintf("dose_%d", 0:10)

cell_factors <- drug_curves %>%
  select(`Cell Line`=cell_line, all_of(data_cols)) %>%
  pivot_longer(-`Cell Line`, names_to="Dose", values_to="viability") %>%
  group_by(`Cell Line`, Dose) %>%
  summarize(Viability=mean(viability))

# most common concentrations for each dose (nM)
drug_conc <- c(0.779, 2.34, 7.01, 21.0, 63.1, 189.0, 568, 1700, 5116, 15300, 46000)

# load cell average viabilities
cell_factors$Dose <- factor(cell_factors$Dose, levels=data_cols, labels=as.character(drug_conc))

# assign colors based on average viability
cell_order <- cell_factors %>%
  group_by(`Cell Line`) %>%
  summarize(mean_viability=mean(Viability)) %>%
  arrange(-mean_viability) %>%
  pull(`Cell Line`)

cell_factors$`Cell Line` <- factor(cell_factors$`Cell Line`, levels=cell_order)

ggplot(cell_factors, aes(x=Dose, y=Viability, group=`Cell Line`, color=`Cell Line`)) +
  geom_line() +
  theme_bw() +
  theme(axis.text=element_text(face="bold"),
        axis.text.x=element_text(hjust=0.5),
        legend.text=element_text(size=8, face="bold")) +
  ggtitle("Average Cell Line Viability") +
  xlab("Concentration (nM)") +
  ylab("Mean cell viability (scaled)")

ggsave(snek@output[[1]], width=3000, height=1690, units="px", dpi=300, device="tiff")

#
# 3b) average drug ac-50 histogram
#
ac50_mat <- read_tsv(snek@input[[2]])

# compute median ac-50 values for each drug
df <- data.frame(x=apply(ac50_mat[, -1], 2, median, na.rm=TRUE))

ggplot(df, aes(x)) +
  geom_histogram() +
  xlab("AC-50 (median)") +
  ylab("Count") +
  ggtitle("Distribution of Drug Median AC-50 Values")

ggsave(snek@output[[2]], width=3000, height=1700, units="px", dpi=300, device="tiff")
