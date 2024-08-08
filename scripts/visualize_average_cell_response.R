#
# Visualizes cell line viability averaged across all drugs
#
library(tidyverse)

drug_curves <- read_tsv(snakemake@input[[1]], show_col_types=FALSE)

data_cols <- sprintf("dose_%d", 0:10)

cell_factors <- drug_curves %>%
  select(`Cell Line`=cell_line, all_of(data_cols)) %>%
  pivot_longer(-`Cell Line`, names_to="Dose", values_to="viability") %>%
  group_by(`Cell Line`, Dose) %>%
  summarize(Viability=median(viability))

cell_factors$Dose <- factor(cell_factors$Dose, levels=paste0("dose_", 0:10), labels=1:11)

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
  theme(axis.text.x=element_text(hjust=0.5)) +
  theme(legend.text=element_text(size=8)) +
  xlab("Dose number") +
  ylab("Median cell viability")

ggsave(snakemake@output[[1]], width=1920, height=1080, units="px", dpi=192)
