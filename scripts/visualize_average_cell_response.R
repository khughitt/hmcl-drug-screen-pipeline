#
# Visualizes cell line viability averaged across all drugs
#
library(tidyverse)

cell_factors <- read_tsv(snakemake@input[[1]], show_col_types=FALSE)

cell_factors$Dose <- factor(cell_factors$Dose, levels=1:11, labels=paste0("dose_", 0:10))

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
  ylab("Mean cell viability")

ggsave(snakemake@output[[1]], width=1920, height=1080, units="px", dpi=192)
