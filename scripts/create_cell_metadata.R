#
# Creates a modified version of the Keats cell line metadata excluding cells not included in the
# experiment and selecting a subset of fields of interest.

# Original metadata can be found at: https://www.keatslab.org/myeloma-cell-lines/hmcl-characteristics
#
library(tidyverse)

snek <- snakemake

# 47 cell lines
cells <- c("AMO1_DSMZ", "ARD_JJKsccE7", "ARP1_JJKsccF8", "Delta47_JCRB", "EJM_DSMZ", "FR4_PLB",
           "H1112_PLB", "INA6_PLB", "JIM1_ECACC", "JIM3_ECACC", "JJN3_DSMZ", "Karpas25_ECACC",
           "Karpas417_ECACC", "Karpas620_DSMZ", "KHM11_PLB", "KMM1_JCRB", "KMS11_JCRBsus",
           "KMS12BM_JCRB", "KMS12PE_JCRB", "KMS20_JCRB", "KMS21BM_JCRB", "KMS26_JCRB", "KMS27_JCRB",
           "KMS28BM_JCRB", "KMS28PE_JCRB", "KMS34_JCRB", "L363_DSMZ", "LP1_DSMZ", "MM1R_ATCC",
           "MM1S_ATCC", "MMM1_PLB", "MOLP8_DSMZ", "NCIH929_DSMZ", "OCIMY1_PLB", "OCIMY5_PLB",
           "OCIMY7_PLB", "OPM1_PLB", "OPM2_DSMZ", "PCM6_RIKEN", "PE2_PLB", "RPMI8226_ATCC",
           "SKMM1_PLB", "U266_ATCC", "UTMC2_PLB", "VP6_DJ", "XG1_PLB", "XG6_PLB")

# load keats hmcl metadata
df <- read_tsv(snek@input[[1]], show_col_types=FALSE)

df <- df[df$Keats_Lab_Name %in% cells, ]

df <- df %>%
  select(cell=Keats_Lab_Name, cell_name=`Public Name`, sex=Sex, ancestry=Ancestry,
         clinical_heavy_chain=`Clinical Heavy Chain`, clinical_light_chain=`Clinial Light Chain`,
         canonical_translocations=Canonical_Translocations,
         kras=KRAS, nras=NRAS, tp53=TP53, traf3=TRAF3)

df$canonical_translocations <- gsub(":", ";", df$canonical_translocations)

# add average cell line viabilities
# cell_viabilities <- read_tsv(snek@input[[4]], show_col_types=FALSE) %>%
#   group_by(cell=`Cell Line`) %>%
#   summarize(mean_viability=mean(Viability))
#
# df <- df %>%
#   left_join(cell_viabilities, by="cell")

df <- df %>%
  select(cell, cell_name, everything()) %>%
  arrange(cell_name)

# save result
df %>%
  write_tsv(snek@output[[1]])
