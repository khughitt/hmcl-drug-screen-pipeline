#
# Creates a modified version of the Keats cell line metadata excluding cells not included in the
# experiment and selecting a subset of fields of interest and creating some new derived fields from
# them.
#
# Original metadata can be found at: https://www.keatslab.org/myeloma-cell-lines/hmcl-characteristics
#
library(tidyverse)

# 47 cell lines
cells <- c("AMO1_DSMZ", "ARD_JJKsccE7", "ARP1_JJKsccF8", "Delta47_JCRB", "EJM_DSMZ", "FR4_PLB",
           "H1112_PLB", "INA6_PLB", "JIM1_ECACC", "JIM3_ECACC", "JJN3_DSMZ", "Karpas25_ECACC",
           "Karpas417_ECACC", "Karpas620_DSMZ", "KHM11_PLB", "KMM1_JCRB", "KMS11_JCRBsus",
           "KMS12BM_JCRB", "KMS12PE_JCRB", "KMS20_JCRB", "KMS21BM_JCRB", "KMS26_JCRB", "KMS27_JCRB",
           "KMS28BM_JCRB", "KMS28PE_JCRB", "KMS34_JCRB", "L363_DSMZ", "LP1_DSMZ", "MM1R_ATCC",
           "MM1S_ATCC", "MMM1_PLB", "MOLP8_DSMZ", "NCIH929_DSMZ", "OCIMY1_PLB", "OCIMY5_PLB",
           "OCIMY7_PLB", "OPM1_PLB", "OPM2_DSMZ", "PCM6_RIKEN", "PE2_PLB", "RPMI8226_ATCC",
           "SKMM1_PLB", "U266_ATCC", "UTMC2_PLB", "VP6_DJ", "XG1_PLB", "XG6_PLB")

df <- read_tsv(snakemake@input[[1]], show_col_types=FALSE)

df <- df[df$Keats_Lab_Name %in% cells, ]

# add sarin et al., (2020) rankings
sarin2020 <- read_tsv(snakemake@input[[2]], show_col_types=FALSE) %>%
  select(cell_line=`Cell Line`, rank=Rank)

mapping <- data.frame(
  keats = c("AMO-1", "ARD", "ARP1", "Delta47", "EJM", "FR4", "H1112", "INA-6", "JIM-1", "JIM-3",
            "JJN-3", "Karpas 25", "Karpas 417", "Karpas 620", "KHM-11", "KMM-1", "KMS-11",
            "KMS-12BM", "KMS-12PE", "KMS-20", "KMS-21BM", "KMS-26", "KMS-27", "KMS-28BM",
            "KMS-28PE", "KMS-34", "L363", "LP-1", "MM-M1", "MM.1R", "MM.1S", "MOLP-8", "NCI-H929",
            "OCI-MY1", "OCI-MY5", "OCI-MY7", "OPM-1", "OPM-2", "PCM-6", "PE2", "RPMI-8226",
            "SKMM-1", "U-266", "UTMC-2", "VP-6", "XG-1", "XG-6"),
  sarin = c("AMO1", NA, NA, "Delta47", "EJM", "FR4", "H1112", "INA6", "JIM1", "JIM3", "JJN3",
            "Karpas25", "Karpas417", "Karpas620", "KHM11", "KMM1", "KMS11_Sus", "KMS12BM",
            "KMS12PE", "KMS20", "KMS21BM", "KMS26", "KMS27", "KMS28BM", "KMS28PE", "KMS34",
            "L363", "LP1", "MMM1", "MM1R", "MM1S", "MOLP8", "NCIH929", "OCIMY1", "OCIMY5",
            "OCIMY7", "OPM1", "OPM2", "PCM6", "PE2", "RPMI8226", "SKMM1", "U266", "UTMC2",
            "VP6", "XG1", "XG6"))

sarin_ids <- mapping$sarin[match(df$cell_name, mapping$keats)]
df$sarin2020_rank <- sarin2020$rank[match(sarin_ids, sarin2020$cell_line)]

df <- df %>%
  select(cell=Keats_Lab_Name, cell_name=`Public Name`, sex=Sex, ancestry=Ancestry, 
         clinical_heavy_chain=`Clinical Heavy Chain`, clinical_light_chain=`Clinial Light Chain`,
         culture_additives=`Culture Additives`, canonical_translocations=Canonical_Translocations, 
         kras=KRAS, nras=NRAS, tp53=TP53, traf3=TRAF3, sarin2020_rank)

df$canonical_translocations <- gsub(':', ';', df$canonical_translocations)

# derived fields
df <- df %>%
  mutate(multiple_translocations=str_detect(canonical_translocations, '\\+')) %>%
  mutate(il6=str_detect(culture_additives, 'IL6')) %>%
  mutate(transloc_11_14=str_detect(canonical_translocations, '\\(11;14')) %>%
  mutate(transloc_4_14=str_detect(canonical_translocations, '\\(4;14')) %>%
  mutate(wt_kras=str_detect(kras, 'Wt')) %>%
  mutate(wt_nras=str_detect(nras, 'Wt')) %>%
  mutate(wt_tp53=str_detect(tp53, 'Wt')) %>%
  mutate(wt_traf3=str_detect(traf3, 'Wt')) %>%
  select(-culture_additives)

df %>%
  write_tsv(snakemake@output[[1]])
