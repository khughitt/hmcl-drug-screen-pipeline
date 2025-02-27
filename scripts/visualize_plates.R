#
# visualize_plates.R
#
library(gridExtra)
library(lattice)
library(tidyverse)
library(viridis)

# constants
PLATE_NUM_ROWS <- 32

# column names to use for plate images; each row contains four drugs
cnames <- c("DMSO", "Pos1", "Pos2", "Neg",
            paste0(c("Drug4_dose", "Drug3_dose"), rep(10:0, each=2)),
            paste0(c("Drug2_dose", "Drug1_dose"), rep(10:0, each=2)))

# legend size
legend_params <- list(labels=list(cex=2))

# load plate data
raw_combined <- read_tsv(snakemake@input[[1]], show_col_types=FALSE) %>%
  as.matrix()

normed_combined <- read_tsv(snakemake@input[[2]], show_col_types=FALSE) %>%
  as.matrix()

bgadj_combined <- read_tsv(snakemake@input[[3]], show_col_types=FALSE) %>%
  as.matrix()

background_plate <- read_tsv(snakemake@input[[4]], show_col_types=FALSE) %>%
  as.matrix()

colnames(background_plate) <- cnames

# load metadata
plate_mdata <- read_tsv(snakemake@input[[5]], show_col_types=FALSE)

out_dir <- dirname(snakemake@output$indiv[[1]])

# iterate over plates and create combined figures
for (plate_id in snakemake@params[["plate_ids"]]) {
  # 1d -> 2d
  raw_plate <- matrix(raw_combined[, plate_id], nrow=PLATE_NUM_ROWS)
  normed_plate <- matrix(normed_combined[, plate_id], nrow=PLATE_NUM_ROWS)
  bgadj_plate <- matrix(bgadj_combined[, plate_id], nrow=PLATE_NUM_ROWS)

  colnames(raw_plate) <- cnames
  colnames(normed_plate) <- cnames
  colnames(bgadj_plate) <- cnames

  outfile <- file.path(out_dir, sprintf("%s.jpg", plate_id))

  # include cell line in figure title
  cell_line <- plate_mdata$cell_line[match(plate_id, plate_mdata$plate)]
  plt_title <- sprintf("%s (%s)", plate_id, cell_line)

  jpeg(outfile, width=1920, height=1080)

  # generate plots of plate at each step of processing
  plts <- list()

  plts[[1]] <- levelplot(t(raw_plate), col.regions=plasma(500),
                         scales=list(x=list(rot=90, fontface="bold"),
                                     y=list(fontface="bold")),
                         xlab="", ylab="",
                         colorkey=legend_params,
                         main="Raw")

  plts[[2]] <- levelplot(t(normed_plate), col.regions=plasma(500),
                         scales=list(x=list(rot=90, fontface="bold"),
                                     y=list(fontface="bold")),
                         xlab="", ylab="",
                         colorkey=legend_params,
                         main="Normed")

  plts[[3]] <- levelplot(t(bgadj_plate), col.regions=plasma(500),
                         scales=list(x=list(rot=90, fontface="bold"),
                                     y=list(fontface="bold")),
                         xlab="", ylab="",
                         colorkey=legend_params,
                         main="Background-adjusted")

  grid.arrange(grobs=plts, nrow=2, top=plt_title)
  dev.off()
}

# visualize mean & median of stacked raw plate images
mean_plate <- matrix(apply(raw_combined, 1, mean), nrow=PLATE_NUM_ROWS)

colnames(mean_plate) <- cnames

png(snakemake@output$mean, width=1920, height=1080)
levelplot(t(mean_plate), col.regions=plasma(500),
          scales=list(x=list(rot=90, cex=1.4, fontface="bold"), y=list(cex=1.4, fontface="bold")),
          xlab="", ylab="",
          colorkey=legend_params,
          main=list(label="Raw plates (mean)", cex=2.5))
dev.off()

median_plate <- matrix(apply(raw_combined, 1, median), nrow=PLATE_NUM_ROWS)
colnames(median_plate) <- cnames

png(snakemake@output$median, width=1920, height=1080)
levelplot(t(median_plate), col.regions=plasma(500),
          scales=list(x=list(rot=90, cex=1.4, fontface="bold"), y=list(cex=1.4, fontface="bold")),
          xlab="", ylab="",
          colorkey=legend_params,
          main=list(label="Raw plates (median)", cex=2.5))
dev.off()

# visualize background image
png(snakemake@output$background, width=1920, height=1080)
levelplot(t(background_plate), col.regions=plasma(500),
          scales=list(x=list(rot=90, cex=1.4, fontface="bold"), y=list(cex=1.4, fontface="bold")),
          xlab="", ylab="",
          colorkey=legend_params,
          main=list(label="Background Image", cex=2.5))
dev.off()
