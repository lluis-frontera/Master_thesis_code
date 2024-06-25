#!/usr/bin/env Rscript

# To fix the Matrix package in Seurat objects
Csparse_validate = 'CsparseMatrix_validate'

### Script to create plots based on merged RDS files using DVG from different immune cell types in OneK1K database ###

# Activate required libraries

suppressMessages(library(dplyr))
suppressMessages(library(purrr))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(UpSetR))
suppressMessages(library(Seurat))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(forcats))
suppressMessages(library(grid))
suppressMessages(library(ggpubr))
suppressMessages(library(cowplot))
suppressMessages(library(Cairo))

# Setting working directory (cluster or local)

path_cluster <- "/gpfs/projects/bsc83/"
path_local <- "/home/lfronter/Documents/"

if(file.exists(path_cluster)){
  setwd(paste(path_cluster))
}else if(file.exists(path_local)){
  setwd(paste(path_local))
}

# Setting path to stored CCV files (same as the output path)

CCV_path_cluster <- "Projects/scRNAseq/lluisf/03_Results/DVA/"
CCV_path_local <- "03_Results/DVA/"

# Setting output path

outpath_cluster <- paste0(path_cluster, CCV_path_cluster)
outpath_local <- paste0(path_local, CCV_path_local)

########################### Plots using 35, 40 and 45 y.o. age thresholds using Age_cat_all #######################

# Extracting DVGs from dataframes using different age thresholds
DVGs_thresholds <- lapply(c(35, 40, 45), function(i){ 
  file <- paste0("/home/lfronter/Documents/03_Results/OneK1K/Age_cat_all/Sex.assignment.date/threshold", i, "/cell_type/CD4_Naive/CD4_Naive.CCV_diff.rds")
  if(file.exists(file)){
    reanalysis_S <- readRDS(file)$dec_all
    reanalysis_S$age_threshold <- i
    reanalysis_S <- reanalysis_S %>%
      filter(z_score.fdr < 0.05)
    print(paste("There are", nrow(reanalysis_S), "DVGs for age threshold", i))
    return(as.data.frame(reanalysis_S))
  }
})

DVGs_thresholds <- DVGs_thresholds %>% discard(is.null)
DVGs_thresholds <- do.call(rbind.data.frame, DVGs_thresholds)

########################### Plots using extremes of the data (<40yo and >60yo) using Age_cat #######################

# Extracting DVGs from dataframes using different age thresholds
filename <- paste0("/home/lfronter/Documents/03_Results/OneK1K/Age_cat/Sex.assignment.date/threshold40/cell_type/CD4_Naive/CD4_Naive.CCV_diff.rds")
DVGs_extremes <- readRDS(filename)$dec_all
DVGs_extremes$age_threshold <- "Extremes (<40 & >60)"
DVGs_extremes <- DVGs_extremes %>%
  filter(z_score.fdr < 0.05) 
print(paste("There are", nrow(DVGs_extremes), "DVGs for age threshold 40"))

########################## Merging and plotting #################################

# Merging

results_df <- rbind(DVGs_thresholds, DVGs_extremes)

# Plotting results
colors_vc <- colorRampPalette(brewer.pal(4, "Set1"))(length(unique(results_df$age_threshold))) # Get one color for each cell type

barplot_all <- ggplot(results_df, aes(x = fct_inorder(as.factor(age_threshold)), fill = as.factor(age_threshold))) +
  geom_bar(stat = "count") +
  labs(x = "Age threshold", y = "Number of DVGs") +
  scale_fill_manual(values=colors_vc, labels="") +
  theme_bw() +
  theme(legend.position="none",
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x=element_text(family = "serif", size = 12),
        axis.title.y=element_text(),
        axis.ticks.y = element_blank(),
        axis.text.x=element_text(), 
        plot.title = element_text(hjust = 0.5))

# Saving plot
barplot.fn <- paste0("/home/lfronter/Documents/Barplot_age_thresholds.png")
CairoPNG(barplot.fn, width = 10, height = 8, units = "in", res = 300)
print(barplot_all)
dev.off()