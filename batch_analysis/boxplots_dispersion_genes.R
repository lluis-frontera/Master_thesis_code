#!/usr/bin/env Rscript

# To fix the Matrix package in Seurat objects
Csparse_validate = 'CsparseMatrix_validate'

########## SETTING LIBRARIES ############### 

shhh <- suppressPackageStartupMessages
shhh(library(optparse))
shhh(library(dplyr))
shhh(library(purrr))
shhh(library(ggplot2))
shhh(library(RColorBrewer))
shhh(library(Seurat))
shhh(library(SingleCellExperiment))
shhh(library(forcats))
shhh(library(grid))
shhh(library(ggpubr))
shhh(library(stringi))
shhh(library(tidyr))
shhh(library(cowplot))
shhh(library(ggsignif))
shhh(library(reshape2))
shhh(library(scales))
shhh(library(modeest))
shhh(library(tibble))
shhh(library(Cairo))

options(digits = 4) # Modify number of digits printed

################# PARSER ARGUMENTS #############

option_list = list(
  make_option(c("--cell_type"), action="store", default=NA, type='character',
              help="Cell types in low resolution (predicted.celltype.l1) or high resolution (cell_type)"),
  make_option(c("--cell_level"), action="store", default=NA, type='character',
              help="From Azimuth: low resolution (predicted.celltype.l1) or high resolution (cell_type)"))
opt = parse_args(OptionParser(option_list=option_list))

########## SETTING DIRECTORIES ############### 

path_cluster <- "/gpfs/projects/bsc83/Projects/scRNAseq/lluisf"
path_local <- "/home/lfronter/Documents"

if(startsWith(getwd(), path_local)){
  main.dir <- paste0(path_local, "/03_Results/Donor_variability/cell_type/")
} else {
  main.dir <- paste0(path_cluster, "/03_Results/Donor_variability/cell_type/")
}

###################################################################################################
### Boxplot of dispersion values separating by batches for OneK1K separating by gene expression ###
###################################################################################################

# Read the dispersion and get number of genes tested
dispersion_df <- readRDS(paste0(main.dir, opt$cell_type, "/Dispersion_updated_5threshold.rds"))
genes_tested <- rownames(dispersion_df)

# Read raw count data ansd extract counts for tested genes
celltype_df <- readRDS(paste0("/gpfs/projects/bsc83/Data/scRNAseq/Yazar2022/sce_data_objects/", opt$cell_type, "_", opt$cell_level, "_sceraw.rds"))
counts_df <- counts(celltype_df)
counts_df <- counts_df[genes_tested, , drop = FALSE]

# Calculate expression levels (e.g., sum across samples)
expression_levels <- rowSums(counts_df)

# Divide into quantiles
quantiles <- quantile(expression_levels, probs = c(0, 0.33, 0.66, 1))

# Assign categories based on quantiles
category <- cut(expression_levels, breaks = quantiles, labels = c("Low", "Mid", "High"), include.lowest = TRUE)

# Create a dataframe with gene names and their corresponding categories
gene_categories <- data.frame(gene = rownames(counts_df), category = category)

# Read the meta data
meta <- readRDS(paste0("/gpfs/projects/bsc83/Data/scRNAseq/Yazar2022/metadata_processed.rds"))

# Ensure the genes are in the same order in both dataframes for proper matching
gene_categories_sorted <- gene_categories[match(rownames(dispersion_df), gene_categories$gene), ]

# Set the rownames for gene_categories_sorted to match dispersion_df
rownames(gene_categories_sorted) <- rownames(dispersion_df)

# Create separate dataframes for each category
low_expression_df <- dispersion_df[gene_categories_sorted$category == "Low", , drop = FALSE]
mid_expression_df <- dispersion_df[gene_categories_sorted$category == "Mid", , drop = FALSE]
high_expression_df <- dispersion_df[gene_categories_sorted$category == "High", , drop = FALSE]

# Create list to iterate over it
list_dispersion_df <- list(low_expression_df, mid_expression_df, high_expression_df)

plots_per_batch <- lapply(list_dispersion_df, function(dispersion_df_subset) {
  
  if (nrow(dispersion_df_subset) == 0) {return(NA)} # If dataframe is empty, return NA
    
  # Extract dispersion values for one gene
  gene <- rownames(dispersion_df_subset)[1]
  expression_onegene <- as.data.frame(dispersion_df_subset[gene, ]) # Transpose to make gene a column
  colnames(expression_onegene) <- "value" # Rename the single column to "value"
  expression_onegene$assignment <- rownames(expression_onegene) # Create assignment column from rownames
  
  # Subset and clean the meta data
  meta_clean <- meta %>%
    select(assignment, date) %>%
    unique()
  rownames(meta_clean) <- NULL
  
  # Create a named vector of assignments to dates
  assignment_to_date <- setNames(meta_clean$date, meta_clean$assignment)
  
  # Extract the assignments
  assignments <- expression_onegene$assignment
  
  # Map the assignments to dates
  col_dates <- assignment_to_date[assignments]
  
  # Create a data frame that includes assignments and their corresponding dates
  col_df <- data.frame(assignment = assignments, date = col_dates, stringsAsFactors = FALSE)
  
  # Create a nested dataframe where each batch (date) contains its corresponding dispersion_df columns
  nested_dispersion <- col_df %>%
    group_by(date) %>%
    nest() %>%
    mutate(data = lapply(data, function(x) {
      expression_onegene %>%
        filter(assignment %in% x$assignment) %>%
        select(assignment, value)
    }))
  
  # Unnest the nested_dispersion object
  unnested_dispersion <- nested_dispersion %>%
    unnest(cols = c(data)) %>% 
    filter(!is.na(value)) %>%
    select(date, value)
  
  unnested_dispersion$group <- "batch"
  
  median_batches <- unnested_dispersion %>%
    group_by(date) %>% 
    summarise(median=median(value))
  
  median_batches <- mean(median_batches$median)
  print(median_batches)

  ## Boxplot of dispersion values for whole OneK1K (comparison) 
  
  dispersion_all_onegene <- as.data.frame(dispersion_df_subset[rownames(dispersion_df_subset) == gene, ])
  colnames(dispersion_all_onegene) <- "value" # Rename the single column to "value"
  dispersion_all_onegene$date <- "total" # Create grouping column for x axis
  dispersion_all_onegene$group <- "total" # Create grouping column for coloring
  dispersion_all_onegene <- dispersion_all_onegene %>% filter(!is.na(value)) # Filter NA values
  median_total <- median(dispersion_all_onegene$value)
  print(median_total)

  ## Merge both dataframes and plot 
  
  # Merge both dataframes and get unique means
  result_df <- rbind(unnested_dispersion, dispersion_all_onegene)

  ## Plotting comparison boxplot
  comparison_boxplot <- ggplot(result_df, aes(x=as.factor(date), y=log2(value))) + geom_violin(aes(fill = as.factor(group)), alpha = .75) +
    geom_boxplot(width=0.2, color = "black", coef = 1.5, outlier.alpha = 0, position = position_dodge(width = 0.9)) + scale_fill_brewer(palette = "Dark2") +
    xlab(NULL) + ylab("Dispersion (log2)") + ggtitle(gene) + labs(fill = "Group") + theme_bw() +
    geom_hline(aes(yintercept = log2(0.02043)), linetype = "solid", size = 1, color = "#1B9E77", show.legend = FALSE) +
    geom_hline(aes(yintercept = log2(median_total)), linetype = "solid", size = 1, color = "#D95F02", show.legend = FALSE) + 
    theme(plot.title = element_text(hjust = 0.5, size=30, face = "bold"), axis.text= element_text(size=20, color="black"), legend.title = element_text(size = 20),
          legend.position = "right", legend.text = element_text(size = 15), axis.title = element_text(size=25), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text.x = element_text(angle = -90, hjust = 0.5, vjust = 0.5, size = 12))
  
  return(comparison_boxplot)

})

names(plots_per_batch) <- c("Low", "Mid", "High")

### Saving comparison boxplot
lapply(c("Low", "Mid", "High"), function(i) {
  figfn <- paste0(path_cluster, "/Median_comparison_by_batch_", i, "_expression.png")
  CairoPNG(figfn, width = 18, height = 10, units = "in", res = 300)
  print(plots_per_batch[[i]])
  dev.off()
  print(paste("Boxplot for batch variation saved in", figfn))
})