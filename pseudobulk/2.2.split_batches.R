#!/usr/bin/env Rscript

########## SETTING DIRECTORIES AND LIBRARIES ############### 

# To fix the Matrix package in Seurat objects
Csparse_validate = 'CsparseMatrix_validate'

### Script to perform DEA on Age and Sex variables by splitting whatever cell type into its batches based on the "date" columns ###
  ## Barplots for Age and Sex distribution, as well as barplots for absolute difference between groups (O/Y, M/Y) are done ## 

# Activate required libraries

shhh <- suppressPackageStartupMessages
shhh(library(SingleCellExperiment))
shhh(library(dplyr))
shhh(library(purrr))
shhh(library(ggplot2))
shhh(library(RColorBrewer))
shhh(library(UpSetR))
shhh(library(Seurat))
shhh(library(forcats))
shhh(library(grid))
shhh(library(tidyr))
shhh(library(ggpubr))
shhh(library(reshape2))
shhh(library(prismatic))
shhh(library(optparse))

############################## OPTIONS PARSER ###################################### 

option_list = list(
  make_option(c("--cell_level"), action="store", default=NA, type='character',
              help="From Azimuth: low resolution (predicted.celltype.l1) or high resolution (cell_type)"),
  make_option(c("--cell_type"), action="store", default=NA, type='character',
              help="Cell types in low resolution (predicted.celltype.l1) or high resolution (cell_type)"))
opt = parse_args(OptionParser(option_list=option_list))

# Setting working directory (cluster or local)

path_cluster <- "/gpfs/projects/bsc83/Projects/scRNAseq/lluisf/"
path_local <- "/home/lfronter/Documents/"

if(file.exists(path_cluster)){
  setwd(paste(path_cluster))
}else if(file.exists(path_local)){
  setwd(paste(path_local))
}

# Setting path to stored CCV files (same as the output path)

CCV_path_cluster <- "03_Results/DVA/"
CCV_path_local <- "03_Results/DVA/"

# Setting output path

outpath_cluster <- paste0(path_cluster, CCV_path_cluster)
outpath_local <- paste0(path_local, CCV_path_local)

# Creating directory output (if needed)

if(!dir.exists(outpath_cluster) & getwd() == path_cluster){dir.create(outpath_cluster, recursive = T)}
if(!dir.exists(outpath_local) & getwd() == path_local){dir.create(outpath_local, recursive = T)}

############ ANALYSIS ############### 

Import cell type raw data and check distribution of Age and Sex across batches

sce <- readRDS(paste0("/gpfs/projects/bsc83/Data/scRNAseq/Yazar2022/sce_data_objects/", opt$cell_type, "_", opt$cell_level, "_sceraw.rds"))
sce <- colData(sce) # Getting metadata
sce$Age_cat_all <- ifelse(sce$Age<=40, 'Y', 'O') # Adding age_cat_all column
sce <- as.data.frame(sce) %>% tidyr::separate(date, into = c("remainder", "date"), sep = "_", extra = "merge") # Separate batch column to get only the ID number
sce <- sce[,!names(sce) %in% "remainder"] # Remove remainder column

##### PLOTTING DISTRIBUTION OF AGE AND SEX FOR CELL TYPE #####

## For Age_cat ##

age_cat_distr <- ggplot(sce, aes(x = as.factor(date), fill = Age_cat_all)) +
  geom_bar(stat = "count", position = position_dodge(width = 0.75), color = "black", width = 0.7) +
  labs(x = NULL, y = "Count") +
  scale_fill_manual(values = c("O" = "#B19CD9", "Y" = "#FFB6C1")) +
  ggtitle("Age distribution by batch") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 6, vjust = 0.5, angle = 90),
        text = element_text(family = "sans"),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        legend.key.size = unit(0.5, "cm")) +
  guides(fill = guide_legend(title = "Age status"))

if(file.exists(path_cluster)){
  outpath <- paste0(outpath_cluster, "/Exploratory_analysis/", opt$cell_type, "/")
  if(!dir.exists(outpath)){dir.create(outpath, recursive = T)}
  outpath <- paste0(outpath, "Barplot_batches_age_cat.pdf")
} else if (file.exists(path_local)){
  outpath <- paste0(outpath_local, "/Exploratory_analysis/", opt$cell_type, "/")
  if(!dir.exists(outpath)){dir.create(outpath, recursive = T)}
  outpath <- paste0(outpath, "Barplot_batches_age_cat.pdf")
}

pdf(file = outpath)
print(age_cat_distr)
dev.off()
print(paste("Barplot for age distribution for each batch in", opt$cell_type, "saved in", outpath))

## For Sex ##

sex_distr <- ggplot(sce, aes(x = as.factor(date), fill = Gender)) +
  geom_bar(stat = "count", position = position_dodge(width = NULL), color = "black", width = 0.7) +
  labs(x = NULL, y = "Count") +
  scale_fill_manual(values = c("M" = "#B19CD9", "F" = "#FFB6C1")) +
  ggtitle("Sex distribution by batch") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 6, vjust = 0.5, angle = 90),
        text = element_text(family = "sans"),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        legend.key.size = unit(0.5, "cm")) +
  guides(fill = guide_legend(title = "Sex"))


if(file.exists(path_cluster)){
  outpath <- paste0(outpath_cluster, "/Exploratory_analysis/", opt$cell_type, "/")
  if(!dir.exists(outpath)){dir.create(outpath, recursive = T)}
  outpath <- paste0(outpath, "Barplot_batches_sex.pdf")
} else if (file.exists(path_local)){
  outpath <- paste0(outpath_local, "/Exploratory_analysis/", opt$cell_type, "/")
  if(!dir.exists(outpath)){dir.create(outpath, recursive = T)}
  outpath <- paste0(outpath, "Barplot_batches_sex.pdf")
}

pdf(file = outpath)
print(sex_distr)
dev.off()
print(paste("Barplot for sex distribution for each batch in", opt$cell_type, "saved in", outpath))

##### PLOTTING ABSOLUTE DIFFERENCE BETWEEEN OLD AND YOUNG GROUPS FOR EACH BATCH #####

## FOR AGE_CAT ##

age_balance_summary <- sce %>% as.data.frame() %>% group_by(date) %>% summarise(age_balance = abs(sum(Age_cat_all == "O") - sum(Age_cat_all == "Y"))) # Calculate absolute difference between old and young groups
age_balance_summary <- age_balance_summary %>% arrange(desc(age_balance)) # Sort by decreasing order

colors_vc <- colorRampPalette(brewer.pal(9, "Pastel1"))(length(age_balance_summary$age_balance)) # Get one color for each batch

abs_diff_age_barplot <- ggplot(age_balance_summary, aes(x = fct_rev(fct_inorder(date)), y = age_balance, fill = date)) +
  geom_bar(stat = "identity", color="black") +
  labs(title = "Old vs young absolute age difference across batches", x = NULL, y = "Abs difference") +
  scale_fill_manual(values = colors_vc) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 7, vjust = 0.5, angle = 90)) +
  guides(fill = "none")

if(file.exists(path_cluster)){
  outpath <- paste0(outpath_cluster, "/Exploratory_analysis/", opt$cell_type, "/")
  if(!dir.exists(outpath)){dir.create(outpath, recursive = T)}
  outpath <- paste0(outpath, "Barplot_balance_diff_batches_age_cat.pdf")
} else if (file.exists(path_local)){
  outpath <- paste0(outpath_local, "/Exploratory_analysis/", opt$cell_type, "/")
  if(!dir.exists(outpath)){dir.create(outpath, recursive = T)}
  outpath <- paste0(outpath, "Barplot_balance_diff_batches_age_cat.pdf")
}

pdf(file = outpath)
print(abs_diff_age_barplot)
dev.off()
print(paste("Barplot for absolute difference between O and Y for each batch using age_cat in", opt$cell_type, "saved in", outpath))

## FOR SEX ##

sex_balance_summary <- sce %>% as.data.frame() %>% group_by(date) %>% summarise(sex_balance = abs(sum(Gender == "M") - sum(Gender == "F"))) # Calculate absolute difference between old and young groups
sex_balance_summary <- sex_balance_summary %>% arrange(desc(sex_balance)) # Sort by decreasing order

colors_vc <- colorRampPalette(brewer.pal(9, "Pastel1"))(length(sex_balance_summary$sex_balance)) # Get one color for each batch

abs_diff_sex_barplot <- ggplot(sex_balance_summary, aes(x = fct_rev(fct_inorder(date)), y = sex_balance, fill = date)) +
  geom_bar(stat = "identity", color="black") +
  labs(title = "Old vs young absolute sex difference across batches", x = NULL, y = "Abs difference") +
  scale_fill_manual(values = colors_vc) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 7, vjust = 0.5, angle = 90)) +
  guides(fill = "none")

if(file.exists(path_cluster)){
  outpath <- paste0(outpath_cluster, "/Exploratory_analysis/", opt$cell_type, "/")
  if(!dir.exists(outpath)){dir.create(outpath, recursive = T)}
  outpath <- paste0(outpath, "Barplot_balance_diff_batches_sex.pdf")
} else if (file.exists(path_local)){
  outpath <- paste0(outpath_local, "/Exploratory_analysis/", opt$cell_type, "/")
  if(!dir.exists(outpath)){dir.create(outpath, recursive = T)}
  outpath <- paste0(outpath, "Barplot_balance_diff_batches_sex.pdf")
}

pdf(file = outpath)
print(abs_diff_sex_barplot)
dev.off()
print(paste("Barplot for absolute difference between O and Y for each batch using sex in", opt$cell_type, "saved in", outpath))

## FOR AGE AND SEX TOGETHER ##

age_sex_balance_summary <- merge(age_balance_summary, sex_balance_summary, by = "date", all = TRUE) %>%
  mutate(sum_values = ifelse(is.na(age_balance), sex_balance, ifelse(is.na(sex_balance), age_balance, age_balance + sex_balance))) %>%
  select(date, sum_values) # Getting sum of absolute difference in age and sex for each batch
age_sex_balance_summary <- age_sex_balance_summary %>% arrange(desc(sum_values)) # Sort by decreasing order

colors_vc <- colorRampPalette(brewer.pal(9, "Pastel1"))(length(age_sex_balance_summary$sum_values)) # Get one color for each batch

abs_diff_all_barplot <- ggplot(age_sex_balance_summary, aes(x = fct_rev(fct_inorder(date)), y = sum_values, fill = date)) +
  geom_bar(stat = "identity", color="black") +
  labs(title = "Old vs young absolute age & sex difference across batches", x = NULL, y = "Abs difference") +
  scale_fill_manual(values = colors_vc) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 7, vjust = 0.5, angle = 90)) +
  guides(fill = "none")

if(file.exists(path_cluster)){
  outpath <- paste0(outpath_cluster, "/Exploratory_analysis/", opt$cell_type, "/")
  if(!dir.exists(outpath)){dir.create(outpath, recursive = T)}
  outpath <- paste0(outpath, "Barplot_balance_diff_batches_all.pdf")
} else if (file.exists(path_local)){
  outpath <- paste0(outpath_local, "/Exploratory_analysis/", opt$cell_type, "/")
  if(!dir.exists(outpath)){dir.create(outpath, recursive = T)}
  outpath <- paste0(outpath, "Barplot_balance_diff_batches_all.pdf")
}

pdf(file = outpath)
print(abs_diff_all_barplot)
dev.off()
print(paste("Barplot for absolute difference between O and Y for each batch using age and sex in", opt$cell_type, "saved in", outpath))

## Retrieving all cells belonging to each cell type and saving them into a RDS file ##

sce <- readRDS(paste0("/gpfs/projects/bsc83/Data/scRNAseq/Yazar2022/sce_data_objects/", opt$cell_type, "_", opt$cell_level, "_sceraw.rds"))

for (pool in unique(sort(sce$date))){
  
  meta_pool <- sce[, sce$date == pool]
  
  if(file.exists(path_cluster)){
    outpath <- paste0(path_cluster, "03_Results/Donor_variability/Splitted_batches/", opt$cell_level, "/", opt$cell_type, "/")
    if(!dir.exists(outpath)){dir.create(outpath, recursive = T)}
    outpath <- paste0(outpath, pool, "_", opt$cell_type, ".rds")
    saveRDS(meta_pool, outpath)
  } else if (file.exists(path_local)){
    outpath <- paste0(path_local, "03_Results/Donor_variability/Splitted_batches/", opt$cell_level, "/", opt$cell_type, "/")
    if(!dir.exists(outpath)){dir.create(outpath, recursive = T)}
    outpath <- paste0(outpath, pool, "_", opt$cell_type, ".rds")
    saveRDS(meta_pool, outpath)
  }
  
  print(paste("Subset for", pool, "using", opt$cell_type, "saved in", outpath))
  
}
