#!/usr/bin/env Rscript

### Script to merge Scran RDS output files from different cell types (but same phenotype) into one single file ###
                                          
## Merging can be made at two levels: ##
#  1. Based on phenotype used #
#  2. Based on merging only DVGs (using significance status) or all genes (regardless of the significance status) #

# Activate required libraries
suppressMessages(library(dplyr))
suppressMessages(library(purrr))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))

# Setting working directory
path_cluster <- '/gpfs/projects/bsc83/'
path_local <- "/home/lfronter/Documents/"

# Setting path to stored CCV files
CCV_path_cov <- "Projects/scRNAseq/lluisf/03_Results/OneK1K/"

# Setting path to cell types list
path_celltypes <- paste0(path_cluster, 'Projects/scRNAseq/lluisf/01_Scripts/order_cells.rds')

# Setting output path
outpath_cov <- paste0(path_cluster,'Projects/scRNAseq/lluisf/03_Results/DVA/')

################ FUNCTIONS ################## 

# Creating function to merge all celltypes for each phenotype and store them into a RDS file, selecting the DVGs based on z_score.fdr and CCV_diff

get_DVG_cov <- function(cell_level, cell_type, pheno){
  if(pheno=="Sex"){
    covariates_path <- "/Age.assignment.date/"
  } else {
    covariates_path <- "/Sex.assignment.date/"
  }
  file <- paste0(path_local, CCV_path_cov, pheno, covariates_path, cell_level, "/", cell_type, "/", cell_type, ".CCV_diff.rds")
  print(file)
  if(file.exists(file)){
    reanalysis_S <- readRDS(file)$dec_all
    reanalysis_S$celltype <- cell_type
    reanalysis_S$phenotype <- pheno
    reanalysis <- subset(reanalysis_S, z_score.fdr < 0.05) # Get all genes in which z_score.fdr is below 0.05
    print("Only genes with a Z_score FDR below 0.05 have been selected")
    reanalysis <- reanalysis %>%
      mutate(regulation = case_when(CCV_diff <= 0 ~ "Down", CCV_diff > 0 ~ "Up")) # Add status of down- or upregulated
    print("Adding status column to distinguish between downregulated or upregulated")
    return(reanalysis)
  }
}

# Creating function to merge all celltypes for each phenotype and store them into a RDS file, selecting all genes and flagging them based on z_score.fdr and CCV_diff

get_all_cov <- function(cell_level, cell_type, pheno){
  if(pheno=="Sex"){
    covariates_path <- "/Age.assignment.date/"
  } else {
    covariates_path <- "/Sex.assignment.date/"
  }
  file <- paste0(path_cluster, CCV_path_cov, pheno, covariates_path, cell_level, "/", cell_type, "/", cell_type, ".CCV_diff.rds")
  if(file.exists(file)){
    reanalysis_S <- readRDS(file)$dec_all
    reanalysis_S$celltype <- cell_type
    reanalysis_S$phenotype <- pheno
    reanalysis_S <- reanalysis_S %>%
      mutate(regulation = case_when(z_score.fdr < 0.05 & CCV_diff < 0 ~ "Down",
                               z_score.fdr < 0.05 & CCV_diff > 0 ~ "Up",
                               z_score.fdr >= 0.05 ~ "No_DVG"
                               ))
    print("Adding status column to distinguish between downregulated, upregulated or no differentially expressed")
    return(reanalysis_S)
  }
}

################ ANALYSIS ################## 

order_cells <- readRDS(path_celltypes)
cells <- order_cells[["cell_type"]]

## WITH COVARIATES (AGE_CAT) ## 

# For only DVG
temp_df_DVG_cov_age_cat <- lapply(cells, function(cell_type) get_DVG_cov("cell_type", cell_type, "Age_cat"))
temp_df_DVG_cov_age_cat <- temp_df_DVG_cov_age_cat %>% discard(is.null)
result_DVG_cov_age_cat <- do.call(rbind.data.frame, temp_df_DVG_cov_age_cat)

# For all genes
temp_df_all_cov_age_cat <- lapply(cells, function(cell_type) get_all_cov("cell_type", cell_type, "Age_cat"))
temp_df_all_cov_age_cat <- temp_df_all_cov_age_cat %>% discard(is.null)
result_all_cov_age_cat <- do.call(rbind.data.frame, temp_df_all_cov_age_cat)

# Saving results
if(!dir.exists(outpath_cov)) {
 dir.create(outpath_cov, recursive = T)
}
outfile_DVG_cov <- paste0(outpath_cov, "DVA_merged_DVG_cell_type_Age_cat.rds")
saveRDS(result_DVG_cov_age_cat, outfile_DVG_cov)
print(paste("DVG for Age_cat using cell_type saved in", outfile_DVG_cov))
outfile_all_cov <- paste0(outpath_cov, "DVA_merged_all_cell_type_Age_cat.rds")
saveRDS(result_all_cov_age_cat, outfile_all_cov)
print(paste("All genes for Age_cat using cell_type saved in", outfile_all_cov))

## WITH COVARIATES (SEX) ## 

# For only DVG
temp_df_DVG_cov_sex <- lapply(cells, function(cell_type) get_DVG_cov("cell_type", cell_type, "Sex"))
temp_df_DVG_cov_sex <- temp_df_DVG_cov_sex %>% discard(is.null)
result_DVG_cov_sex <- do.call(rbind.data.frame, temp_df_DVG_cov_sex)

# For all genes
temp_df_all_cov_sex <- lapply(cells, function(cell_type) get_all_cov("cell_type", cell_type, "Sex"))
temp_df_all_cov_sex <- temp_df_all_cov_sex %>% discard(is.null)
result_all_cov_sex <- do.call(rbind.data.frame, temp_df_all_cov_sex)

# Saving results
if(!dir.exists(outpath_cov)) {
  dir.create(outpath_cov, recursive = T)
}
outfile_DVG_cov <- paste0(outpath_cov, "DVA_merged_DVG_cell_type_Sex.rds")
saveRDS(result_DVG_cov_sex, outfile_DVG_cov)
print(paste("DVG for Sex using cell_type saved in", outfile_DVG_cov))
outfile_all_cov <- paste0(outpath_cov, "DVA_merged_all_cell_type_Sex.rds")
saveRDS(result_all_cov_sex, outfile_all_cov)
print(paste("All genes for Sex using cell_type saved in", outfile_DVG_cov))
