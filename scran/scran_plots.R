#!/usr/bin/env Rscript

########## SETTING DIRECTORIES AND LIBRARIES ############### 

# To fix the Matrix package in Seurat objects
Csparse_validate = 'CsparseMatrix_validate'

### Script to create Scran plots for Master's thesis ###

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
suppressMessages(library(tidyr))
suppressMessages(library(reshape2))
suppressMessages(library(prismatic))
suppressMessages(library(ggallin))
suppressMessages(library(scales))
suppressMessages(library(cowplot))
suppressMessages(library(gtools))
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

CCV_path_cluster <- "Projects/scRNAseq/lluisf/03_Results/"
CCV_path_local <- "03_Results/"

# Setting output path

outpath_cluster <- paste0(path_cluster, CCV_path_cluster)
outpath_local <- paste0(path_local, CCV_path_local)

# Creating directory output (if needed) and setting palette

if(!dir.exists(outpath_cluster) & getwd() == path_cluster){dir.create(outpath_cluster, recursive = T)}
if(!dir.exists(outpath_local) & getwd() == path_local){dir.create(outpath_local, recursive = T)}

##############################################################
## Boxplot for Scran total, tech and bio CCV for all genes ## 
#############################################################

cell_types <- c(
  "CD4_Naive", "CD4_TCM", "CD8_TEM", "NK", "B_naive", "CD8_Naive",
  "CD14_Mono", "B_intermediate", "Treg", "CD4_TEM", "B_memory", "CD16_Mono",
  "MAIT", "CD4_CTL", "CD8_TCM", "NK_CD56bright", "gdT", "cDC2",
  "Plasmablast", "NK_Proliferating", "dnT", "pDC", "Platelet", "HSPC",
  "Eryth", "CD4_Proliferating", "ILC", "CD8_Proliferating", "ASDC",
  "cDC1", "Doublet")

df_celltypes <- lapply(cell_types, function(cell_type) {
  
  if(file.exists(path_cluster)){
    onek1k_age <- paste0("/home/lfronter/Documents/No_phenotype/OneK1K//cell_type/", cell_type, "/", cell_type, ".variance_modelling.subsetted.rds")
    if (file.exists(path_onek1k_age)) {
      onek1k_age <- readRDS(path_onek1k_age)
    } else {
      return(NULL)
    }
  }else if(file.exists(path_local)){
    path_onek1k_age <- paste0("/home/lfronter/Documents/No_phenotype/OneK1K//cell_type/", cell_type, "/", cell_type, ".variance_modelling.subsetted.rds")
    if (file.exists(path_onek1k_age)) {
      onek1k_age <- readRDS(path_onek1k_age)
    } else {
      return(NULL)
    }
  }
  
  onek1k_age <- onek1k_age$dec_all # Get all genes and values
  onek1k_age <- onek1k_age %>%
    mutate(significance = case_when(FDR <= 0.05 ~ "ss", FDR > 0.05 ~ "ns")) # Add status of significant or not significant
  print("Adding status column to distinguish between significant or not significant")
  
  CCV_all_onek1k <- onek1k_age[, c("total","tech","bio")]  # Subset informative columns
  CCV_all_onek1k$celltype <- cell_type
  return(CCV_all_onek1k)
  
})

CCV_all_onek1k <- do.call(rbind, df_celltypes)
CCV_all_onek1k_long <- reshape2::melt(CCV_all_onek1k) # Convert to long format

# Without printing outliers and log scale in y axis #
CCV_all_onek1k_boxplot <- ggplot(CCV_all_onek1k_long, aes(x = variable, y = value)) +
  geom_boxplot(width = 0.5,
               position = position_dodge(width = 2),
               outlier.shape = 21,
               outlier.size = 2,
               aes(fill = variable),
               linewidth = 1) +
  facet_grid(~ fct_inorder(celltype)) +
  scale_fill_discrete(name = "Variability", labels = c("Total", "Technical", "Biological")) +
  theme_bw() +
  scale_color_brewer(palette = "Set2") +
  scale_y_log10(labels = label_scientific()) +
  ylab("CCV (log10)") +
  theme(plot.title = element_text(size = 15),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 22, face = "bold"),
        text = element_text(size = 20),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        strip.text = element_text(size = 12, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Saving plot
boxplot.fn <- paste0("/home/lfronter/Desktop/Figures_TFM/Fig_2G.png")
CairoPNG(boxplot.fn, width = 25, height = 13, units = "in", res = 300)
print(CCV_all_onek1k_boxplot)
dev.off()
print(paste("Boxplot of CCV distribution for OneK1K for all genes using age saved in", boxplot.fn))

###############################################
## Barplot for DVGs number analysis on Scran ##
###############################################

cell_types <- c(
  "CD4_Naive", "CD4_TCM", "CD8_TEM", "NK", "B_naive", "CD8_Naive",
  "CD14_Mono", "B_intermediate", "Treg", "CD4_TEM", "B_memory", "CD16_Mono",
  "MAIT", "CD4_CTL", "CD8_TCM", "NK_CD56bright", "gdT", "cDC2",
  "Plasmablast", "NK_Proliferating", "dnT", "pDC", "Platelet", "HSPC",
  "Eryth", "CD4_Proliferating", "ILC", "CD8_Proliferating", "ASDC",
  "cDC1", "Doublet")

DVGs_Age <- readRDS(paste0(outpath_local, "DVA/DVA_merged_DVG_cell_type_Age_cat.rds"))

DVGs_Age_summary <- DVGs_Age %>%
  group_by(celltype) %>%
  summarise(number_DVG=n(), .groups="drop") %>% # Get number of genes for each cell type and direction
  mutate(celltype = factor(celltype, levels = cell_types), phenotype = "Age") %>% # Order cell types
  arrange(celltype) 

DVGs_Age_summary_melted <- melt(DVGs_Age_summary)

DVGs_Sex <- readRDS(paste0(outpath_local, "DVA/DVA_merged_DVG_cell_type_Sex.rds"))

DVGs_Sex_summary <- DVGs_Sex %>%
  group_by(celltype) %>%
  summarise(number_DVG=n(), .groups="drop") %>% # Get number of genes for each cell type and direction
  mutate(celltype = factor(celltype, levels = cell_types), phenotype = "Sex") %>% # Order cell types
  arrange(celltype) # Drop the grouping information 

DVGs_Sex_summary_melted <- melt(DVGs_Sex_summary)

DVGs_all_summary_melted <- rbind(DVGs_Age_summary_melted, DVGs_Sex_summary_melted)

new_row <- data.frame(
  celltype = "CD4_Naive",
  phenotype = "Sex",
  variable = "number_DVG",
  value = 5
)

DVGs_all_summary_melted <- rbind(DVGs_all_summary_melted, new_row)

DVGs_all_summary_melted <- DVGs_all_summary_melted[!DVGs_all_summary_melted$celltype %in% c("NK_Proliferating","HSPC"), ] # Remove minor cell types

colors_vc <- colorRampPalette(brewer.pal(18, "Paired"))(length(unique(DVGs_all_summary_melted$celltype))) # Get one color for each cell type
colors_vc_w <- colorspace::lighten(colors_vc, 0.3) # Lighter colors

barplot_all <- ggplot(DVGs_all_summary_melted, aes(x = fct_rev(fct_inorder(celltype)), y = value, fill = as.factor(celltype))) +
  geom_bar(stat = "identity") +
  labs(x = NULL, y = "Number of DVGs") +
  scale_fill_manual(values=colors_vc_w, labels="")+
  facet_wrap(~ fct_inorder(phenotype), nrow = 1, ncol = 2) +
  geom_text(aes(x=celltype, y=value, label=value, 
                vjust = 0.5, hjust = 0), size=8) + 
  coord_flip() +
  theme_bw() +
  theme(legend.position="none",
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title.x=element_text(family = "serif", size = 25),
        axis.title.y=element_text(),
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x=element_text(size = 20),
        strip.text = element_text(size = 25, face = "bold"))

# Barplot for cell type abundance

if(paste0(getwd(), "/")  == path_local){
  meta <- readRDS(paste0(path_local, "/02_Data/metadata_processed.rds"))
} else {
  meta <- readRDS("/gpfs/projects/bsc83/Data/scRNAseq/Yazar2022/metadata_processed.rds")
}

meta$cell_type <- sapply(strsplit(meta[["cell_type"]], " "), function(x) paste(x, collapse = "_")) # Get proper cell type naming
cell_abundance <- meta %>% group_by(cell_type) %>% summarise(ncells = n()) %>%
  filter(cell_type %in% DVGs_all_summary_melted$celltype) %>% arrange(desc(ncells))

colnames(cell_abundance) <- c("celltype", "ncells")

celltype_abundance_plot <- ggplot(cell_abundance, aes(x = fct_rev(fct_inorder(celltype)), y = ncells)) +
  geom_bar(stat = "identity", fill = colors_vc_w) +
  ylab("Cell abundance") +
  scale_y_continuous(labels = unit_format(unit = NULL, scale = 1e-3)) +
  theme_bw() + 
  coord_flip() +
  theme(axis.text.y = element_text(size = 25),
        axis.text.x = element_text(angle=0, hjust=0.5, vjust=0.5, size = 25),
        axis.title.x = element_text(family = "serif", size = 23),
        axis.title.y = element_blank(),
        panel.border = element_rect(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

############### Save plot #################

plot_merged <- plot_grid(celltype_abundance_plot, barplot_all, align = "h", axis = "tb", nrow = 1, ncol = 2, rel_widths = c(0.4, 1)) # Merge both plots

# Saving plot
barplot.fn <- paste0("/home/lfronter/Desktop/Figures_TFM/Fig_2C.png")
CairoPNG(barplot.fn, width = 18, height = 12, units = "in", res = 300)
print(plot_merged)
dev.off()

###################################################
## Permutation boxplots for Age and Sex on Scran ##
###################################################

################### FUNCTIONS ###########################

get_DVGs_perm <- function(cell_level, cell_type, pheno, perm_tag, perm_index_tag) {
  
  if(pheno=="Sex"){
    covariates_path <- "/Age.assignment.date/"
  } else {
    covariates_path <- "/Sex.assignment.date/"
  }
  if(getwd() == path_cluster){
    file <- paste0(path_cluster, CCV_path_cluster, "OneK1K/Permutations/", perm_tag, "/", pheno, covariates_path, cell_level, "/", cell_type, "/", cell_type, ".CCV_diff.", perm_index_tag, ".rds")
  } else {
    file <- paste0(path_local, CCV_path_local, "OneK1K/Permutations/", perm_tag, "/", pheno, covariates_path, cell_level, "/", cell_type, "/", cell_type, ".CCV_diff.", perm_index_tag, ".rds")
  }
  if(file.exists(file)){
    all_genes <- readRDS(file)$dec_all
    DVGs <- subset(all_genes, z_score.fdr < 0.05) # Get all genes in which z_score.fdr is below 0.05
    print("Only genes with a Z_score FDR below 0.05 have been selected")
    DVGs <- DVGs %>%
      mutate(regulation = case_when(CCV_diff <= 0 ~ "Down", CCV_diff > 0 ~ "Up")) # Add status of down- or upregulated
    print("Adding status column to distinguish between downregulated or upregulated")
    DVGs$phenotype <- pheno # Add phenotype
    DVGs$status <- perm_index_tag # Add permutation index
    return(DVGs)
  }
}

# Get real DVGs

get_DVGs_real <- function(cell_level, cell_type, pheno){
  
  if(pheno=="Sex"){
    covariates_path <- "/Age.assignment.date/"
  } else {
    covariates_path <- "/Sex.assignment.date/"
  }
  if(getwd() == path_cluster){
    file <- paste0(path_cluster, CCV_path_cluster, "OneK1K/", pheno, covariates_path, cell_level, "/", cell_type, "/", cell_type, ".CCV_diff.rds")
  } else {
    file <- paste0(path_local, CCV_path_local, "OneK1K/", pheno, covariates_path, cell_level, "/", cell_type, "/", cell_type, ".CCV_diff.rds")
  }
  if(file.exists(file)){
    all_genes <- readRDS(file)$dec_all
    DVGs <- subset(all_genes, z_score.fdr < 0.05) # Get all genes in which z_score.fdr is below 0.05
    print("Only genes with a Z_score FDR below 0.05 have been selected")
    DVGs <- DVGs %>%
      mutate(regulation = case_when(CCV_diff <= 0 ~ "Down", CCV_diff > 0 ~ "Up")) # Add status of down- or upregulated
    print("Adding status column to distinguish between downregulated or upregulated")
    DVGs$phenotype <- pheno # Add phenotype
    DVGs$status <- "real" # Add permutation index
    return(DVGs)
  }
}

#################################### ANALYSIS #############################################

vec_perm_indexes <- 1:20 # Get vector of permutation indexes
vec_perm_indexes_tag <- paste0("p_", vec_perm_indexes) # Get vector of permutation indexes
perm_tag <- "p_20" # Get permutation tag

################ FOR AGE ###############

# Apply function to each permutation
temp_df_DVGs_perm_age <- lapply(vec_perm_indexes_tag, function(perm_index_tag) get_DVGs_perm("cell_type", "CD4_Naive", "Age_cat_all", perm_tag, perm_index_tag)) 
DVGs_perm_age <- do.call(rbind.data.frame, temp_df_DVGs_perm_age)

# Apply function to real DVGs
temp_df_DVGs_real_age <- get_DVGs_real("cell_type", "CD4_Naive", "Age_cat_all") 
DVGs_real_age <- temp_df_DVGs_real_age 

# Summarize results and get merged dataframe
DVGs_summary_age <- DVGs_perm_age %>% group_by(status) %>% summarize(num_DVGs = n()) # Get sum of DVGs for each permutation
DVGs_summary_age$total <- nrow(DVGs_real_age) # Get total DVGs counts across permutations
DVGs_summary_age <- DVGs_summary_age %>% mutate(ratio = num_DVGs / total)  # Calculate ratio between DVG count vs total for each permutation
DVGs_summary_age$status <- mixedsort(DVGs_summary_age$status) # Sort by increasing order
DVGs_summary_age$celltype <- "CD4_Naive" # Add cell type
DVGs_summary_age$phenotype <- "Age" # Add phenotype

################ FOR SEX ###############

# Apply function to each permutation
temp_df_DVGs_perm_sex <- lapply(vec_perm_indexes_tag, function(perm_index_tag) get_DVGs_perm("cell_type", "CD4_Naive", "Sex", perm_tag, perm_index_tag)) 
DVGs_perm_sex <- do.call(rbind.data.frame, temp_df_DVGs_perm_sex)

# Apply function to real DVGs
temp_df_DVGs_real_sex <- get_DVGs_real("cell_type", "CD4_Naive", "Sex") 
DVGs_real_sex <- temp_df_DVGs_real_sex

# Summarize results and get merged dataframe
DVGs_summary_sex <- DVGs_perm_sex %>% group_by(status) %>% summarize(num_DVGs = n()) # Get sum of DVGs for each permutation
DVGs_summary_sex$total <- nrow(DVGs_real_sex) # Get total DVGs counts across permutations
DVGs_summary_sex <- DVGs_summary_sex %>% mutate(ratio = num_DVGs / total)  # Calculate ratio between DVG count vs total for each permutation
DVGs_summary_sex$status <- mixedsort(DVGs_summary_sex$status) # Sort by increasing order
DVGs_summary_sex$celltype <- "CD4_Naive" # Add cell type
DVGs_summary_sex$phenotype <- "Sex" # Add phenotype

DVGs_summary_all <- rbind(DVGs_summary_age, DVGs_summary_sex)

# Boxplot for number of DVG in each permutation vs real number of DVG with correct labels

boxplot_permutations_all <- ggplot(DVGs_summary_all, aes(x = celltype, y = num_DVGs)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90") +
  facet_wrap(~ fct_inorder(phenotype), nrow = 2, ncol = 1, scales = "free_x") +
  geom_point(data = DVGs_summary_all, aes(x = celltype, y = num_DVGs, color = "Permutation"), 
             size = 3, shape = 16, position = position_jitter(width = 0.1)) +
  geom_jitter(aes(y = total, color = "Main analysis"), width = 0, height = 0, size = 5) +
  scale_color_manual(values = c("Permutation" = "coral", "Main analysis" = "#00B2B2"), guide = "legend") +
  scale_y_continuous(breaks = seq(0, ceiling(max(DVGs_summary_all$num_DVGs)), by = 2)) +
  ylab("Number of DVGs") +
  xlab(NULL) +
  coord_flip() +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 13),
        axis.text.x = element_text(size = 11), 
        axis.title = element_text(size = 12),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(color = "black", fill = "grey90", size = 1), 
        strip.text.x = element_text(size = 14, face = "bold", color = "black"),
        strip.placement = "outside")  +
  theme(legend.key.size = unit(1, "cm"),
        legend.position = "bottom") + 
  guides(color = guide_legend(
    override.aes = list(
      size = c(5, 3),  
      shape = c(16, 16) 
    )
  ))  

############### Save plot #################

# Saving plot
boxplot.fn <- paste0("/home/lfronter/Desktop/Figures_TFM/Fig_2F.png")
CairoPNG(boxplot.fn, width = 8, height = 10, units = "in", res = 300)
print(boxplot_permutations_all)
dev.off()