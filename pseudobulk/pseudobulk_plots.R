#!/usr/bin/env Rscript

########## SETTING DIRECTORIES AND LIBRARIES ############### 

# To fix the Matrix package in Seurat objects
Csparse_validate = 'CsparseMatrix_validate'

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


### Script to create pseudobulk plots for Master's thesis ###

path_cluster <- "/gpfs/projects/bsc83/Projects/scRNAseq/lluisf"
path_local <- "/home/lfronter/Documents"

if(startsWith(getwd(), path_local)){
  main.dir <- paste0(path_local, "/03_Results/Donor_variability/cell_type/")
} else {
  main.dir <- paste0(path_cluster, "/03_Results/Donor_variability/cell_type/")
}

#####################################################
## Barplot for DVGs number analysis on pseudobulk ## 
#####################################################

################# Read data ###################

DVGs_df_Age <- readRDS(paste0(main.dir, "DVGs_Age_all_cell_type_5threshold_test.rds"))
DVGs_sig_Age <- DVGs_df_Age %>%
  filter(p_value < 0.05) # Without beta threshold filtering 

DVGs_df_Sex <- readRDS(paste0(main.dir, "DVGs_Sex_all_cell_type_5threshold_test.rds"))
DVGs_sig_Sex <- DVGs_df_Sex %>%
  filter(p_value < 0.05) # Without beta threshold filtering 

################# Functions ###################

Mybinom <- function(subdf){ # Calculate binomial test for each cell type
  n1 <- subdf %>% filter(direction==1) %>% dplyr::pull(ngene) # Get number of downregulated genes 
  n2 <- subdf %>% filter(direction==2) %>% dplyr::pull(ngene) # Get number of upregulated genes 
  if(is.integer(n1) && length(n1) == 0L || is.integer(n2) && length(n2) == 0L) {return(NA)} # If no up- or downregulated genes, add NA
  ngene <- c(n1, n2)
  if(n1>n2){
    res <- binom.test(ngene, 0.5, alternative="greater")
  }else{
    res <- binom.test(ngene, 0.5, alternative="less") 
  }
  res$p.value
}

Mysymb <- function(pval){ # Add significance symbol to each cell type
  if (is.na(pval)) {symb <- NA}
  else if (pval<0.001) {symb <- "***"}
  else if(pval>=0.001 & pval<0.01) {symb <- "**"}
  else if (pval>=0.01 & pval<0.05) {symb <- "*"}
  else if(pval>0.05) {symb <- ""}
  symb
}

Mypos <- function(subdf){ # Get position for each cell type
  ny <- subdf %>% filter(direction==1) %>% dplyr::pull(ngene)
  ny
}

desired_order <- c(
    "CD4_Naive", "CD4_TCM", "CD8_TEM", "NK", "B_naive", "CD8_Naive",
    "CD14_Mono", "B_intermediate", "Treg", "CD4_TEM", "B_memory", "CD16_Mono",
    "MAIT", "CD4_CTL", "CD8_TCM", "NK_CD56bright", "gdT", "cDC2",
    "Plasmablast", "NK_Proliferating", "dnT", "pDC", "Platelet", "HSPC",
    "Eryth", "CD4_Proliferating", "ILC", "CD8_Proliferating", "ASDC",
    "cDC1", "Doublet")

############# Processing for Age ################

sigs_Age <- DVGs_sig_Age %>%
  mutate(direction=ifelse(beta > 0, "1", "2")) %>% # Get direction
  group_by(direction, celltype) %>%
  summarise(ngene=n(), .groups="drop") %>% # Get number of genes for each cell type and direction
  mutate(ngene2=ifelse(direction == 2,-ngene, ngene), # Add negative symbol
         celltype = factor(celltype, levels = desired_order),
         phenotype = "Age") %>% # Order cell types
  arrange(direction, celltype) # Drop the grouping information 

anno_df_Age <- sigs_Age %>% group_by(celltype) %>% nest() %>% # Add symbol and positions
  mutate(pval = map_dbl(data, Mybinom), 
         symb = map_chr(pval, Mysymb),
         ypos = map_dbl(data, Mypos)) 

############# Processing for Sex ################

sigs_Sex <- DVGs_sig_Sex %>%
  mutate(direction=ifelse(beta > 0, "1", "2")) %>% # Get direction
  group_by(direction, celltype) %>%
  summarise(ngene=n(), .groups="drop") %>% # Get number of genes for each cell type and direction
  mutate(ngene2=ifelse(direction == 2,-ngene, ngene), # Add negative symbol
         celltype = factor(celltype, levels = desired_order),
         phenotype = "Sex") %>% # Order cell types
  arrange(direction, celltype) # Drop the grouping information 

anno_df_Sex <- sigs_Sex %>% group_by(celltype) %>% nest() %>% # Add symbol and positions
  mutate(pval = map_dbl(data, Mybinom), 
         symb = map_chr(pval, Mysymb),
         ypos = map_dbl(data, Mypos)) 

############# Merging and plotting ################

sigs_all <- rbind(sigs_Age, sigs_Sex)

## Colors
colors_vc <- colorRampPalette(brewer.pal(18, "Paired"))(length(unique(sigs_Age$celltype))) # Get one color for each cell type
colors_vc_w <- colorspace::lighten(colors_vc, 0.3) # Lighter colors

## Create plot
barplot_all <- ggplot(sigs_all, aes(x=fct_rev(celltype), y=ngene2)) +
  geom_bar(aes(fill=celltype),stat="identity") +
  coord_flip() +
  facet_wrap(~ fct_inorder(phenotype), nrow = 1, ncol = 2, scales = "free_x") +
  ylab("DVGs_number") +
  scale_fill_manual(values=colors_vc_w, labels="")+
  geom_hline(yintercept=0, color="grey60")+
  geom_text(aes(x=celltype, y=ngene2, label=abs(ngene2), 
                vjust=ifelse(direction==2, 0.5, 0.5), hjust=ifelse(direction==2, 1, 0)), size=8) + 
  scale_y_continuous(labels = function(x) abs(x)) +
  theme_bw() +
  theme(legend.position="none",
        axis.title.x=element_text(family = "serif", size = 25),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_text(angle=0, hjust=0.5, vjust=0.5, size = 20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(size = 25, face = "bold"),
        panel.background = element_blank())

# Barplot for cell type abundance

if(getwd()  == path_local){
  meta <- readRDS(paste0(path_local, "/02_Data/metadata_processed.rds"))
} else {
  meta <- readRDS("/gpfs/projects/bsc83/Data/scRNAseq/Yazar2022/metadata_processed.rds")
}

meta$cell_type <- sapply(strsplit(meta[["cell_type"]], " "), function(x) paste(x, collapse = "_")) # Get proper cell type naming
cell_abundance <- meta %>% group_by(cell_type) %>% summarise(ncells = n()) %>%
  filter(cell_type %in% sigs_all$celltype) %>% arrange(desc(ncells))

celltype_abundance_plot <- ggplot(cell_abundance, aes(x = fct_rev(fct_inorder(cell_type)), y = ncells)) +
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

plot_merged <- plot_grid(celltype_abundance_plot, barplot_all, align = "h", axis = "tb", nrow = 1, ncol = 2, rel_widths = c(0.3, 1)) # Merge both plots

# Saving plot
barplot.fn <- paste0("/home/lfronter/Desktop/Figures_TFM/Fig_3C.png")
CairoPNG(barplot.fn, width = 30, height = 16, units = "in", res = 300)
print(plot_merged)
dev.off()

########################################################
## Permutation boxplots for Age and Sex on pseudobulk ##
########################################################

#################################### FUNCTIONS #############################################

get_DVGs_perm <- function(cell_level, cell_type, perm_idx.tag, phenotype) {
  
  file <- paste0(main.dir, perm_idx.tag, "/", cell_level, "/", cell_type, "/Results_", phenotype, "/1.DVG_results_5threshold_test.rds")
  if(file.exists(file)){
    all_genes <- readRDS(file)
    DVGs <- subset(all_genes, p_value < 0.05) # Get all genes in which qval is below 0.05
    print("Only genes with a pval below 0.05 have been selected")
    DVGs <- DVGs %>%
      mutate(regulation = case_when(beta <= 0 ~ "Down", beta > 0 ~ "Up")) # Add status of down- or upregulated
    print("Adding status column to distinguish between downregulated or upregulated")
    DVGs$status <- perm_idx.tag # Add permutation index
    DVGs$celltype <- cell_type
    return(DVGs)
  }
}

get_DVGs_real <- function(cell_level, cell_type, phenotype){
  
  if(startsWith(getwd(), path_local)){
    file <- paste0(path_local, "/03_Results/Donor_variability/cell_type/", cell_type, "/Results_", phenotype, "/1.DVG_results_5threshold_test.rds") 
  } else {
    file <- paste0(path_cluster, "/03_Results/Donor_variability/cell_type/", cell_type, "/Results_", phenotype, "/1.DVG_results_5threshold_test.rds")
  }
  if(file.exists(file)){
    if(file.exists(file)){
      all_genes <- readRDS(file)
      DVGs <- subset(all_genes, p_value < 0.05) # Get all genes in which qval is below 0.05
      print("Only genes with a pval below 0.05 have been selected")
      DVGs <- DVGs %>%
        mutate(regulation = case_when(beta <= 0 ~ "Down", beta > 0 ~ "Up")) # Add status of down- or upregulated
      print("Adding status column to distinguish between downregulated or upregulated")
      DVGs$status <- "real" # Add permutation index
      DVGs$celltype <- cell_type
      return(DVGs)
    }
  }
}

#################################### ANALYSIS #############################################

main.dir <- paste0(path_local, "/03_Results/Donor_variability/Permutations/p_20/")

perm_num.tag <- paste0("p_20") # Get permutation tag
vec_perm_indexes <- 1:10 # Get vector of permutation indexes
vec_perm_indexes_tag <- paste0("p_", vec_perm_indexes) # Get vector of permutation indexes

################## For Age ######################

# Apply function to each permutation
temp_df_DVGs_perm_Age <- lapply(c("CD4_TCM", "CD4_Naive"), function(cell_type) {
  lapply(vec_perm_indexes_tag, function(perm_index_tag) get_DVGs_perm("cell_type", cell_type, perm_index_tag, "Age"))
})
temp_df_DVGs_perm_Age <- temp_df_DVGs_perm_Age 
DVGs_perm_Age <- bind_rows(temp_df_DVGs_perm_Age)

# Apply function to real DVGs
temp_df_DVGs_real_Age <- lapply(c("CD4_TCM", "CD4_Naive"), function(cell_type) {get_DVGs_real("cell_type", cell_type, "Age")})
temp_df_DVGs_real_Age <- temp_df_DVGs_real_Age
DVGs_real_Age <- bind_rows(temp_df_DVGs_real_Age)

# Summarize results for permuted and real dataframes and merge together 
DVGs_summary_Age <- DVGs_perm_Age %>% group_by(status, celltype) %>% summarize(num_DVGs = n()) # Get sum of DVGs for each permutation
num_DVGs_real_Age <- DVGs_real_Age %>% group_by(celltype) %>% summarize(real = n())
DVGs_summary_Age <- merge(DVGs_summary_Age, num_DVGs_real_Age, by = "celltype") # Get total DVGs counts across permutations for each celltype
DVGs_summary_Age <- DVGs_summary_Age %>% mutate(ratio = num_DVGs / real)  # Calculate ratio between DVG count vs total for each permutation
DVGs_summary_Age$phenotype <- "Age"

################## For Sex ######################

# Apply function to each permutation
temp_df_DVGs_perm_Sex <- lapply(c("CD4_TCM", "CD4_Naive"), function(cell_type) {
  lapply(vec_perm_indexes_tag, function(perm_index_tag) get_DVGs_perm("cell_type", cell_type, perm_index_tag, "Sex"))
})
temp_df_DVGs_perm_Sex <- temp_df_DVGs_perm_Sex
DVGs_perm_Sex <- bind_rows(temp_df_DVGs_perm_Sex)

# Apply function to real DVGs
temp_df_DVGs_real_Sex <- lapply(c("CD4_TCM", "CD4_Naive"), function(cell_type) {get_DVGs_real("cell_type", cell_type, "Sex")})
temp_df_DVGs_real_Sex <- temp_df_DVGs_real_Sex 
DVGs_real_Sex <- bind_rows(temp_df_DVGs_real_Sex)

# Summarize results for permuted and real dataframes and merge together 
DVGs_summary_Sex <- DVGs_perm_Sex %>% group_by(status, celltype) %>% summarize(num_DVGs = n()) # Get sum of DVGs for each permutation
num_DVGs_real_Sex <- DVGs_real_Sex %>% group_by(celltype) %>% summarize(real = n())
DVGs_summary_Sex <- merge(DVGs_summary_Sex, num_DVGs_real_Sex, by = "celltype") # Get total DVGs counts across permutations for each celltype
DVGs_summary_Sex <- DVGs_summary_Sex %>% mutate(ratio = num_DVGs / real)  # Calculate ratio between DVG count vs total for each permutation
DVGs_summary_Sex$phenotype <- "Sex"

############### Merging and plotting ################### 

DVGs_summary_all <- rbind(DVGs_summary_Age, DVGs_summary_Sex)

boxplot_permutations <- ggplot(DVGs_summary_all, aes(x = celltype, y = num_DVGs)) +
  geom_boxplot(outlier.shape = NA, fill = "grey95", color = "black", width = 0.5) +
  facet_wrap(~ fct_inorder(phenotype), nrow = 2, ncol = 1, scales = "free_x") +
  geom_point(aes(color = "Permutation"), size = 3, shape = 16, position = position_jitter(width = 0.1)) +
  geom_jitter(aes(y = real, color = "Main analysis"), size = 5, width = 0, height = 0, alpha = 1) +
  scale_color_manual(values = c("Permutation" = "coral", "Main analysis" = "#00B2B2"), guide = "legend") +
  scale_y_continuous(breaks = seq(0, ceiling(max(DVGs_summary_all$num_DVGs)), by = 10)) +
  labs(y = "Number of DVGs", x = NULL) +
  coord_flip() +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text = element_text(color = "black", size = 13),
    axis.title = element_text(size = 12),
    axis.text.x = element_text(size = 11), 
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.border = element_rect(),
    panel.background = element_blank(),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    strip.background = element_rect(color = "black", fill = "grey90", size = 1),
    strip.text.x = element_text(size = 10, face = "bold", color = "black"), 
    strip.placement = "outside",
    legend.key.size = unit(1, "cm")
  )

# Saving plot
boxplot.fn <- paste0("/home/lfronter/Desktop/Figures_TFM/Fig_3F.png")
CairoPNG(boxplot.fn, width = 8, height = 10, units = "in", res = 300)
print(boxplot_permutations)
dev.off()


############################################################################
## Boxplot for pseudobulked total, tech and bio variability for all genes ## 
############################################################################

cell_types <- c(
  "CD4_Naive", "CD4_TCM", "CD8_TEM", "NK", "B_naive", "CD8_Naive",
  "CD14_Mono", "B_intermediate", "Treg", "CD4_TEM", "B_memory", "CD16_Mono",
  "MAIT", "CD4_CTL", "CD8_TCM", "NK_CD56bright", "gdT", "cDC2",
  "Plasmablast", "NK_Proliferating", "dnT", "pDC", "Platelet", "HSPC",
  "Eryth", "CD4_Proliferating", "ILC", "CD8_Proliferating", "ASDC",
  "cDC1", "Doublet")

df_celltypes <- lapply(cell_types, function(cell_type) {
  
  # Get datasets
  res_dis.fn <- paste0("/home/lfronter/Documents/03_Results/Donor_variability/cell_type/", cell_type, "/Dispersion_5threshold.rds")
  if (file.exists(res_dis.fn)) {
    res_dis <- readRDS(res_dis.fn)
  } else {
    return(NULL)
  }
  
  res_dis_upd.fn <- paste0("/home/lfronter/Documents/03_Results/Donor_variability/cell_type/", cell_type, "/Dispersion_updated_5threshold_filtered.rds")
  if (file.exists(res_dis_upd.fn)) {
    res_dis_upd <- readRDS(res_dis_upd.fn)
  } else {
    return(NULL)
  }
  
  # Get mean for each gene
  res_dis_mean <- rowMeans(res_dis, na.rm = TRUE)
  res_dis_upd_mean <- rowMeans(res_dis_upd, na.rm = TRUE)
  
  # Get only rownames of the dataframe and add means 
  res_dis <- res_dis[, ncol(res_dis)]
  res_dis_df <- as.data.frame(res_dis)
  res_dis_df$mean <- res_dis_mean
  res_dis_df$type <- "res_dis"
  res_dis_df <- res_dis_df[, -1]
  
  res_dis_upd <- res_dis_upd[, ncol(res_dis_upd)]
  res_dis_upd_df <- as.data.frame(res_dis_upd)
  res_dis_upd_df$mean <- res_dis_upd_mean
  res_dis_upd_df$type <- "res_dis_upd" 
  res_dis_upd_df <- res_dis_upd_df[, -1]
  
  # Merge dataframes 
  df_merged <- rbind(res_dis_df, res_dis_upd_df)
  df_merged$celltype <- cell_type
  
  return(df_merged)
  
})

dis_all_onek1k <- do.call(rbind, df_celltypes)

# Plotting results
boxplot_comparison <- ggplot(dis_all_onek1k, aes(x = type, y = as.numeric(mean))) +
  geom_boxplot(width = 0.5,
               position = position_dodge(width = 2),
               outlier.shape = 21,
               outlier.size = 2,
               aes(fill = type),
               linewidth = 1) +
  facet_grid( ~ fct_inorder(celltype)) +
  theme_bw() +
  scale_y_log10(labels = scales::label_number(accuracy = 1)) +
  ylab("Dispersion (log10)") +
  scale_fill_manual(name = "Dispersion", 
                    labels = c("Biological", "Technical"),
                    values = c("res_dis_upd" = "#377EB8", "res_dis" = "#4DAF4A")) +
  theme(plot.title = element_text(size = 15),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=13, face = "bold"),
        text = element_text(size = 12),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        strip.text = element_text(size = 12, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Saving plot
boxplot.fn <- paste0("/home/lfronter/Desktop/Figures_TFM/Fig_3G.png")
CairoPNG(boxplot.fn, width = 25, height = 13, units = "in", res = 300)
print(boxplot_comparison)
dev.off()

