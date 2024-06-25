#!/usr/bin/env Rscript

########## SETTING DIRECTORIES AND LIBRARIES ############### 

# To fix the Matrix package in Seurat objects
Csparse_validate = 'CsparseMatrix_validate'

### Script to create boxplots based on raw expression files using different immune cell types from OneK1K database ###
                      ## Expression boxplots can be made using Age or Sex as grouping factor ##

# Activate required libraries

suppressMessages(library(dplyr))
suppressMessages(library(purrr))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(Seurat))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(forcats))
suppressMessages(library(grid))
suppressMessages(library(ggpubr))
suppressMessages(library(Cairo))
suppressMessages(library(optparse))
suppressMessages(library(cowplot))

############################## OPTIONS PARSER ###################################### 

option_list = list(
  make_option(c("--cell_level"), action="store", default=NA, type='character',
              help="From Azimuth: low resolution (predicted.celltype.l1) or high resolution (cell_type)"),
  make_option(c("--phenotype"), action="store", default="none", type='character',
              help="Age or Sex"),
  make_option(c("--cell_type"), action="store", default=NA, type='character',
              help="Cell types in low resolution (predicted.celltype.l1) or high resolution (cell_type)"))
opt = parse_args(OptionParser(option_list=option_list))

# Setting working directory (cluster or local)

path_cluster <- "/gpfs/projects/bsc83/"
path_local <- "/home/lfronter/Documents/"

if(file.exists(path_cluster)){
  setwd(paste(path_cluster))
}else if(file.exists(path_local)){
  setwd(paste(path_local))
}

# Setting path to stored CCV files (same as the output path)

CCV_path_cluster_cov <- "Projects/scRNAseq/lluisf/03_Results/DVA/"
CCV_path_local_cov <- "03_Results/DVA/"

# Setting output path

outpath_cluster_cov <- paste0(path_cluster, CCV_path_cluster_cov)
outpath_local_cov <- paste0(path_local, CCV_path_local_cov)

# Creating directory output (if needed)

if(!dir.exists(outpath_cluster_cov) & getwd() == path_cluster){dir.create(outpath_cluster_cov, recursive = T)}
if(!dir.exists(outpath_local_cov) & getwd() == path_local){dir.create(outpath_local_cov, recursive = T)}

##################################################
## Expression boxplots for Age and Sex on Scran ##
##################################################

############ FUNCTIONS ##############

plot_expr_sex <- function(gene, cell_level, cell_type, CCV, fdr, pval, pal, dataset = "OneK1K"){
  expr <- readRDS(paste0("/gpfs/projects/bsc83/Data/scRNAseq/Yazar2022/sce_data_objects/", cell_type, "_", cell_level, "_sceraw.rds"))

  # Adding "Sex" column from expression file
  males <- colnames(expr)[expr$Gender == "M"]
  females <- colnames(expr)[expr$Gender == "F"]

  # Preparing data
  expr_M_gene <- as.data.frame(t(as.matrix(counts(expr[gene,males]))))
  expr_M_gene$Sex <- "M"
  expr_F_gene <- as.data.frame(t(as.matrix(counts(expr[gene,females]))))
  expr_F_gene$Sex <- "F"
  expr_gene <- rbind(expr_F_gene, expr_M_gene)
  expr_gene$Sex <- factor(expr_gene$Sex, levels=c("M", "F"))
  count_F <- nrow(subset(expr_gene, Sex == "F"))
  count_M <- nrow(subset(expr_gene, Sex == "M"))
  colnames(expr_gene)[1] <- "gene_expression"

  # Plotting
  p <- ggplot(expr_gene, aes(x=Sex, y=gene_expression, fill=Sex)) + geom_violin(aes(fill=Sex), alpha = .75) +
    geom_boxplot(width=0.2, fill = "white", color = "black", coef = 0, outlier.alpha = 0) +
    ylab("Expression (log10)") + ggtitle(paste0(gene, " ", "(CD4_Naive)") ) + theme_bw() + scale_fill_manual(values = pal) +
    theme(plot.title = element_text(hjust = 0.5, size=50, face = "bold"), axis.text= element_text(size=40, color="black"),
          legend.position = "none", axis.title = element_text(size=45)) +
    scale_x_discrete(labels=c(paste0("F\n(n=", count_F, ")"), paste0("M\n(n=", count_M, ")")))

  # Adding metrics as text within plot
  median_F <- median(expr_gene$gene_expression[expr_gene$Sex == "F"])
  x_position <- "F"
  y_position <-  median_F - 1
  CCV_diff <- sprintf("CCV diff = %.3g", CCV)
  fdr <- sprintf("FDR = %.3g", fdr)
  pval <- sprintf("pval = %.3g", pval)
  label_text <- paste(CCV_diff, fdr, pval, sep = "\n")
  p <- p + geom_text(aes(x = "F", y = 1, hjust = 1, vjust = 0, label = label_text, fontface = "bold"), size = 13)
  return(p)
}

plot_expr_age_cat_all <- function(gene, cell_level, cell_type, CCV, fdr, pval, pal, dataset = "OneK1K"){
  expr <- readRDS(paste0("/gpfs/projects/bsc83/Data/scRNAseq/Yazar2022/sce_data_objects/", cell_type, "_", cell_level, "_sceraw.rds"))
  meta <- readRDS("/gpfs/projects/bsc83/Data/scRNAseq/Yazar2022/metadata_processed.rds")

  if(file.exists(path_cluster)){
    data <- readRDS(paste0(outpath_cluster_cov, "DVA_merged_DVG_cell_type_Age_cat.rds"))
  } else if (file.exists(path_local)){
    data <- readRDS(paste0(outpath_local_cov, "DVA_merged_DVG_cell_type_Age_cat.rds"))
  }

  # Adding "Age_cat_all" column from metadata to "sceraw" file
  inter_cells <- intersect(rownames(meta), colnames(expr))
  column_age_cat_all <- meta[rownames(meta) %in% inter_cells, "Age_cat_all"]
  colData(expr)[inter_cells, "Age_cat_all"] <- column_age_cat_all
  expr_Y <- colnames(expr)[expr$Age_cat_all == "1"] # 1 = Young
  expr_O <- colnames(expr)[expr$Age_cat_all == "2"] # 2 = Old

  # Preparing data
  expr_Y_gene <- as.data.frame(t(as.matrix(counts(expr[gene,expr_Y]))))
  expr_Y_gene$Age <- "Y"
  expr_O_gene <- as.data.frame(t(as.matrix(counts(expr[gene,expr_O]))))
  expr_O_gene$Age <- "O"
  expr_gene <- rbind(expr_Y_gene, expr_O_gene)
  expr_gene$Age <- factor(expr_gene$Age, levels=c("Y", "O"))
  count_Y <- nrow(subset(expr_gene, Age == "Y"))
  count_O <- nrow(subset(expr_gene, Age == "O"))
  colnames(expr_gene)[1] <- "gene_expression"

  # Plotting
  p <- ggplot(expr_gene, aes(x=Age, y=gene_expression, fill=Age)) + geom_violin(aes(fill=Age), alpha = .75) +
    geom_boxplot(width=0.2, fill = "white", color = "black", coef = 0, outlier.alpha = 0) +
    ylab("Expression (log10)") + ggtitle(paste0(gene, " ", "(",cell_type,")") ) + theme_bw() + scale_fill_manual(values = pal) +
    theme(plot.title = element_text(hjust = 0.5, size=50, face = "bold"), axis.text= element_text(size=40, color="black"),
          legend.position = "none", axis.title = element_text(size=45)) +
    scale_x_discrete(labels=c(paste0("Y\n(n=", count_Y, ")"), paste0("O\n(n=", count_O, ")")))

  # Adding metrics as text within plot
  median_O <- median(expr_gene$gene_expression[expr_gene$Age == "O"])
  x_position <- "O"
  y_position <-  median_O - 1
  CCV_diff <- sprintf("CCV diff = %.3g", CCV)
  fdr <- sprintf("FDR = %.3g", fdr)
  pval <- sprintf("pval = %.3g", pval)
  label_text <- paste(CCV_diff, fdr, pval, sep = "\n")
  p <- p + geom_text(aes(x = "Y", y = 1, hjust = 1, vjust = 0, label = label_text, fontface = "bold"), size = 13)
  return(p)
}

# Creating function to save previous boxplots (distinguishing between log-scale or not)

save_plot <- function(gene, cell, cell_level, pheno, CCV, fdr, pval, log = T, pal=c("#B19CD9", "#FFB6C1")){
  if(pheno == "Sex"){
    p <- plot_expr_sex(gene, cell_level, cell, CCV, fdr, pval, pal=c("#6A89E1", "#E08FE0"), dataset = dataset)
    if(log == T){
      return(p + scale_y_log10() + ylab("Expression (log10)"))
    } else {
      return(p + ylab("Expression"))
    }
  }else{
    p <- plot_expr_age_cat_all(gene, cell_level, cell, CCV, fdr, pval, pal, dataset = dataset)
    if(log == T){
      return(p + scale_y_log10() + ylab("Expression (log10)"))
    }else{
      return(p + ylab("Expression"))
    }
  }
}

############ ANALYSIS ###############

####################### WITH COVARIATES FOR AGE_CAT #############################

if(file.exists(path_cluster)){
  file <- paste0(outpath_cluster_cov, "DVA_merged_DVG_cell_type_Age_cat.rds")
} else if (file.exists(path_local)){
  file <- paste0(outpath_local_cov, "DVA_merged_DVG_cell_type_Age_cat.rds")
}

df_cov_age_cat <- readRDS(file)

## FACET FOR EXPRESSION BOXPLOT FOR DVGs in CD4_NAIVE SORTED BY CCV DIFFERENCE ##

df_CD4_Naive_cov_age <- subset(df_cov_age_cat, celltype == "CD4_Naive")
df_CD4_Naive_cov_age <- df_CD4_Naive_cov_age[order(df_CD4_Naive_cov_age$CCV_diff),]
df_CD4_Naive_cov_age <- head(df_CD4_Naive_cov_age, 3)

boxplot_list <- list()

for(row in 2:nrow(df_CD4_Naive_cov_age)){
  p <- save_plot(df_CD4_Naive_cov_age[row, 1], df_CD4_Naive_cov_age[row, 18], "cell_type", "Age_cat_all", df_CD4_Naive_cov_age[row, 14], df_CD4_Naive_cov_age[row, 17], df_CD4_Naive_cov_age[row, 16])
  boxplot_list[[df_CD4_Naive_cov_age[row, 1]]] <- p
}

## Saving first plot
if(file.exists(path_cluster)){
  plot.fn <- paste0(outpath_cluster_cov, "Boxplot_immune_genes/Age_cat/cell_type/CD4_Naive_boxplots_1.png")
}else if(file.exists(path_local)){
  plot.fn <- paste0(outpath_local_cov, "Boxplot_immune_genes/Age_cat/cell_type/CD4_Naive_boxplots_1.png")
}

CairoPNG(plot.fn, width = 20, height = 20, units = "in", res = 300)
print(boxplot_list[[1]])
dev.off()
print(paste("Boxplots for cell_type using Age_cat saved in", plot.fn))

## Saving second plot
if(file.exists(path_cluster)){
  plot.fn <- paste0(outpath_cluster_cov, "Boxplot_immune_genes/Age_cat/cell_type/CD4_Naive_boxplots_2.png")
}else if(file.exists(path_local)){
  plot.fn <- paste0(outpath_local_cov, "Boxplot_immune_genes/Age_cat/cell_type/CD4_Naive_boxplots_2.png")
}

CairoPNG(plot.fn, width = 20, height = 20, units = "in", res = 300)
print(boxplot_list[[2]])
dev.off()
print(paste("Boxplots for cell_type using Age_cat saved in", plot.fn))


####################### WITH COVARIATES FOR SEX #############################

if(file.exists(path_cluster)){
  file <- paste0(outpath_cluster_cov, "DVA_merged_DVG_cell_type_Sex.rds")
} else if (file.exists(path_local)){
  file <- paste0(outpath_local_cov, "DVA_merged_DVG_cell_type_Sex.rds")
}

df_cov_sex <- readRDS(file)

## FACET FOR EXPRESSION BOXPLOT FOR DVGs in CD4_NAIVE SORTED BY CCV DIFFERENCE ##

df_CD4_Naive_cov_sex <- subset(df_cov_sex, celltype == "CD4_Naive")
df_CD4_Naive_cov_sex <- df_CD4_Naive_cov_sex[order(df_CD4_Naive_cov_sex$CCV_diff),]
df_CD4_Naive_cov_sex <- head(df_CD4_Naive_cov_sex, 2)

boxplot_list <- list()

for(row in 1:nrow(df_CD4_Naive_cov_sex)){
  p <- save_plot(df_CD4_Naive_cov_sex[row, 1], df_CD4_Naive_cov_sex[row, 18], "cell_type", "Sex", df_CD4_Naive_cov_sex[row, 14], df_CD4_Naive_cov_sex[row, 17], df_CD4_Naive_cov_sex[row, 16])
  boxplot_list[[df_CD4_Naive_cov_sex[row, 1]]] <- p
}

## Saving first plot
if(file.exists(path_cluster)){
  plot.fn <- paste0(outpath_cluster_cov, "Boxplot_immune_genes/Sex/cell_type/CD4_Naive_boxplots_1.png")
}else if(file.exists(path_local)){
  plot.fn <- paste0(outpath_local_cov, "Boxplot_immune_genes/Sex/cell_type/CD4_Naive_boxplots_2.png")
}

CairoPNG(plot.fn, width = 20, height = 20, units = "in", res = 300)
print(boxplot_list[[1]])
dev.off()
print(paste("Boxplots for cell_type using Sex saved in", plot.fn))

## Saving second plot
if(file.exists(path_cluster)){
  plot.fn <- paste0(outpath_cluster_cov, "Boxplot_immune_genes/Sex/cell_type/CD4_TCM_boxplots_2.png")
}else if(file.exists(path_local)){
  plot.fn <- paste0(outpath_local_cov, "Boxplot_immune_genes/Sex/cell_type/CD4_TCM_boxplots_2.png")
}

CairoPNG(plot.fn, width = 20, height = 20, units = "in", res = 300)
print(boxplot_list[[2]])
dev.off()
print(paste("Boxplots for cell_type using Sex saved in", plot.fn))

#######################################################
## Expression boxplots for Age and Sex on pseudobulk ##
#######################################################

########## SETTING DIRECTORIES ############### 

path_cluster <- "/gpfs/projects/bsc83/Projects/scRNAseq/lluisf"
path_local <- "/home/lfronter/Documents"

if(startsWith(getwd(), path_local)){
  main.dir <- paste0(path_local, "/03_Results/Donor_variability/", opt$cell_level, "/")
} else {
  main.dir <- paste0(path_cluster, "/03_Results/Donor_variability/", opt$cell_level, "/")
}

results_dir <- paste0("Plots_", opt$phenotype, "/")

outdir <- paste0(main.dir, results_dir)
if(!dir.exists(outdir)){dir.create(outdir, recursive = T)}

######################
### 4.1. Functions ###
######################

plot_expr_age_cat_all <- function(gene, cell_level, cell_type, fdr, phenotype, logFC, pal){
  
  if(startsWith(getwd(), path_local)){
    expression <- readRDS(paste0("/home/lfronter/Downloads/scRNAseq/Yazar2022/sce_data_objects/", cell_type, "_", cell_level, "_sceraw.rds"))
  } else {
    expression <- readRDS(paste0("/gpfs/projects/bsc83/Data/scRNAseq/Yazar2022/sce_data_objects/", cell_type, "_", cell_level, "_sceraw.rds"))
  }
  
  if(startsWith(getwd(), path_local)){
    meta <- readRDS(paste0("/home/lfronter/Downloads/scRNAseq/Yazar2022/sce_data_objects/", cell_type, "_", cell_level, "_metadata.rds"))
  } else {
    meta <- readRDS(paste0("/gpfs/projects/bsc83/Data/scRNAseq/Yazar2022/sce_data_objects/", cell_type, "_", cell_level, "_metadata.rds"))
  }
  
  if(phenotype == "Age"){
    
    # Preparing metadata
    meta$Age_cat_all <- ifelse(meta$Age<=40, 'Y', 'O') # Categorize age
    inter_cells <- intersect(rownames(meta), colnames(expression))
    column_age_cat_all <- meta[rownames(meta) %in% inter_cells, "Age_cat_all"]
    SummarizedExperiment::colData(expression)[inter_cells, "Age_cat_all"] <- column_age_cat_all
    expr_Y <- as.vector(colnames(expression)[expression$Age_cat_all == "Y" ]) # Y = Young
    expr_O <- as.vector(colnames(expression)[expression$Age_cat_all == "O" ]) # O = Old
    
    # Preparing expression data
    expr_Y_gene <- base::as.data.frame(base::t(base::as.matrix(counts(expression[as.character(gene), expr_Y])))) 
    expr_Y_gene$Age <- "Y"
    expr_Y_gene$num_age <- meta$Age[match(expr_Y, rownames(meta))]
    expr_O_gene <- base::as.data.frame(base::t(base::as.matrix(counts(expression[as.character(gene), expr_O]))))
    expr_O_gene$Age <- "O"
    expr_O_gene$num_age <- meta$Age[match(expr_O, rownames(meta))]
    expr_gene <- rbind(expr_Y_gene, expr_O_gene)
    expr_gene$Age <- factor(expr_gene$Age, levels=c("Y", "O"))
    colnames(expr_gene)[1] <- "gene_expression"
    
    donors <- character(nrow(expr_gene))  # Get cells for which we have expression
    
    for (i in 1:nrow(expr_gene)) {
      barcode <- rownames(expr_gene)[i]  # Get the barcode from the current row name
      match_row <- match(barcode, rownames(meta))  # Find the corresponding row in metadata
      
      if (!is.na(match_row)) {
        donors[i] <- meta$assignment[match_row]  # Retrieve the donor value and store it in the vector
      } else {
        next
      }
    }
    
    expr_gene$donor <- donors # Add vector of donors
    
    expr_gene <- expr_gene[expr_gene$donor %in% donors_available, ] # Filter by donors which have variability values (from the other function)
    
    expr_gene <- expr_gene %>%
      mutate(gene_expression = as.numeric(gene_expression),
             Age = expr_gene$Age,
             num_age = as.numeric(num_age),
             donor = as.character(donor)) # Make class columns appropiate
    
    expr_gene_donor <- expr_gene %>%
      group_by(donor) %>%
      mutate(total_expression = base::sum(as.numeric(gene_expression))) %>%
      ungroup() # Get total expression per donor
    
    expr_gene_donor$gene_expression <- NULL # Remove cell expression
    expr_gene_donor <- unique(expr_gene_donor) # Get unique donors
    expr_gene_donor <- expr_gene_donor %>% group_by(Age) %>%
      mutate(cumulative_expression = cumsum(total_expression)) # Get cumulative expression for each group
    expr_gene_donor$Age <- factor(expr_gene_donor$Age, levels = rev(levels(expr_gene_donor$Age))) # Relevel Age column
    
    count_Y <- nrow(subset(expr_gene_donor, Age == "Y")) # Get number of "Y" donors
    count_O <- nrow(subset(expr_gene_donor, Age == "O")) # Get number of "O" donors
    
    # Plotting
    p <- ggplot(expr_gene_donor, aes(x=Age, y=log2(total_expression), fill=Age)) + geom_violin(aes(fill=Age), alpha = .75) +
      geom_boxplot(width=0.2, fill = "grey95", color = "black", coef = 1.5, outlier.alpha = 0) +
      xlab("Age") + ggtitle(paste0(gene, " ", "(",cell_type,")") ) + theme_bw() + scale_fill_manual(values = pal[[1]]) +
      theme(plot.title = element_text(hjust = 0.5, size=50, face = "bold"), axis.text= element_text(size=40, color="black"),
            legend.position = "none", axis.title = element_text(size=45),  plot.tag = element_text(), panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), panel.background = element_blank()) +
      scale_x_discrete(labels=c(paste0("O\n(n=", count_O, ")"), paste0("Y\n(n=", count_Y, ")")))
    
    # Adding metrics as text within plot
    median_O <- median(expr_gene_donor$total_expression[expr_gene_donor$Age == "O"])
    x_position <- "O"
    y_position <-  median_O - 1
    fdr <- paste("fdr =", signif(as.numeric(fdr), digits = 3))
    logFC <- paste("logFC =", signif(as.numeric(logFC), digits = 4))
    label_text <- paste(fdr, logFC, sep = "\n")
    p <- p + geom_text(aes(x = "Y", y = 1, hjust = 1, vjust = 0, label = label_text, fontface = "bold"), size = 13)
    
    # Plotting cumulative plot 
    cum_plot <- ggplot(expr_gene_donor, aes(x = cumulative_expression, color = Age)) +
      stat_ecdf(geom = "step", size = 1) +
      xlab("Total expression") +
      ylab("f(Total expression)") +
      ggtitle(paste0(gene, " ", "(",cell_type,")")) +
      scale_color_manual(values = pal[[1]]) +
      theme_bw() +
      theme(axis.text = element_text(size=12, color="black"),
            plot.title = element_text(hjust = 0.5, size=20, face = "bold"),
            axis.title = element_text(size=16),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank())
    
    return(list(boxplot = p, cumulative_plot = cum_plot))
    
  } else {
    
    # Preparing metadata
    colnames(meta)[colnames(meta) == "Gender"] <- "Sex"
    inter_cells <- intersect(rownames(meta), colnames(expression))
    column_sex <- meta[rownames(meta) %in% inter_cells, "Sex"]
    SummarizedExperiment::colData(expression)[inter_cells, "Sex"] <- column_sex
    expr_M <- as.vector(colnames(expression)[expression$Sex == "M" ]) # M = Male
    expr_F <- as.vector(colnames(expression)[expression$Sex == "F" ]) # F = Female
    
    # Preparing expression data
    expr_M_gene <- base::as.data.frame(base::t(base::as.matrix(counts(expression[as.character(gene), expr_M])))) 
    expr_M_gene$Sex <- "M"
    expr_F_gene <- base::as.data.frame(base::t(base::as.matrix(counts(expression[as.character(gene), expr_F]))))
    expr_F_gene$Sex <- "F"
    expr_gene <- rbind(expr_M_gene, expr_F_gene)
    expr_gene$Sex <- factor(expr_gene$Sex, levels=c("M", "F"))
    colnames(expr_gene)[1] <- "gene_expression"
    
    donors <- character(nrow(expr_gene))  # Get vector on donors (empty)
    
    for (i in 1:nrow(expr_gene)) {
      barcode <- rownames(expr_gene)[i]  # Get the barcode from the current row name
      match_row <- match(barcode, rownames(meta))  # Find the corresponding row in metadata
      
      if (!is.na(match_row)) {
        donors[i] <- meta$assignment[match_row]  # Retrieve the donor value and store it in the vector
      } else {
        next
      }
    }
    
    expr_gene$donor <- donors # Add vector of donors
    
    expr_gene <- expr_gene[expr_gene$donor %in% donors_available, ] # Filter by donors which have variability values (from the other function)
    
    expr_gene <- expr_gene %>%
      mutate(gene_expression = as.numeric(gene_expression),
             Sex = expr_gene$Sex,
             donor = as.character(donor)) # Make class columns appropiate
    
    expr_gene_donor <- expr_gene %>%
      group_by(donor) %>%
      mutate(total_expression = base::sum(as.numeric(gene_expression))) %>%
      ungroup() # Get total expression per donor
    
    expr_gene_donor$gene_expression <- NULL # Remove cell expression
    expr_gene_donor <- unique(expr_gene_donor) # Get unique donors
    expr_gene_donor <- expr_gene_donor %>% group_by(Sex) %>%
      mutate(cumulative_expression = cumsum(total_expression)) # Get cumulative expression for each group
    
    count_M <- nrow(subset(expr_gene_donor, Sex == "M")) # Get number of "M" donors
    count_F <- nrow(subset(expr_gene_donor, Sex == "F")) # Get number of "F" donors
    
    # Plotting
    p <- ggplot(expr_gene_donor, aes(x=Sex, y=log2(total_expression), fill=Sex)) + geom_violin(aes(fill=Sex), alpha = .75) +
      geom_boxplot(width=0.2, fill = "grey95", color = "black", coef = 1.5, outlier.alpha = 0) +
      xlab("Sex") + ggtitle(paste0(gene, " ", "(",cell_type,")") ) + theme_bw() + scale_fill_manual(values = pal[[2]]) +
      theme(plot.title = element_text(hjust = 0.5, size=50, face = "bold"), axis.text= element_text(size=40, color="black"),
            legend.position = "none", axis.title = element_text(size=45),  plot.tag = element_text(), panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), panel.background = element_blank()) +
      scale_x_discrete(labels=c(paste0("M\n(n=", count_M, ")"), paste0("F\n(n=", count_F, ")")))
    
    # Adding metrics as text within plot
    median_M <- median(expr_gene_donor$total_expression[expr_gene_donor$Sex == "M"])
    x_position <- "M"
    y_position <-  median_M - 1
    fdr <- paste("fdr =", signif(as.numeric(fdr), digits = 3))
    logFC <- paste("logFC =", signif(as.numeric(logFC), digits = 4))
    label_text <- paste(fdr, logFC, sep = "\n")
    p <- p + geom_text(aes(x = "F", y = 1, hjust = 1, vjust = 0, label = label_text, fontface = "bold"), size = 13)
    
    # Plotting cumulative plot 
    cum_plot <- ggplot(expr_gene_donor, aes(x = cumulative_expression, color = Sex)) +
      stat_ecdf(geom = "step", size = 1) +
      xlab("Total expression") +
      ylab("f(Total expression)") +
      ggtitle(paste0(gene, " ", "(",cell_type,")")) +
      scale_color_manual(values = pal[[2]]) +
      theme_bw() +
      theme(axis.text = element_text(size=12, color="black"),
            plot.title = element_text(hjust = 0.5, size=20, face = "bold"),
            axis.title = element_text(size=16),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank())
    
    return(list(boxplot = p, cumulative_plot = cum_plot))
    
  }
}

plot_var_age_cat_all <- function(gene, cell_level, cell_type, fdr, phenotype, beta, pal){
  
  var_data <- readRDS(paste0(main.dir, cell_type, "/Dispersion_updated_5threshold.rds")) # Get residual dispersion data
  var_data_gene <- subset(var_data, rownames(var_data) == gene) # Filter data from gene of interest
  var_data_gene <- as.data.frame(t(var_data_gene)) # Trasverse dataframe
  var_data_gene$donor <- rownames(var_data_gene)
  rownames(var_data_gene) <- NULL
  var_data_gene <- var_data_gene[complete.cases(var_data_gene), ]
  donors <- var_data_gene$donor # Get donors 
  donors_available <<- var_data_gene$donor
  
  # Add age status
  if(startsWith(getwd(), path_local)){
    meta <- readRDS(paste0("/home/lfronter/Downloads/scRNAseq/Yazar2022/sce_data_objects/", cell_type, "_", cell_level, "_metadata.rds"))
  } else {
    meta <- readRDS(paste0("/gpfs/projects/bsc83/Data/scRNAseq/Yazar2022/sce_data_objects/", cell_type, "_", cell_level, "_metadata.rds"))
  }
  
  if(phenotype == "Age"){
    
    age_donors <- as.data.frame(meta) %>% # Get dataframe containing donor and age 
      filter(assignment %in% donors) %>%
      distinct(assignment, Age) %>%
      mutate(donor = assignment, age = Age) %>%
      select(donor, age)
    rownames(age_donors) <- NULL
    var_gene_donors <- merge(age_donors, var_data_gene, by = "donor")
    colnames(var_gene_donors)[colnames(var_gene_donors) == gene] <- "total_expression" # Change colname of total expression column
    var_gene_donors$Age_cat_all <- ifelse(var_gene_donors$age<=40, 'Y', 'O') # Categorize age
    var_gene_donors <- var_gene_donors %>% group_by(Age_cat_all) %>%
      mutate(cumulative_expression = cumsum(total_expression)) # Get cumulative expression for each group
    
    count_Y <- nrow(subset(var_gene_donors, Age_cat_all == "Y")) # Get number of "Y" donors
    count_O <- nrow(subset(var_gene_donors, Age_cat_all == "O")) # Get number of "O" donors
    
    # Plotting
    p <- ggplot(var_gene_donors, aes(x=Age_cat_all, y=log2(total_expression), fill=Age_cat_all)) + geom_violin(aes(fill=Age_cat_all), alpha = .75) +
      geom_boxplot(width=0.2, fill = "grey95", color = "black", coef = 1.5, outlier.alpha = 0) +
      xlab("Age") + ggtitle(paste0(gene, " ", "(",cell_type,")") ) + theme_bw() +
      scale_fill_manual(values = pal[[1]]) + 
      theme(plot.title = element_text(hjust = 0.5, size=50, face = "bold"), axis.text= element_text(size=40, color="black"),
            legend.position = "none", axis.title = element_text(size=45),  plot.tag = element_text(), panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), panel.background = element_blank()) +
      scale_x_discrete(labels=c(paste0("O\n(n=", count_O, ")"), paste0("Y\n(n=", count_Y, ")")))
    
    # Adding metrics as text within plot
    median_O <- median(var_gene_donors$beta[var_gene_donors$Age_cat_all == "O"])
    x_position <- "O"
    y_position <-  median_O - 1
    fdr <- paste("fdr =", signif(as.numeric(fdr), digits = 3))
    beta <- paste("beta =", signif(as.numeric(beta), digits = 4))
    label_text <- paste(fdr, beta, sep = "\n")
    p <- p + geom_text(aes(x = "Y", y = 1, hjust = 1, vjust = 0, label = label_text, fontface = "bold"), size = 13)
    
    # Plotting cumulative plot 
    cum_plot <- ggplot(var_gene_donors, aes(x = cumulative_expression, color = Age_cat_all)) +
      stat_ecdf(geom = "step", size = 1) +
      xlab("Total variability") +
      ylab("f(Total variability)") +
      ggtitle(paste0(gene, " ", "(",cell_type,")")) +
      scale_color_manual(values = pal[[1]]) +
      theme_bw() +
      theme(axis.text = element_text(size=12, color="black"),
            plot.title = element_text(hjust = 0.5, size=20, face = "bold"),
            axis.title = element_text(size=16),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank())
    
    return(list(boxplot = p, cumulative_plot = cum_plot))
    
  } else {
    
    colnames(meta)[colnames(meta) == "Gender"] <- "Sex"
    sex_donors <- as.data.frame(meta) %>% # Get dataframe containing donor and sex 
      filter(assignment %in% donors) %>%
      distinct(assignment, Sex) %>%
      mutate(donor = assignment, sex = Sex) %>%
      select(donor, sex)
    rownames(sex_donors) <- NULL
    var_gene_donors <- merge(sex_donors, var_data_gene, by = "donor")
    colnames(var_gene_donors)[colnames(var_gene_donors) == gene] <- "total_expression" # Change colname of total expression column
    var_gene_donors <- var_gene_donors %>% group_by(sex) %>%
      mutate(cumulative_expression = cumsum(total_expression), sex = as.factor(sex)) # Get cumulative expression for each group
    var_gene_donors$sex <- factor(var_gene_donors$sex, levels = rev(levels(var_gene_donors$sex))) # Relevel sex column
    
    count_M <- nrow(subset(var_gene_donors, sex == "M")) # Get number of "M" donors
    count_F <- nrow(subset(var_gene_donors, sex == "F")) # Get number of "F" donors
    
    # Plotting
    p <- ggplot(var_gene_donors, aes(x=sex, y=log2(total_expression), fill=sex)) + geom_violin(aes(fill=sex), alpha = .75) +
      geom_boxplot(width=0.2, fill = "grey95", color = "black", coef = 1.5, outlier.alpha = 0) +
      xlab("Sex") + ggtitle(paste0(gene, " ", "(",cell_type,")") ) + theme_bw() +
      scale_fill_manual(values = pal[[2]]) + 
      theme(plot.title = element_text(hjust = 0.5, size=50, face = "bold"), axis.text= element_text(size=40, color="black"),
            legend.position = "none", axis.title = element_text(size=45),  plot.tag = element_text(),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
      scale_x_discrete(labels=c(paste0("M\n(n=", count_M, ")"), paste0("F\n(n=", count_F, ")")))
    
    # Adding metrics as text within plot
    median_M <- median(var_gene_donors$beta[var_gene_donors$sex == "M"])
    x_position <- "M"
    y_position <-  median_M - 1
    fdr <- paste("fdr =", signif(as.numeric(fdr), digits = 3))
    beta <- paste("beta =", signif(as.numeric(beta), digits = 4))
    label_text <- paste(fdr, beta, sep = "\n")
    p <- p + geom_text(aes(x = "F", y = 1, hjust = 1, vjust = 0, label = label_text, fontface = "bold"), size = 13)
    
    # Plotting cumulative plot 
    cum_plot <- ggplot(var_gene_donors, aes(x = cumulative_expression, color = sex)) +
      stat_ecdf(geom = "step", size = 1) +
      xlab("Total variability") +
      ylab("f(Total variability)") +
      ggtitle(paste0(gene, " ", "(",cell_type,")")) +
      scale_color_manual(values = pal[[2]]) +
      theme_bw() +
      theme(axis.text = element_text(size=12, color="black"),
            plot.title = element_text(hjust = 0.5, size=20, face = "bold"),
            axis.title = element_text(size=16),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank())
    
    return(list(boxplot = p, cumulative_plot = cum_plot))
  }
}

palette_age <- c("#B19CD9", "#FFB6C1")
palette_sex <- c("#6A89E1", "#E08FE0")
palette_plotting <- list(palette_age, palette_sex)

save_plot_expr <- function(gene, cell_level, cell, fdr, phenotype, logFC, log, pal=palette_plotting){
  result <- plot_expr_age_cat_all(gene, cell_level, cell, fdr, phenotype, logFC, pal)
  p <- result$boxplot
  cum_plot <- result$cumulative_plot
  if(log == T){
    p <- p + ylab("Expression (log2)")
    return(list(boxplot = p, cum_plot = cum_plot))
  }else{
    p <- p + ylab("Expression")
    return(list(boxplot = p, cum_plot = cum_plot))
  }
}

save_plot_var <- function(gene, cell_level, cell, fdr, phenotype, beta, log, pal=palette_plotting){
  result <- plot_var_age_cat_all(gene, cell_level, cell, fdr, phenotype, beta, pal)
  p1 <- result$boxplot
  cum_plot <- result$cumulative_plot
  if(log == T){
    p1 <- p1 + ylab("Variability (log2)")
    return(list(boxplot = p1, cum_plot = cum_plot))
  }else{
    p1 <- p1 + ylab("Variability")
    return(list(boxplot = p1, cum_plot = cum_plot))
  }
}

# Get significant genes in variability analysis
DVGs_df <- readRDS(paste0("/gpfs/projects/bsc83/Projects/scRNAseq/lluisf/03_Results/Donor_variability/cell_type/DVGs_Sex_all_cell_type_5threshold_test.rds"))
DVGs_sig <- DVGs_df %>%
  filter(p_value < 0.05) # Without beta threshold filtering 
gene_row_var <- DVGs_sig %>% filter(gene == "ATP5C1" & celltype == "CD4_Naive")

# Get significant genes in expression analysis
expr_data <- readRDS(paste0("/gpfs/projects/bsc83/Projects/scRNAseq/lluisf/03_Results/Donor_variability/cell_type/DEG_Sex_cell_type.rds"))
expr_data$celltype <- sapply(strsplit(expr_data$celltype, " "), function(x) paste(x, collapse = "_")) # Get proper cell type naming
gene_row_expr <- expr_data %>% filter(gene == "ATP5C1" & celltype == "CD4_Naive")

# Execute boxplots for variability and expression analysis 
p_var <- save_plot_var(gene_row_var$gene, opt$cell_level, opt$cell_type, gene_row_var$p_value, opt$phenotype, gene_row_var$beta, log = T) ## Plotting variability boxplot 
p_expr <- save_plot_expr(gene_row_expr$gene, opt$cell_level, opt$cell_type, gene_row_expr$adj.P.Val, opt$phenotype, gene_row_expr$logFC, log = T) ## Plotting expression boxplot

# Saving boxplots
figfn <- paste0(outdir, "Boxplots_definitive_log2.png")
CairoPNG(figfn, width = 20, height = 20, units = "in", res = 300)
plot_dvg_not_deg <- plot_grid(p_expr$boxplot, p_var$boxplot, labels = "AUTO", align = "h") # Merge both plots
title <- ggdraw() +
  draw_label("DVG not DEG", size = 28, fontface = 'bold')
plot_merged <- plot_grid(title, plot_dvg_not_deg, ncol = 1, rel_heights = c(0.1, 1)) # Merge plots and title
print(plot_merged)
dev.off()