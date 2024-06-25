#!/usr/bin/env Rscript

# To fix the Matrix package in Seurat objects
Csparse_validate = 'CsparseMatrix_validate'

# Setting working directory (cluster or local)
path_cluster <- '/gpfs/projects/bsc83/'
path_em <- '/home/aripol1/Desktop/bsc/'
path_opensuse <- '/home/bscuser/bsc/'

if(file.exists(path_cluster)){
  setwd(paste(path_cluster))
}else if(file.exists(path_em)){
  setwd(paste(path_em))
}else if(file.exists(path_opensuse)){
  setwd(paste(path_opensuse))
}

# Options parser
shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list = list(
  make_option(c("--dataset"), action="store", default=NA, type='character',
              help="pilot3/v2/v3/OneK1K"),
  make_option(c("--cell_level"), action="store", default=NA, type='character',
              help="From Azimuth: low resolution (predicted.celltype.l1) or high resolution (cell_type)"),
  make_option(c("--cell_type"), action="store", default=NA, type='character',
              help="Cell types in low resolution (predicted.celltype.l1) or high resolution (cell_type)"),
  make_option(c("--phenotype"), action="store", default="none", type='character',
              help="Age_cat/Age_cat_all/Sex/CMV_status"),
  make_option(c("--covariates"), action="store", default="none", type='character',
              help="'24.cell_to_cell_variability.covariates.tab' if defined, 'Sex/Age/assignment/date' if not defined, must be comma separated"),
  make_option(c("--out_dir"), action="store", default='24.cell_to_cell_variability', type='character',
              help="Output main directory"))
opt = parse_args(OptionParser(option_list=option_list))

# Packages
shhh(library(SingleCellExperiment))
shhh(library(SummarizedExperiment))
shhh(library(scran))
shhh(library(scuttle))
shhh(library(Seurat))
shhh(library(reshape2))
shhh(library(stringi))
shhh(library(stringr))
shhh(library(ggplot2))
shhh(library(ggpubr))

################################## Set variables and load data ################################## 

# Main directories
pilot3_v2_v3.dir <- 'Projects/scRNAseq/aripol1/wijst-2020-hg19/v1/aging/00.so_split_by_celltype/' 
onek1k.dir <- '/gpfs/projects/bsc83/Data/scRNAseq/Yazar2022/sce_data_objects'
onek1kagecorrected.dir <- 'Data/scRNAseq/Yazar2022/sca_data_objects'

# opt$dataset <- 'v2' 
# opt$dataset <- 'v3'
# opt$dataset <- 'pilot3'
# opt$dataset <- 'OneK1K'
# opt$dataset <- 'OneK1Kagecorrected'

# opt$phenotype <- 'Age_cat' #3 categories
# opt$phenotype <- 'Age_cat_all' #2 categories
# opt$phenotype <- 'Sex'
# opt$phenotype <- 'CMV_status'

# cell_level
# opt$cell_level <- 'predicted.celltype.l1'
# opt$cell_level <- 'cell_type' #'predicted.celltype.l2' (in OneK1K)

# cell_type
# opt$cell_type <- 'CD4_T' #if cell_level=predicted.celltype.l1
# opt$cell_type <- 'CD8_TEM' #if cell_level=cell_type
# opt$cell_type <- 'CD4_CTL' #if cell_level=cell_type
# opt$cell_type <- 'CD4_TCM' #if cell_level=cell_type
# opt$cell_type <- 'Mono' #if cell_level=predicted.celltype.l1
# opt$cell_type <- 'CD14_Mono' #if cell_level=cell_type
# opt$cell_type <- 'B' #if cell_level=predicted.celltype.l1
# opt$cell_type <- 'other_T' #if cell_level=predicted.celltype.l1
# opt$cell_type <- 'CD8_TCM'
# opt$cell_type <- 'NK_CD56bright'

# opt$covariates <- '24.cell_to_cell_variability.covariates.tab' # If not defined
# opt$covariates <- 'Sex/Age/assignment/date' # If defined, and comma-separated

# NOTE: 

# If opt$phenotype is "Sex", there is no need to add "Sex" as a covariate. If opt$phenotype is "Age_cat" or "Age_cat_all", there is no need to add "Age" as a covariate.

if(grepl('OneK1K',opt$dataset)){
  if(opt$dataset=='OneK1K'){in.dir <- onek1k.dir}
  if(opt$dataset=='OneK1Kagecorrected'){in.dir <- onek1kagecorrected.dir}
  cell_type <- opt$cell_type 
  if(opt$cell_level=='predicted.celltype.l1'){
    cell_type.fn <- paste0(cell_type, "_", opt$cell_level, '_sceraw.rds') 
  }else{
    cell_type.fn <- paste0(cell_type, "_", opt$cell_level, '_sceraw.rds')
  }
  cell_type.fn <- paste0(in.dir, '/', cell_type.fn)
}else{
  cell_type <- opt$cell_type
  in.dir <- pilot3_v2_v3.dir
  cell_type.fn <- paste0(in.dir, '/', opt$dataset, '/', opt$cell_level, '/', cell_type, '.rds') 
}

# Output
main.dir <- 'Projects/scRNAseq/aripol1/wijst-2020-hg19/v1/aging/'
out.dir <- paste0(path_cluster, "Projects/scRNAseq/lluisf/", opt$out_dir, "/", opt$dataset, "/")
covariates_vec <- NULL
if(opt$phenotype!='none'){
  out.dir <- paste0(out.dir, opt$phenotype, '/')
  if(grepl("\\.tab$",opt$covariates)){ 
    covariates.fn <- paste0(main.dir, 'scripts/', opt$covariates)
    covariates_vec <- read.table(covariates.fn)$V1
    covariates_tag <- paste(covariates_vec, collapse = '.')
    cov_to_remove <- NULL
    if(opt$phenotype=='Sex'){cov_to_remove <- opt$Sex} 
    if(opt$phenotype%in%c('Age_cat', 'Age_cat_all')){cov_to_remove <- 'Age'}
    covariates_vec <- covariates_vec[!covariates_vec%in%cov_to_remove]
    out.dir <- paste0(out.dir, covariates_tag, '/')
  } else {
    covariates_vec <- unlist(strsplit(opt$covariates, split = ","))
    covariates_tag <- paste(covariates_vec, collapse = ".") 
    cov_to_remove <- NULL
    if(opt$phenotype=='Sex'){cov_to_remove <- opt$Sex}
    if(opt$phenotype%in%c('Age_cat', 'Age_cat_all')){cov_to_remove <- 'Age'}
    covariates_vec <- covariates_vec[!covariates_vec%in%cov_to_remove]
    out.dir <- paste0(out.dir, covariates_tag, '/')
  }
}

out.dir <- paste0(out.dir, opt$cell_level, '/', opt$cell_type, "/")
if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}
print(paste0('Main output directory: ',out.dir))

# Report
print(paste0('Dataset: ', opt$dataset))
print(paste0('Cell level: ', opt$cell_level))
print(paste0('Cell type: ', opt$cell_type))
print(paste0('Phenotype: ', opt$phenotype))
print(paste0('Covariates: ', opt$covariates))
print(paste0('Input dir: ', in.dir))
print(paste0('Output dir: ', out.dir))

################################## Functions ##################################

# Filter out lowly expressed genes 
# USED: In their paper (more relaxed) -->  based on proportion of 0's (5% threshold) 
# sce <- sce_object
# na.rm = TRUE
# freq_expressed = 0.05
filter_by_prop_zeros <- function(sce, na.rm = TRUE, freq_expressed = 0.05){
  # Extract logcounts matrix from SCE object
  print('Extracting logcounts matrix from SCE object...')
  counts <- logcounts(sce)
  
  # Select genes
  sc_mat <- t(counts) > 0
  out <- sparseMatrixStats::colMeans2(sc_mat, rows = NULL, cols = NULL, na.rm = na.rm)
  names(out) <- colnames(sc_mat)
  genes_kept <- names(out[out>freq_expressed])
  
  # Filter SCE
  sce_filt <- sce[rownames(sce)%in%genes_kept,]
  
  return(sce_filt)
}

# Scran: SCE log normalized + variance modelling
# sce <- sce_object
# sce <- sce_object.list[[1]]
# covs <- covariates_vec
# phenotype = opt$phenotype
# block_var = 'assignment'
variance_modelling <- function(sce, covs = NULL, phenotype = opt$phenotype, block_var = 'assignment'){
  # Check nCells per donor
  md <- droplevels(as.data.frame(colData(sce)))
  ncells_per_donor <- sort(table(md[[block_var]]),decreasing = T)
  ncells_per_donor.summary <- summary(as.numeric(ncells_per_donor))
  donors_kept <- names(ncells_per_donor[ncells_per_donor>5])
  md.donors_kept <- droplevels(md[md[[block_var]]%in%donors_kept,])
  bc_kept <- rownames(md.donors_kept)
  nDonors_total <- length(ncells_per_donor)
  nDonors_kept <- length(donors_kept)
  propDonors_kept <- (nDonors_kept/nDonors_total)*100
  nDonors_lost <- nDonors_total - nDonors_kept
  propDonors_lost <- (nDonors_lost/nDonors_total)*100
  nBarcodes_total <- length(colnames(sce))
  nBarcodes_kept <- length(bc_kept)
  propBarcodes_kept <- (nBarcodes_kept/nBarcodes_total)*100
  nBarcodes_lost <- nBarcodes_total - nBarcodes_kept
  propBarcodes_lost <- (nBarcodes_lost/nBarcodes_total)*100
  
  print('Summary:')
  print(ncells_per_donor.summary)
  print(paste0('nDonors total: ', nDonors_total))
  print(paste0('nDonors kept: ', nDonors_kept))
  print(paste0('nCells total: ', nBarcodes_total))
  print(paste0('nCells kept: ', nBarcodes_kept))
  print(paste0('propCells kept: ', round((length(bc_kept)/length(colnames(sce)))*100,2)))
  
  nCells_per_donor.df <- data.frame(nCells_per_Donor.min = min(ncells_per_donor),
                                    nCells_per_Donor.max = max(ncells_per_donor),
                                    nCells_per_Donor.mean = mean(ncells_per_donor),
                                    nCells_per_Donor.median = median(ncells_per_donor),
                                    nDonors_total = nDonors_total,
                                    nDonors_kept = nDonors_kept,
                                    propDonors_kept = propDonors_kept,
                                    nDonors_lost = nDonors_lost,
                                    propDonors_lost = propDonors_lost,
                                    nBarcodes_total = nBarcodes_total,
                                    nBarcodes_kept = nBarcodes_kept,
                                    propBarcodes_kept = propBarcodes_kept,
                                    nBarcodes_lost = nBarcodes_lost,
                                    propBarcodes_lost = propBarcodes_lost)
  
  print('Subsetting nCells...')
  sce <- sce[,colnames(sce)%in%bc_kept]
  colData(sce)$Sex <- droplevels(as.factor(colData(sce)$Gender))
  colData(sce)$date <- droplevels(as.factor(colData(sce)$date))
  colData(sce)$assignment <- droplevels(as.factor(colData(sce)$assignment))
  
  # DF with Variance estimates
  print('modelGeneVar..')
  if(!is.null(covs)){
    covs_to_regress <- covs[covs!=phenotype]
    form_covs.string <- paste0('~',paste(covs_to_regress,collapse='+'))
    form_covs <- as.formula(form_covs.string)
    print(form_covs)
    metadata <- droplevels(as.data.frame(colData(sce)))
    design_mat <- model.matrix(form_covs, metadata)
    n_cols.full <- ncol(design_mat)
    dec <- try(modelGeneVar(sce, design=design_mat))
    if(class(dec)=='try-error'){
      # Change design_mat
      ## Error in .ranksafeQR(design) : design matrix is not of full rank 
      ## .ranksafeQR() (https://rdrr.io/bioc/scuttle/src/R/fitLinearModel.R) called in modelGeneVar() (https://rdrr.io/bioc/scran/src/R/modelGeneVar.R), inside .model_gene_var >  .compute_mean_var
      print('Remove linearly dependent columns from the design matrix...')
      design_mat <- design_mat[,qr(design_mat)$pivot[seq_len(qr(design_mat)$rank)]]
      n_cols.subset <- ncol(design_mat)
      print(paste0('Final columns: ', as.character(n_cols.subset), '/', as.character(n_cols.full)))
      dec <- try(modelGeneVar(sce, design=design_mat))
    }
  }else{
    dec <- try(modelGeneVar(sce, block=sce[[block_var]], equiweight=FALSE))
  }
  if(class(dec)=='try-error'){dec <- modelGeneVar(sce)}
  
  # Get genes
  print('getTopHVGs...')
  
  # Get the top 10% of genes.
  top.hvgs_top10 <- getTopHVGs(dec, prop=0.1)
  
  # Get the top 2000 genes.
  top.hvgs_top2000 <- getTopHVGs(dec, n=2000)
  
  # Get all genes with positive biological components.
  top.hvgs_bio <- getTopHVGs(dec, var.threshold=0)
  
  # Get all genes with FDR below 5%.
  top.hvgs_fdr <- getTopHVGs(dec, fdr.threshold=0.05)
  
  # Output
  dec.df <- as.data.frame(dec[,c(1:6)])
  dec.df$gene_id <- rownames(dec.df)
  dec_fdr.df <- dec.df[rownames(dec.df)%in%top.hvgs_fdr,]
  dec_pos.df <- dec.df[rownames(dec.df)%in%top.hvgs_bio,]
  dec_top10prop.df <- dec.df[rownames(dec.df)%in%top.hvgs_top10,]
  dec_top2000abs.df <- dec.df[rownames(dec.df)%in%top.hvgs_top2000,]
  
  out <- list(dec = dec,
              dec_all = dec.df,
              dec_fdr = dec_fdr.df,
              dec_pos = dec_pos.df,
              dec_top10prop = dec_top10prop.df,
              dec_top2000abs = dec_top2000abs.df,
              nCells_per_donor = nCells_per_donor.df)
  return(out)
}

# Scran: CCV difference (old vs. young) --> zscore --> p-values (only if !is.null(opt$phenotype))
# vm_list = variance_modelling.subsetted
# group1 = 'Y'
# group2 = 'O'
# tag = 'dec_all'
# tag = 'dec_fdr'
# tag = 'dec_pos'

CCV_diff <- function(vm_list, group1, group2, tag){
  # Prepare data
  print('Preparing data...')
  print(paste0('Tag: ', tag))
  print(paste0('Group1: ', group1))
  print(paste0('Group2: ', group2))
  df_1 <- vm_list[[group1]][[tag]]
  colnames(df_1) <- paste0(colnames(df_1), '.', group1)
  group1_bio <- paste0('bio.', group1)
  
  df_2 <- vm_list[[group2]][[tag]]
  colnames(df_2) <- paste0(colnames(df_2), '.', group2)
  group2_bio <- paste0('bio.', group2)
  
  df_1$gene_id <- rownames(df_1)
  df_2$gene_id <- rownames(df_2)
  
  # CCV diff --> z-score --> p-value
  print('Merging group1 and group2 variance estimates...')
  system.time(df_merged <- merge(df_1, df_2, by = 'gene_id'))
  
  df_merged.selected <- NULL
  if(nrow(df_merged)>0){
    print('Computing CCV diff --> z-score --> p-values...')
    df_merged$CCV_diff <- df_merged[[group2_bio]] - df_merged[[group1_bio]]
    df_merged$z_score <- (df_merged$CCV_diff-mean(df_merged$CCV_diff))/sd(df_merged$CCV_diff)
    df_merged$z_score.pval <- pnorm(df_merged$z_score)
    df_merged$z_score.fdr <- p.adjust(df_merged$z_score.pval, 'fdr')
    
    # Variables of interest
    common_metrics <- 'gene_id'
    group_metrics <- c('mean', 'total', 'tech', 'bio', 'p.value', 'FDR')
    group1_metrics <- paste0(group_metrics, '.', group1)
    group2_metrics <- paste0(group_metrics, '.', group2)
    diff_metrics <- c('CCV_diff', 'z_score', 'z_score.pval', 'z_score.fdr')
    selected_metrics <- c(common_metrics, group1_metrics, group2_metrics, diff_metrics)
    df_merged.selected <- df_merged[,selected_metrics]
    df_merged.selected <- df_merged.selected[order(df_merged.selected$z_score.pval),]
    df_merged.selected <- as.data.frame(df_merged.selected)
    
    # Check
    df_merged.selected.pval_ss <- df_merged.selected[df_merged.selected$z_score.pval<0.05,]
    df_merged.selected.fdr_ss <- df_merged.selected[df_merged.selected$z_score.fdr<0.05,]
    vars <- c(group1_metrics, group2_metrics, diff_metrics)
    sapply(vars, function(i) summary(df_merged.selected[[i]]), simplify = FALSE)
    
    n.sign_pval <- nrow(df_merged.selected.pval_ss)
    n.sign_fdr <- nrow(df_merged.selected.fdr_ss)
    print(paste0('# of total genes in both groups: ', nrow(df_merged.selected)))
    print(paste0('# of significant (pval) genes: ', n.sign_pval))
    print(paste0('# of significant (fdr) genes: ', n.sign_fdr))
  }
 
  return(df_merged.selected)
}

# CV (sd/mean*100) on SCE log normalized
CV <- function(x){
  (sd(x)/mean(x))*100
}
func_to_list <- function(l, f){
  vapply(l, f, FUN.VALUE=0.0)
}

cv.func <- function(sce, margin.idx){
  # Extract logcounts matrix from SCE object
  print('Extracting logcounts matrix from SCE object...')
  sce_logcounts <- logcounts(sce)
  
  # Split matrix into list by rows
  print('Splitting matrix to list by rows...')
  system.time(gene_vec.list <- asplit(sce_logcounts, margin.idx))
  
  # Apply CV function to the list of gene vectors
  print('Applying CV function to the list of gene vectors...')
  system.time(out <- func_to_list(gene_vec.list, CV))
  out <- sort(out, decreasing=T)
  
  return(out)
}

################################## Analyses ##################################

print('Read data...')
sce_object <- readRDS(cell_type.fn)

# Categorize Age
colData(sce_object)$Age_cat <- ifelse(colData(sce_object)$Age<=40, 'Y',
                        	ifelse(colData(sce_object)$Age>=60, 'O', 'M'))
colData(sce_object)$Age_cat_all <- ifelse(colData(sce_object)$Age<=40, 'Y', 'O')

colData(sce_object)$Sex <- as.factor(colData(sce_object)$Gender)
colData(sce_object)$date <- as.factor(colData(sce_object)$date)
colData(sce_object)$assignment <- as.factor(colData(sce_object)$assignment)

if(opt$phenotype=='CMV_status'){
  # CMV metadata
  ## Read data
  cmv_fn <- 'Data/scRNAseq/LifeLinesDeep_CMV/Metadata_matchLLD_CMV.tsv'
  cmv_df <- read.delim(cmv_fn, check.names = FALSE)
  
  ## Subset by sc-donors
  aggregate_metadata.original <- as.data.frame(colData(sce_object))
  sc_donors <- unique(aggregate_metadata.original$assignment)
  print(paste0('# of sc donors: ', length(sc_donors)))
  cmv_df.sc <- droplevels(cmv_df[cmv_df$LLDEEP_ID%in%sc_donors,])
  print(paste0('# of match donors: ', nrow(cmv_df.sc)))
  cmv_df.sc <- cmv_df.sc[!is.na(cmv_df.sc$CMV_Baseline) | !is.na(cmv_df.sc$CMV_Followup),]
  print(paste0('# of match donors (without NAs in CMV_Baseline and CMV_Followup): ', nrow(cmv_df.sc)))
  
  ## CMV status
  cmv_df.sc$CMV_status <- ifelse(!is.na(cmv_df.sc$CMV_Followup), cmv_df.sc$CMV_Followup, cmv_df.sc$CMV_Baseline)
  cmv_df.sc$Age_diff_SCvsLLD <- ifelse(!is.na(cmv_df.sc$CMV_Followup), cmv_df.sc$Age_diff_SCvsLLD2, cmv_df.sc$Age_diff_SCvsLLD1)
  table(cmv_df.sc$CMV_status)
  summary(cmv_df.sc$Age_diff_SCvsLLD)
  
  ## Add CMV metadata to seurat metadata
  aggregate_metadata <- merge(aggregate_metadata.original, cmv_df.sc, by.x = 'assignment', by.y = 'LLDEEP_ID')
  aggregate_metadata$CMV_status <- as.factor(as.character(aggregate_metadata$CMV_status))
  rownames(aggregate_metadata) <- aggregate_metadata$bare_barcode_lane
  aggregate_metadata <- droplevels(aggregate_metadata)
  
  ## Subset sce_object by final cells with CMV_status annotation (rownames of aggregate_metadata)
  bc <- rownames(aggregate_metadata)
  sce_object <- sce_object[,colnames(sce_object)%in%bc]

  ## add proper metadata
  cnames_diff <- setdiff(colnames(aggregate_metadata), colnames(colData(sce_object)))
  aggregate_metadata.subset <- aggregate_metadata[,cnames_diff]
  aggregate_metadata.subset$bare_barcode_lane <- rownames(aggregate_metadata.subset)
  sce_object.md <- merge(colData(sce_object), aggregate_metadata.subset, by = 'bare_barcode_lane')
  sce_object.md <- sce_object.md[match(colnames(sce_object), sce_object.md$bare_barcode_lane),]
  identical(sce_object.md$bare_barcode_lane, colnames(sce_object))
  rownames(sce_object.md) <- sce_object.md$bare_barcode_lane
  identical(rownames(sce_object.md), colnames(sce_object))
  colData(sce_object) <- sce_object.md
}

# Filter out lowly expressed genes (in their paper)
print('Filtering out lowly expressed genes (proportion of 0s>0.05)...')
system.time(sce_object <- filter_by_prop_zeros(sce_object))

# Save tested genes
fn <- paste0(out.dir, opt$cell_type, '.bg.txt')
print(paste0('Saving background genes (tested) in: ', fn))
write.table(as.data.frame(rownames(sce_object)), fn, quote = FALSE, row.names = FALSE, col.names = FALSE)

# Cell-to-cell variability metrics
if(opt$phenotype!='none'){
  print(paste0('Phenotype: ', opt$phenotype))

  # Split SCE object by opt$phenotype
  print('Splitting SCE object by metadata variable...')
  barcodes_by_phenotypes.list <- lapply(split(colData(sce_object), colData(sce_object)[[opt$phenotype]]), function(x) rownames(x))
  sce_object.list <- lapply(barcodes_by_phenotypes.list, function(x){sce_subset <- sce_object[,which(colnames(sce_object)%in%x)]})
  
  ## 1. Regression modelling (Scran)
  ### 1.1. Variance modelling on SCE log normalized
  print('Scran: variance modelling on SCE log normalized...')
  system.time(variance_modelling.list <- lapply(sce_object.list, function(x) variance_modelling(x, covariates_vec)))
  variance_modelling.all <- lapply(variance_modelling.list, function(x) x[['dec']])
  variance_modelling.subsetted <- lapply(variance_modelling.list, function(x) x[c('dec_all', 'dec_fdr', 'dec_pos', 'dec_top10prop', 'dec_top2000abs')])
  variance_modelling.nCells_per_donor <- lapply(variance_modelling.list, function(x) x[['nCells_per_donor']])
  variance_modelling.nCells_per_donor.df <- do.call("rbind", variance_modelling.nCells_per_donor)
  variance_modelling.nCells_per_donor.df[[opt$phenotype]] <- str_split_fixed(rownames(variance_modelling.nCells_per_donor.df), '\\.', 2)[,1]
  
  ### 1.2. CCV_diff --> z-scores --> p-values
  #### compute
  if(opt$phenotype=='CMV_status'){Group_order <- c('0','1')}
  if(opt$phenotype=='Sex'){Group_order <- c('M','F')}
  if(opt$phenotype%in%c('Age_cat', 'Age_cat_all')){Group_order <- c('Y','O')}
  group1_var <- Group_order[1]
  group2_var <- Group_order[2]
  tag.vec <- c('dec_all', 'dec_fdr', 'dec_pos', 'dec_top10prop', 'dec_top2000abs')
  system.time(CCV_diff.list <- sapply(tag.vec, function(i) CCV_diff(vm_list = variance_modelling.subsetted, 
                                                                    group1 = group1_var,
                                                                    group2 = group2_var,
                                                                    tag = i), simplify = FALSE))
  lapply(CCV_diff.list, function(x) nrow(x[x$z_score.fdr<=0.05,])) #check
  CCV_diff.fn <- paste0(out.dir, opt$cell_type, '.CCV_diff.rds')
  print(paste0('Saving CCV difference (df) in: ', CCV_diff.fn))
  saveRDS(CCV_diff.list, CCV_diff.fn)
  
  #### Report
  CCV_diff.report <- sapply(names(CCV_diff.list), function(i){
    print(i)
    
    # Save genes
    x <- CCV_diff.list[[i]]
    if(!is.null(x)){
      x$CCV_diff.direction <- ifelse(x$CCV_diff>0, 'increase', 'decrease')
      x_list <- split(x, x$CCV_diff.direction)
      save_ccv.diff_genes <- sapply(names(x_list), function(j){
        print(j)
        x <- x_list[[j]]
        fn <- paste0(out.dir, opt$cell_type, '.', i, '.', j, '.txt')
        write.table(as.data.frame(x$gene_id), fn, quote = FALSE, row.names = FALSE, col.names = FALSE)
      }, simplify = FALSE)
      
      # Wilcoxon test
      group1_bio <- paste0('bio.', group1_var)
      group2_bio <- paste0('bio.', group2_var)
      xx <- x[,c('gene_id', group1_bio, group2_bio)]
      wt <- wilcox.test(xx[[group1_bio]], xx[[group2_bio]])
      xx.melt <- melt(xx, id.vars = 'gene_id')
      y_var <- 'value'
      if(min(xx.melt$value)>0){
        xx.melt$value_log10 <- log10(xx.melt$value)
        y_var <- 'value_log10'
      }
      p <- ggboxplot(xx.melt, x = "variable", y = y_var,
                     color = "variable", palette = "jco",
                     add = "jitter", add.params = list(alpha = 0.4))
      p <- p + stat_compare_means()
      p.fn <- paste0(out.dir, opt$cell_type, '.', i, '.png')
      print(paste0('Saving boxplot (Wilcoxon test) in: ', p.fn))
      ggsave(p.fn, p, width = 4.5, height = 5)
    }else{
      print('No results.')
    }
  
  }, simplify = FALSE)
  
  ## 2. Generics --> CV (sd/mean*100) on SCE log normalized
  print('CV on SCE log normalized...')
  system.time(cv.list <- lapply(sce_object.list, function(x) cv.func(x, 1)))
}else{
  ## 1. Regression modelling --> Scran: variance modelling on SCE log normalized
  print('Scran: variance modelling on SCE log normalized...')
  system.time(variance_modelling.list <- variance_modelling(sce_object))
  variance_modelling.all <- variance_modelling.list[['dec']]
  variance_modelling.subsetted <- variance_modelling.list[c('dec_all', 'dec_fdr', 'dec_pos', 'dec_top10prop', 'dec_top2000abs')]
  variance_modelling.nCells_per_donor.df <- variance_modelling.list[['nCells_per_donor']]
  
  ## 2. Generics --> CV (sd/mean*100) on SCE log normalized
  print('CV on SCE log normalized...')
  system.time(cv.list <- cv.func(sce_object, 1))
}

# Saving common outputs
## Variance modelling: !is.null(opt$phenotype)=list // !is.null(opt$phenotype)=df
variance_modelling.fn <- paste0(out.dir, opt$cell_type, '.variance_modelling.rds')
print(paste0('Saving variance modelling (object) in: ', variance_modelling.fn))
system.time(saveRDS(variance_modelling.all, variance_modelling.fn))

variance_modelling.subsetted.fn <- paste0(out.dir, opt$cell_type, '.variance_modelling.subsetted.rds')
print(paste0('Saving variance modelling (subsetted) in: ', variance_modelling.subsetted.fn))
system.time(saveRDS(variance_modelling.subsetted, variance_modelling.subsetted.fn))

variance_modelling.report.fn <- paste0(out.dir, opt$cell_type, '.variance_modelling.report.txt')
print(paste0('Saving variance modelling (report) in: ', variance_modelling.report.fn))
write.table(variance_modelling.nCells_per_donor.df, variance_modelling.report.fn, quote = FALSE, row.names = FALSE, sep = '\t')

## CV: !is.null(opt$phenotype)=list // !is.null(opt$phenotype)=df
cv.fn <- paste0(out.dir, opt$cell_type, '.cv.rds')
print(paste0('Saving CV in: ', cv.fn))
saveRDS(cv.list, cv.fn)

