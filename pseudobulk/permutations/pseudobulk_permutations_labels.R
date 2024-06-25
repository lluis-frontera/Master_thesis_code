#!/usr/bin/env Rscript

# To fix the Matrix package in Seurat objects
Csparse_validate = 'CsparseMatrix_validate'

# Setting working directory (cluster or local)

path_cluster <- "/gpfs/projects/bsc83/"
path_local <- "/home/lfronter/Documents/"

if(startsWith(getwd(), path_local)){
  main.dir <- path_local
} else {
  main.dir <- path_cluster
}

# Options parser
shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list = list(
  make_option(c("--dataset"), action="store", default=NULL, type='character',
              help="v2, v3, pilot3 or OneK1K"),
  make_option(c("--cell_level"), action="store", default=NULL, type='character',
              help="predicted.celltype.l1 or cell_type"),
  make_option(c("--cell_type"), action="store", default=NULL, type='character',
              help="CD4_T/NK..."),
  make_option(c("--phenotype"), action="store", default=NULL, type='character',
              help="Sex/Age/Age_cat/Age_cat_all"),
  make_option(c("-n", "--n_permutations"), action="store", default=20, type='integer',
              help="Number of permutations"),
  make_option(c("--in_dir"), action="store", default='sce_data_objects/', type='character',
              help="Main directory"),
  make_option(c("--out_dir"), action="store", default='03_Results/', type='character',
              help="Output main directory"))
opt = parse_args(OptionParser(option_list=option_list))

# Packages
shhh(library(plyr))
shhh(library(dplyr))
shhh(library(reshape2))
shhh(library(stringi))
shhh(library(stringr))
shhh(library(ggplot2))
shhh(library(Seurat))
shhh(library(SingleCellExperiment))
shhh(library(SummarizedExperiment))

#################### Set Variables and load Data #################### 
# Dataset
# opt$dataset <- 'v2'
# opt$dataset <- 'v3'
# opt$dataset <- 'pilot3'

# Cell level
# opt$cell_level <- 'predicted.celltype.l1'
# opt$cell_level <- 'cell_type'

# Cell type
# opt$cell_type <- 'CD4_T'
# opt$cell_type <- 'NK'

# Phenotype
# opt$phenotype <- 'CMV_status'
# opt$phenotype <- 'Gender'
# opt$phenotype <- 'Age'
# opt$phenotype <- 'Age_cat'
# opt$phenotype <- 'Age_cat_all'

# Permutations
n_perm.ch <- as.character(opt$n_permutations)
n_perm.tag <- paste0('p_', n_perm.ch)

# Input data
if(opt$dataset == "OneK1K"){
  main.dir.in <- paste0(main.dir, 'Data/scRNAseq/Yazar2022/')
  in.dir <- paste0(main.dir.in, opt$in_dir)
  in.fn <- paste0(in.dir, opt$cell_type, '_', opt$cell_level, '_metadata.rds') 
} else {
  in.fn <- paste0(main.dir, "/Projects/scRNAseq/aripol1/wijst-2020-hg19/v1/aging/seurat/v2_sct_azimuth_TL.lastrelease.rds")
}

# Output directory
out.dir <- paste0(main.dir, 'Projects/scRNAseq/lluisf/')
if(opt$phenotype == "Gender"){ 
  out.dir <- paste0(out.dir, opt$out_dir, "/Permutations/", n_perm.tag, '/', opt$dataset, "/Sex/", opt$cell_level, '/', opt$cell_type, '/')
} else {
  out.dir <- paste0(out.dir, opt$out_dir, "/Permutations/", n_perm.tag, '/', opt$dataset, '/', opt$phenotype, "/", opt$cell_level, '/', opt$cell_type, '/')
  }
if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}
print(paste0('Main results directory: ', out.dir))

# Report
print('############################')
print(paste0('Dataset: ', opt$dataset))
print(paste0('Cell level: ', opt$cell_level))
print(paste0('Cell type: ', opt$cell_type))
print(paste0('Phenotype: ', opt$phenotype))
print(paste0('# of permutations: ', n_perm.ch))
print(paste0('Input data: ', in.fn))
print(paste0('Output directory: ', out.dir))
print('############################')
cat('\n')

#################### Set Variables and load Data ####################
# Read metadata
cell_md <- readRDS(in.fn)

# Categorize Age
cell_md$Age_cat <- ifelse(cell_md$Age<=40, 'Y', ifelse(cell_md$Age>=60, 'O', 'M'))
cell_md$Age_cat_all <- ifelse(cell_md$Age<=40, 'Y', 'O')

# Get informative columns
rnames <- rownames(cell_md)
if(opt$dataset == "OneK1K"){
  donor_vars <- c('assignment','date', 'Gender','Age', 'Age_cat', 'Age_cat_all')
} else {
  donor_vars <- c('assignment','batch', 'Gender','Age', 'Age_cat', 'Age_cat_all')
}
donor_md <- cell_md[,donor_vars]
rownames(donor_md) <- NULL

# New phenotype variable with shuffled labels
perm_vec <- seq(1:opt$n_permutations)
perm_names <- paste0('p_', perm_vec)

#################### Define functions ####################

# i <- perm_vec[1]
# df = donor_md
# phenotype = opt$phenotype
perm.func <- function(i, df = donor_md, phenotype = opt$phenotype){
  print(paste0('# Permutation number: ', as.character(i)))
  
  # Create random vector
  set.seed(i)
  levels_vec <- df[[phenotype]]
  levels_vec.shuffled <- sample(levels_vec)
  
  # Add random vector to the new column
  phe_random <- paste0(phenotype,'.shuffled')
  df[phe_random] <- levels_vec.shuffled
  match.n <- table(df[[phenotype]]==df[[phe_random]]) # Get logical table of previous dataframe and shuffled one
  match.n_df <- as.data.frame(match.n)
  colnames(match.n_df) <- c('match', 'n')
  match.prop <- prop.table(match.n) # Same as before, but with proportions
  match.prop_df <- as.data.frame(match.prop)
  colnames(match.prop_df) <- c('match', 'prop')
  match.df <- merge(match.n_df, match.prop_df)
  out <- list(df_shuffled = df,
              df_match = match.df)
  return(out) # Get two lists, one for the shuffled metadata and the other for the matching with the original metadata
}

# i <- stats_vars[1]
# df = df_match
# title_var = ggtitle_var
bp.func <- function(i, df = df_match, title_var = ggtitle_var){
  yintercept_var <- 0.5
  if(i=='n'){
    n_total <- sum(df_match[df_match$perm=='p_1',]$n)
    n_half <- round(n_total/2)
    yintercept_var <- n_half
  }
  p <- ggplot(df_match, aes(x = match, y = .data[[i]])) +
    geom_boxplot() + theme_bw() + 
    ggtitle(ggtitle_var) + 
    theme(plot.title = element_text(face="bold", hjust = 0.5)) + 
    geom_hline(yintercept=yintercept_var)
  p.fn <- paste0(out.dir, 'match_boxplot.', i, '.png')
  ggsave(p.fn, p, width = 4, height = 5)
}

#################### Analysis ####################

## Apply function

perm.res <- lapply(perm_vec, function(x) perm.func(x))
names(perm.res) <- perm_names

## Extract results

### Shuffled labels
perm.list <- lapply(perm.res, function(x) x$df_shuffled)
perm.list_fn <- paste0(out.dir, 'permutations.donor_metadata.list.rds')
print(paste0('Saving the donor metadata from the permutations in: ', perm.list_fn))
saveRDS(perm.list, perm.list_fn)

### Match
df_match <- do.call("rbind", lapply(perm.res, function(x) x$df_match))
df_match$perm <- str_split_fixed(rownames(df_match), '\\.', 2)[,1]
rownames(df_match) <- NULL
df_match <- df_match[,c('perm','match','n','prop')]
df_match <- df_match[order(df_match$perm, -df_match$prop),]
df_match_fn <- paste0(out.dir, 'match.donor_metadata.df.rds')
print(paste0('Saving the donor metadata match from the permutations in: ', df_match_fn))
saveRDS(df_match, df_match_fn)

## Boxplots

## Parameters
ggtitle_var <- paste0(opt$dataset, ' - ', opt$cell_type, ' - ', opt$phenotype)
original_levels <- table(donor_md[[opt$phenotype]])
original_levels.vec <- paste(names(original_levels), unname(original_levels), sep = ':')
if(length(original_levels.vec)<=10){
  original_levels.tag <- paste(original_levels.vec, collapse=', ')
  ggtitle_var <- paste0(ggtitle_var, '\n(', original_levels.tag, ')')
}

## Apply function
stats_vars <- c('n', 'prop')
bp.res <- lapply(stats_vars, function(x) bp.func(x))
