#!/usr/bin/env Rscript

########## SETTING LIBRARIES ############### 

shhh <- suppressPackageStartupMessages
shhh(library(tidyr))
shhh(library(parallel))
shhh(library(purrr))
shhh(library(reshape2))
shhh(library(qqman))
shhh(library(qvalue))
shhh(library(ggplot2))
shhh(library(cowplot))
shhh(library(grid))
shhh(library(gridExtra))
shhh(library(gtable))
shhh(library(ggsignif))
shhh(library(pheatmap))
shhh(library(corrplot))
shhh(library(RColorBrewer))
shhh(library(viridis))
shhh(library(optparse))
shhh(library(stringi))
shhh(library(dplyr))
shhh(library(ggrastr))
shhh(library(SingleCellExperiment))

############################## OPTIONS PARSER ###################################### 

option_list = list(
  make_option(c("--cell_level"), action="store", default=NA, type='character',
              help="From Azimuth: low resolution (predicted.celltype.l1) or high resolution (cell_type)"),
  make_option(c("--covariates"), action="store", default="none", type='character',
              help="'24.cell_to_cell_variability.covariates.tab' if defined, 'Sex/Age/date' if not defined, must be comma separated"),
  make_option(c("--perm_num"), action="store", default=20, type='integer',
              help="Number of permutations (100, 1000, etc..)"),
  make_option(c("--perm_idx"), action="store", default=NULL, type='integer',
              help="Number of permutation (1:1000)"),
  make_option(c("--phenotype"), action="store", default="none", type='character',
              help="Age or Sex"),
  make_option(c("--cell_type"), action="store", default=NA, type='character',
              help="Cell types in low resolution (predicted.celltype.l1) or high resolution (cell_type)"))
opt = parse_args(OptionParser(option_list=option_list))

options(digits = 4)

########## SETTING DIRECTORIES ############### 

path_cluster <- "/gpfs/projects/bsc83/Projects/scRNAseq/lluisf"
path_local <- "/home/lfronter/Documents"

perm_idx.tag <- paste0('p_', as.character(opt$perm_idx))
perm_num.tag <- paste0('p_', as.character(opt$perm_num))
cell_level <- opt$cell_level
cell_type <- opt$cell_type
phenotype <- opt$phenotype
covariates <- opt$covariates

if(startsWith(getwd(), path_local)){
  main.dir <- paste0(path_local, "/03_Results/Donor_variability/cell_type/B_intermediate/") # CHANGE MANUALLY WHEN LOCAL
} else {
  main.dir <- paste0(path_cluster, "/03_Results/Donor_variability/Permutations/", perm_num.tag, "/",  perm_idx.tag, "/", cell_level, "/", cell_type, "/")
}

outdir <- paste0(main.dir, "Results_", phenotype, "/")
if(!dir.exists(outdir)){dir.create(outdir, recursive = T)}

####################### Gene expression variance and differential analysis using OneK1K pseudobulk ##############################
                    ######## Based on script by Julong Wei, modified by Lluis Frontera ##########

##############################
### Differential functions ###
##############################

### myDE, DEA for a single contrast (batch together)
myDE <- function(y, gene, design_matrix) {
  
  # y: Vector of gene expression values (filtered to avoid NA values)
  # gene: Gene identifier
  # design_matrix: Design matrix for covariates
  
  y <- try(log2(y), silent = TRUE)
  lm0 <- try(lm(y ~ 0 + design_matrix), silent = TRUE) # Fit linear model
  
  summary_fit <- summary(lm0)
  
  # Extract summary statistics from model
  if (class(y) != "try-error" && class(lm0) != "try-error") {
    b <- coef(lm0) # Get coefficients
    if(phenotype == "Age"){
      coef_interest <- "design_matrixAge"
    } else {
      coef_interest <- "design_matrixSexM"
    }
    b <- b[coef_interest] # Extract the one regarding age
    vb <- diag(vcov(lm0)) # Get diagonal of covariance matrix
    vb <- vb[coef_interest] # Extract the one regarding age
    std_error <- summary_fit$coefficients[, "Std. Error"]  # Get std error
    std_error <- std_error[coef_interest] # Extract the one regarding age
    z_score <- summary_fit$coefficients[, "t value"]  # Get t value
    z_score <- z_score[coef_interest] 
    p_value <- summary_fit$coefficients[, "Pr(>|t|)"]  # Get p-value
    p_value <- p_value[coef_interest] # Get p-value of coefficient of interest

  # Extract own summary statistics
    bhat <- b
    sdhat <- sqrt(vb)
    z <- bhat / sdhat
    p <- 2 * pnorm(-abs(z))
    dd <- data.frame(gene = gene, beta = bhat, stderr_own = sdhat, z_score_own = z, p_value_own = p, 
                     stderr = std_error, z_score = z_score, p_value = p_value) # Merge all into dataframe 
  } else {
    dd <- NA
  }
}

#################################################################
### 1. Differential dispersion after removing mean effects ###
#################################################################

####################################
### 1.1 Differential procedure ###
####################################

## Differential residual dispersion

# Read dispersion and expression data

if(startsWith(getwd(), path_local)){
  PhxNew <- readRDS(paste0("/home/lfronter/Documents/03_Results/Donor_variability/", cell_level, "/", cell_type, "/Dispersion_updated_5threshold.rds"))
} else {
  PhxNew <- readRDS(paste0("/gpfs/projects/bsc83/Projects/scRNAseq/lluisf/03_Results/Donor_variability/", cell_level, "/", cell_type, "/Dispersion_updated_5threshold.rds"))
}

if(startsWith(getwd(), path_local)){
  sce_object <- readRDS(paste0("/home/lfronter/Downloads/scRNAseq/Yazar2022/sce_data_objects/", cell_type, "_", cell_level, "_sceraw.rds"))
} else {
  sce_object <- readRDS(paste0("/gpfs/projects/bsc83/Data/scRNAseq/Yazar2022/sce_data_objects/", cell_type, "_", cell_level, "_sceraw.rds"))
}

# Get permuted metadata
perm.fn <- paste0(path_cluster, "/03_Results/Permutations/", perm_num.tag, "/OneK1K/", phenotype, "/",
                  cell_level, '/', cell_type, '/permutations.donor_metadata.list.rds')
perm.list <- readRDS(perm.fn)
perm.md <- perm.list[[opt$perm_idx]]
perm_var <- grep('shuffled', colnames(perm.md), value = TRUE)
donor_permvar.md <- perm.md[,c('assignment', perm_var)]
perm_var.renamed <- gsub('.shuffled', '', perm_var)
colnames(donor_permvar.md)[colnames(donor_permvar.md)==perm_var] <- perm_var.renamed

# Add permuted metadata
df <- colData(sce_object)
df_removed <- df[,-which(colnames(df)==perm_var.renamed)]
df_permvar <- merge(df_removed, donor_permvar.md, by = 'assignment')
rnames_idx <- match(rownames(colData(sce_object)), df_permvar$bare_barcode_lane)
df_permvar <- df_permvar[rnames_idx,]
rownames(df_permvar) <- df_permvar$bare_barcode_lane
identical(rownames(df_permvar), rownames(colData(sce_object)))
colData(sce_object) <- df_permvar
sce_object.Raw <- sce_object

print(paste("Rownames of the orginal dataframe matching the order of the permuted dataframe", identical(rownames(df_permvar), rownames(colData(sce_object)))))

# Create design matrix
covariates_vec <- unlist(strsplit(covariates, split = ",")) # Get covariates vector
form_covs.string <- paste('~', phenotype, "+", paste0(covariates_vec, collapse = " + "))
form_covs <- as.formula(form_covs.string) # Get linear model
print(form_covs) # Print linear model

metadata <- droplevels(as.data.frame(colData(sce_object))) # Get single-cell metadata
donors <- purrr::set_names(levels(as.factor(metadata$assignment))) # Get donors
n_cells <- as.numeric(table(metadata$assignment)) # Get number of cells per donor
m <- match(donors, sce_object$assignment) # Reorder the samples (rows) of the metadata to match the order of sample names in donors vector
metadata_donor <- data.frame(colData(sce_object)[m, ],
                             n_cells, row.names = NULL) # Get pseudobulked metadata
colnames(metadata_donor)[colnames(metadata_donor) == "Gender"] <- "Sex" # Change column name for convenience
rownames(metadata_donor) <- metadata_donor$assignment # Set donors as rownames

design_mat <- model.matrix(form_covs, metadata_donor) # Get pseudobulked design matrix (GENERAL ONE)
donors_available <- intersect(rownames(design_mat), colnames(PhxNew)) # Get available donors for cell type (some of them were trimmed during filtering)
design_mat <- design_mat[donors_available, ] # Get donors available
design_mat <- design_mat[,-1 ] # Remove intercept column

## Checking if colnames from design matrix are sorted as rownames from dispersion file

print(paste("Rownames of design matrix matching colnames of dispersion file:", all(rownames(design_mat) == colnames(PhxNew)))) # TRUE

### 1.1 Estimate differential results by donor
print(paste("1.1 Differential analysis by donor, genes analyzed:", nrow(PhxNew)))
rn <- rownames(PhxNew) # Get gene names
## Start loop by gene
TMP <- mclapply(rn, function(ii){
  y <- PhxNew[ii, ] # Get expression vector for gene
  y <- y[!is.na(y)] # Get only expression values for donors in which we have value (not NA)
  design_mat <- design_mat[names(y), ] # Adapt design matrix to get only donors for which we have expression values
  dd <- myDE(y, ii, design_mat)
  dd
}, mc.cores=1) ### End loop by gene
TMP <- TMP[!sapply(TMP, is.null)] # Filter NULL values
TMP <- as.data.frame(do.call(rbind, TMP)) %>% mutate(celltype =  cell_type)
TMP$p_value_adj_own <- p.adjust(TMP$p_value_own,  method = "BH") # Add own adjusted p-value 
TMP$p_value_adj <- p.adjust(TMP$p_value,  method = "BH")   # Add adjusted p-value 
rownames(TMP) <- NULL

opfn <- paste0(outdir, "1.DVG_results_5threshold_test.rds")
saveRDS(TMP, file = opfn)
print(paste("DVA results for", cell_type, "saved in", opfn))

###########################
### 2. Summary results ###
###########################

print("Summary results")

############################
### 2.1. Histogram plots ###
############################

DVGs <- readRDS(paste0(outdir, "1.DVG_results_5threshold_test.rds")) # Open file
DVGs <- DVGs %>% mutate(significance=ifelse(p_value < 0.05, "sig", "not_sig"))
di <- DVGs %>% arrange(p_value) # Sort by pvalue in ascending order
ngene <- nrow(DVGs) # Get number of genes 
di <- di %>% mutate(observed=-log10(p_value), expected=-log10(ppoints(ngene))) # Get observed and expected p-values for qqplots

## 2.1.1 qq plots
print("Plotting qq plot for expected and observed qvalue")

fig1 <- ggplot(di, aes(x=expected,y=observed)) +
  ggrastr::rasterise(geom_point(size=0.3, colour="grey30"), dpi=300)+
  geom_abline(colour="red")+
  xlab(bquote("Expected"~ -log[10]~"("~italic(plain(P))~")"))+
  ylab(bquote("Observed"~-log[10]~"("~italic(plain(P))~")"))+
  theme_bw()+
  theme(strip.text=element_text(size=12))

## 2.1.2 Histogram
print("Plotting histogram for beta coefficient")

fig2 <- ggplot(di, aes(x = beta)) +
  geom_histogram(fill = "lightblue", colour = "grey30", binwidth = 0.5, show.legend = FALSE) +
  geom_vline(xintercept = c(-0.5, 0.5), color = "red", linetype = "dashed", alpha = 0.5, size = 0.5) +
  xlab(bquote("Beta coefficient")) +
  ylab("Counts") +
  theme_bw() +
  theme(strip.text = element_text(size = 12))

figfn <- paste0(outdir, "Figure1_Summary_Results_test.pdf")
pdf(figfn, width=12, height=6)
print(plot_grid(fig1, fig2, ncol=2, labels="AUTO", label_fontface="plain"))
dev.off()

print(paste("Summary plots for", cell_type, "saved in", figfn))
