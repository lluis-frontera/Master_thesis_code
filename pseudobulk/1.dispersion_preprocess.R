#!/usr/bin/env Rscript

# To fix the Matrix package in Seurat objects
Csparse_validate = 'CsparseMatrix_validate'

########## SETTING LIBRARIES ############### 

shhh <- suppressPackageStartupMessages
shhh(library(corpcor))
shhh(library(tidyr))
shhh(library(Matrix))
shhh(library(MASS))
shhh(library(scales))
shhh(library(dplyr))
shhh(library(parallel))
shhh(library(data.table))
shhh(library(future))
shhh(library(purrr))
shhh(library(Rcpp))
shhh(library(reshape2))
shhh(library(qqman))
shhh(library(Seurat))
shhh(library(ggplot2))
shhh(library(cowplot))
shhh(library(grid))
shhh(library(gridExtra))
shhh(library(RColorBrewer))
shhh(library(gtable))
shhh(library(pheatmap))
shhh(library(corrplot))
shhh(library(viridis))
shhh(library(scales))
shhh(library(optparse))
shhh(library(SingleCellExperiment))
shhh(library(stringi))

############################## OPTIONS PARSER ###################################### 

option_list = list(
  make_option(c("--cell_level"), action="store", default=NA, type='character',
              help="From Azimuth: low resolution (predicted.celltype.l1) or high resolution (cell_type)"),
  make_option(c("--cell_type"), action="store", default=NA, type='character',
              help="Cell types in low resolution (predicted.celltype.l1) or high resolution (cell_type)"))
opt = parse_args(OptionParser(option_list=option_list))

options(digits = 4) # Modify number of digits printed

########## SETTING DIRECTORIES ############### 

path_cluster <- "/gpfs/projects/bsc83/Projects/scRNAseq/lluisf"
path_local <- "/home/lfronter/Documents"

if(startsWith(path_local, getwd())){
  main.dir <- path_local
} else {
  main.dir <- path_cluster
}

outdir <- paste0(main.dir, "/03_Results/Donor_variability/", opt$cell_level, "/", opt$cell_type, "/")

if(!file.exists(outdir)){dir.create(outdir, showWarnings=F, recursive = T)}

####################### Gene expression variance and differential analysis using OneK1K pseudobulk ##############################
                    ######## Based on script by Julong wei, modified by Lluis Frontera ##########

### SECTION 1. Calculate gene variance, dispersion and mean based on negative binomial distribution ### 

######################
### 1. Read data ###
######################

## 1.1 Read data
print("1.1 Read data from Single Cell Experiment object")

sce_object <- readRDS(paste0("/gpfs/projects/bsc83/Data/scRNAseq/Yazar2022/sce_data_objects/", opt$cell_type, "_", opt$cell_level, "_sceraw.rds"))

### 1.2 Getting cell information and using to classify

print("1.2 Getting cell information")

## Filter donors with at least 5 cells
don_distr <- colData(sce_object) %>% as.data.frame() %>% group_by(assignment) %>% summarise(ncell=n()) # Get number of cells for each donor
bti <- don_distr %>% filter(ncell>20) %>% dplyr::select(assignment) %>% unlist() # Filter by ncells: 15 at 5%, 25 at 10%, default 20
donor_diff <- length(don_distr$assignment) - length(bti)

## Get sce object with filtered donors
sce_object <- sce_object[, sce_object$assignment %in% bti]

#################
### FILTERING ###
#################

## Filter genes with less than 20 reads across cells
gene_counts <- Matrix::rowSums(SingleCellExperiment::counts(sce_object))
selected_genes <- names(gene_counts[gene_counts >= 20])
sce_filtered <- sce_object[selected_genes, ]

## Get counts and metadata for the filtered dataframe (gene and cell filtering)
meta <- colData(sce_filtered) # Get metadata for filtered cells
count <- SingleCellExperiment::counts(sce_filtered) # Get counts for filtered cells and selected genes

# Report
cat("\n")
print("################### Report ###################")
print(paste0('Cell level: ', opt$cell_level))
print(paste0('Cell type: ', opt$cell_type))
print(paste0('Output dir: ', outdir))
print(paste("Donors excluded due to low number of cells (< 20):", as.numeric(donor_diff)))
print(paste("Donors included in the following analysis:", length(bti)))
print(paste("Cells included in the following analysis:", nrow(meta)))
print(paste("Genes with more than 20 reads across cells (for all dataset):", nrow(count)))
print("################### End report ###################")
cat("\n")

##############################################
### 2. Estimate gene expression variance ###
##############################################

### 2.1 Estimate gene expression variance based on Negative Binomial model 

## 2.1.1 Calculate for each combination, from the same individual and the same cell type

print("2.1.1 Estimate gene expression variance based on Negative Binomial model for each combination from the same individual and the same cell type")

size_all <-  median(colSums(count)) # Get median for expression counts of genes for current donor (for normalization)

TMP <- lapply(1:length(bti), function(i){ # WITH SECOND FILTERING

   time0 <- Sys.time()
   bti0 <- as.character(bti[i]) # Get current donor
   cell_i <- as.character(meta[meta$assignment==bti0, "bare_barcode_lane"]) # Get cell names for current donor
   size_all <-  median(colSums(count[, cell_i])) # Get median for expression counts of genes for current donor (for normalization)

   # Apply filtering to get only genes in which there are more than 5 counts in at least 5 cells

   count_i <- count[, cell_i] # Extract gene expression data for current donor
   genes <- rownames(count_i) # Get all available genes
   row_sum <- rowSums(count_i) # Sum of counts for each gene
   num_cells_expressing_gene <- rowSums(count_i > 0) # Get number of cells expressing certain gene
   genes <- names(which(row_sum > 5 & num_cells_expressing_gene > 5)) # Get only genes in which there are more than 5 counts in at least 5 cells
   count_i <- count_i[rownames(count_i) %in% genes, ] # Retrieve expression information for selected genes
   cell_size <- colSums(count_i)/size_all # Normalizing

   print(paste("Genes with more than 5 reads in at least 5 cells in donor", i, "(second filtering):", length(genes)))

   genes <- rownames(count_i) # Get all genes
   v <- rep(NA,nrow(count_i)) # Variance parameter
   names(v) <- genes # Assign gene_id as row names
   b <- rep(NA,nrow(count_i)) # Mean parameter
   names(b) <- genes # Assign gene_id as row names
   phi <- rep(NA,nrow(count_i)) # Dispersion parameter
   names(phi) <- genes # Assign gene_id as row names
   se.mu <- rep(NA,nrow(count_i)) # Standard error mean
   names(se.mu) <- genes # Assign gene_id as row names
   se.phi <- rep(NA,nrow(count_i)) # Standard error dispersion
   names(se.phi) <- genes # Assign gene_id as row names

   # Function to calculate negative log-likelihood for negative binomial distribution
   
   nb_llik <- function(theta, x, size){
     ## Negative binomial log likelihood
     mu <- size*exp(theta[1])
     inv_disp <- exp(theta[2])

     res <- x*log(mu/inv_disp)- x*log(1+mu/inv_disp)-
       inv_disp*log(1+mu/inv_disp)+
       lgamma(x+inv_disp)-lgamma(x+1)-lgamma(inv_disp)
     llik <- -sum(res)
   }
   
   genes_converged <- 0

   ## For each gene
   tmp <- mclapply(genes, function(k){ # Put "gene0" if filtering is needed
      xk <- as.numeric(count_i[k, ])

      # Given initial value by mean, var
      mu_0 <- mean(xk/cell_size) # Mean
      va_0 <- var(xk/cell_size) # Variance
      phi_0 <- va_0/mu_0^2 # Dispersion
      theta <- c(log(mu_0),log(1/phi_0)) # Parameter vector

      parm <- try(optim(par=theta, fn=nb_llik, x=xk, size=cell_size,
                        method="L-BFGS-B", hessian=T,
                        lower=c(-Inf,-Inf),upper=c(Inf, Inf)), silent=T) # Execute parameter optimization
    
      # if-1, "try-error"
      if (class(parm)!="try-error"){

         # if-4, it is convergence
         if (parm$convergence==0){
           
           genes_converged <<- genes_converged + 1
            
            parr <- parm$par
            mu <- exp(parr[1])
            phk <- 1/exp(parr[2])
            bk <- mu
            vk <- mu^2*phk
            sek <- sqrt(diag(corpcor::pseudoinverse(parm$hessian)))
            res <- c(vk, bk, phk, sek[1], sek[2])
            if (parm$value<=0) res <- c(NA,NA,NA,NA,NA)
         } else {
            res <- c(NA,NA,NA,NA,NA)
         }
      } else {
         res <- c(NA,NA,NA,NA,NA)
      }
      return(res)
   } , mc.cores=1)
   
   tmp <- do.call(rbind, tmp)

   # Variance
   v1 <- tmp[, 1]
   v[genes] <- v1
   # Mean value
   b1 <- tmp[, 2]
   b[genes] <- b1
   # Dispersion
   phi1 <- tmp[, 3]
   phi[genes] <- phi1
   # Standard error
   se.mu[genes] <- tmp[, 4]
   se.phi[genes] <- tmp[, 5]

   time1 <- Sys.time()
   elapsed <- difftime(time1, time0, units="secs")
   return(list(v, b, phi, se.mu, se.phi, length(genes), genes_converged))
})

### 2.2 Save results for Negative Binomial model

print("2.2.2. Save parameter values in RDS files")

# Variance
v <- lapply(TMP, function(ii) ii[[1]]) # Extract variance of each gene of each donor
Vx <- do.call(cbind, v) # Combine into a common dataframe
colnames(Vx) <- as.character(bti) # Set column names based on batch_donor
var.fn <- paste0(outdir, "Variance_5threshold.rds")
saveRDS(Vx, file = var.fn)
print(paste("Variance values for each gene and donor saved in", var.fn))

# Mean
b <- lapply(TMP, function(ii) ii[[2]]) # Extract mean of each gene of each donor
Bx <- do.call(cbind, b) # Combine into a common dataframe
colnames(Bx) <- as.character(bti) # Set column names based on batch_donor
mean.fn <- paste0(outdir, "Mean_5threshold.rds")
saveRDS(Bx, file = mean.fn)
print(paste("Mean values for each gene and donor saved in", mean.fn))

# Dispersion
phx <- lapply(TMP, function(ii) ii[[3]]) # Extract dispersion of each gene of each donor
Phx <- do.call(cbind, phx) # Combine into a common dataframe
colnames(Phx) <- as.character(bti) # Set column names based on batch_donor
dis.fn <- paste0(outdir, "Dispersion_5threshold.rds")
saveRDS(Phx, file = dis.fn)
print(paste("Dispersion values for each gene and donor saved in", dis.fn))

# Standard error mean
Sx.mu <- lapply(TMP, function(ii) ii[[4]]) # Extract standard error mean of each gene of each donor
Sx.mu <- do.call(cbind, Sx.mu) # Combine into a common dataframe
colnames(Sx.mu) <- as.character(bti) # Set column names based on batch_donor
se_mean.fn <- paste0(outdir, "SE_Mean_5threshold.rds")
saveRDS(Sx.mu, file = se_mean.fn)
print(paste("Standard error of mean values for each gene and donor saved in", se_mean.fn))

# Standard error dispersion
Sx.phi <- lapply(TMP, function(ii) ii[[5]]) # Extract standard error dispersion of each gene of each donor
Sx.phi <- do.call(cbind, Sx.phi) # Combine into a common dataframe
colnames(Sx.phi) <- as.character(bti) # Set column names based on batch_donor
se_dis.fn <- paste0(outdir, "SE_Dispersion_5threshold.rds")
saveRDS(Sx.phi, file = se_dis.fn)
print(paste("Standard error of dispersion values for each gene and donor saved in", se_dis.fn))

# Number genes
ngene <- lapply(TMP, function(ii) ii[[6]]) # Extract number of genes analyzed
ngene <- do.call(c,ngene) # Combine into a common vector
names(ngene) <- as.character(bti) # Set column names based on batch_donor
ngene.fn <- paste0(outdir, "Number_genes_5threshold.rds")
saveRDS(ngene, file = ngene.fn)
print(paste("Number of tested genes for each donor saved in", ngene.fn))

# Number of converged genes
ngene_converged <- lapply(TMP, function(ii) ii[[7]]) # Extract number of genes analyzed
ngene_converged <- do.call(c,ngene_converged) # Combine into a common vector
names(ngene_converged) <- as.character(bti) # Set column names based on batch_donor
ngene_converged.fn <- paste0(outdir, "Number_genes_converged_5threshold.rds")
saveRDS(ngene_converged, file = ngene_converged.fn)
print(paste("Number of tested genes for each donor saved in", ngene_converged.fn))

### 2.3 Calculate residual dispersion

print("2.3. Calculate residual dispersion")

# Load negative binomial parameters 
Vx <- readRDS(paste0(outdir, "Variance_5threshold.rds"))
Bx <- readRDS(paste0(outdir, "Mean_5threshold.rds"))
Phx <- readRDS(paste0(outdir, "Dispersion_5threshold.rds"))
Sx.mu <-readRDS(paste0(outdir, "SE_Mean_5threshold.rds"))
Sx.phi <-readRDS(paste0(outdir, "SE_Dispersion_5threshold.rds"))

### 2.3.1 Regression

print("2.3.1 Calculate regression")

d1 <- melt(Vx)
d2 <- melt(Bx)
d3 <- melt(Phx)
d4 <- melt(Sx.mu)
d5 <- melt(Sx.phi)

ddx <- data.frame(genes=d1$Var1, donor=d1$Var2, va=d1$value, mu=d2$value,  # Var1 = gene_id, Var2 = donor
                  phi=d3$value, se.mu=d4$value, se.phi=d5$value) %>% 
                  drop_na(va, mu, phi, se.mu, se.phi)

dd3 <- ddx %>% mutate(x=log2(mu), y=log2(phi)) # Make log2 of mean and dispersion

lmx <- dd3 %>% group_by(donor) %>% # Group by donor group
       nest() %>% # Each group represented as list column
       mutate(lmr=map(data, ~lm(y~x, data=.x))) # Execute linear regression

lmr <- lmx$lmr
names(lmr) <- lmx$donor

### 2.3.2 Residual Phx

print("2.3.2 Calculate residual dispersion")

Phx_log <- log2(Phx)
Bx_log <- log2(Bx)

PhxNew <- Phx_log # Clone log dispersion for analysis

for (i in bti){ # For each donor...
  Phx0 <- Phx_log[,i] # Get dispersions
  Bx0 <- Bx_log[,i] # Get means

  lm0 <- coef(lmr[[i]]) # Fit a linear model
  a <- lm0[1] # Extract coefficients
  b1 <- lm0[2] # Extract coefficients
  PhxNew[,i] <- Phx0-(a+b1*Bx0) # Substract modeled effect from original dispersion values
}

PhxNew2 <- 2^PhxNew # Back-transform dispersion file

saveRDS(PhxNew2, file = paste0(outdir, "Dispersion_updated_5threshold.rds"))
print(paste("Residual dispersion saved in", paste0(outdir, "Dispersion_updated_5threshold.rds")))

### Additional filtering: filter out genes with dispersion values in less than 50% of donors 

threshold_value <- round(length(bti)/2, 0) # Get filtering threshold
PhxNew2_filtered <- PhxNew2[rowSums(!is.na(PhxNew2)) >= threshold_value, ] # Recover genes having less than 50% NAs

saveRDS(PhxNew2_filtered, file = paste0(outdir, "Dispersion_updated_5threshold_filtered.rds"))
print(paste("Residual dispersion filtered saved in", paste0(outdir, "Dispersion_updated_5threshold_filtered.rds")))

##############################
### 3. Summary parameters ###
##############################

# Create dataframe for plotting
Vx <- readRDS(paste0(outdir, "Variance.rds"))
Bx <-readRDS(paste0(outdir, "Mean.rds"))
Phx <-readRDS(paste0(outdir, "Dispersion.rds"))
PhxNew <-readRDS(paste0(outdir, "Dispersion_updated.rds"))
Sx.mu <-readRDS(paste0(outdir, "SE_Mean.rds"))
Sx.phi <-readRDS(paste0(outdir, "SE_Dispersion.rds"))

d1 <- melt(Vx)
d2 <- melt(Bx)
d3.1 <- melt(Phx)
d3.2 <- melt(PhxNew)
d4 <- melt(Sx.mu)
d5 <- melt(Sx.phi)

ddx <- data.frame(genes=d1$Var1, donor=d1$Var2, va=d1$value, mu=d2$value,  # Var1 = gene_id, Var2 = donor
                  phi=d3.1$value, phiNew=d3.2$value, se.mu=d4$value, se.phi=d5$value) %>%
       drop_na(va, mu, phi, se.mu, se.phi) # Remove NA values

############# FUNCTIONS DEFINITION ###################

# Defining function to extract R-squared value from linear model

fmod <- function(df){
  lm0 <- lm(y~x, data=df)
  r2 <- round(summary(lm0)$r.squared, digits=3)
  eq <- bquote(italic(R)^2==.(r2))
  as.character(as.expression(eq))
}

# Defining function to extract formatted equations indicating significance levels based on p-value       

feq <- function(x){
  r <- round(as.numeric(x$estimate),digits=3)
  p <- x$p.value
  if(p<0.001) symb <- "***"
  if(p>=0.001 & p<0.01) symb <- "**"
  if (p>=0.01 & p<0.05) symb <- "*"
  if(p>0.05) symb <- "NS"
  
  eq <- bquote(italic(R)==.(r)~","~.(symb))
  as.character(as.expression(eq)) 
} 

# Defining function to extract the numeric estimate without the formatted equation

feq2 <- function(x){
  r <- round(as.numeric(x$estimate),digits=3)
  p <- x$p.value
  if(p<0.001) symb <- "***"
  if(p>=0.001 & p<0.01) symb <- "**"
  if (p>=0.01 & p<0.05) symb <- "*"
  if(p>0.05) symb <- "NS"
  r 
} 

### 3.1 Density plot for variance, mean and dispersion

print("3.1 Calculate density plots for variance, mean and dispersion")

# Variance at 10% percentage
p1 <- ggplot(ddx, aes(x=log10(va+1e-04))) +
     geom_density() + 
     xlim(-5,5) + 
     scale_y_continuous(labels=label_number(accuracy = 0.01)) +
     xlab("log10(Var+1e-04)") +
     ylab(NULL) +
     theme_bw() + 
     ggtitle("Variance") + 
     theme(plot.title=element_text(hjust=0.5))

# Mean 
p2 <- ggplot(ddx, aes(x=log10(mu))) + 
      xlim(-5,5) + 
      scale_y_continuous(labels=label_number(accuracy = 0.01)) +
      geom_density()+
      xlab("log10(Mean)") +
      ylab(NULL) +
      theme_bw() + 
      ggtitle("Mean") + 
      theme(plot.title=element_text(hjust=0.5))
        
# Dispersion at 10%  
p3 <- ggplot(ddx, aes(x=log10(phi+1e-02))) +
     geom_density() + 
     scale_y_continuous(labels=label_number(accuracy = 0.01)) +
     xlim(-5,5) +
     xlab("log10(Dis+0.01)") +
     ylab(NULL) +
     theme_bw() + 
     ggtitle("Dispersion") + 
     theme(plot.title=element_text(hjust=0.5))

figfn <- paste0(outdir, "Summary_parameters/")
if(!file.exists(figfn)){dir.create(figfn, showWarnings=F)}

figfn_1 <- paste0(figfn, "1.1.Density_basic_par.png")
png(figfn_1, width=1000, height=500, res=120)
print(plot_grid(p1, p2, p3, ncol=3))
invisible(dev.off())

### 3.2 Distribution of standard error of mean and dispersion

print("3.2 Calculate distribution plots for standard error of mean and dispersion")

col_num <- length(levels(ddx$donor))
col_per_iteration <- 6
num_iterations <- ceiling(col_num / col_per_iteration) # Calculate the number of iterations needed

# Use lapply to iterate over subsets of columns
lapply(seq(1, num_iterations), function(i) {
  start_col <- (i - 1) * col_per_iteration + 1
  end_col <- min(i * col_per_iteration, col_num)
  
  # Extract the columns for the current iteration
  group_donors <- ddx %>% filter(donor %in% levels(donor)[start_col:end_col])
  
  # Standard error of mean
  fig1 <- ggplot(group_donors, aes(x=se.mu)) +
    geom_histogram(fill="grey70", color="grey40", bins=50) +
    facet_wrap(~donor, nrow=2, scales="free") +
    ylab(NULL) +
    xlab("SE Mean") +
    theme_bw()
  
  figfn_2 <- paste0(figfn, "SE_mean_distr/")
  if(!file.exists(figfn_2)){dir.create(figfn_2, showWarnings=F)}
  figfn_2_donors <- paste0(figfn_2, "1.2.Density_SE_mean_", start_col, "_", end_col, ".png")
  png(figfn_2_donors, width=600, height=500, res=120)
  print(fig1)
  invisible(dev.off())
  
  # Standard error of dispersion
  fig2 <- ggplot(group_donors, aes(x=se.phi))+
    geom_histogram(fill="grey70", color="grey40", bins=50)+
    facet_wrap(~donor, nrow=2, scales="free") +
    ylab(NULL) +
    xlab("SE Dispersion") +
    theme_bw()
  
  figfn_3 <- paste0(figfn, "SE_dispersion_distr/")
  if(!file.exists(figfn_3)){dir.create(figfn_3, showWarnings=F)}
  figfn_3_donors <- paste0(figfn_3, "1.2.Density_SE_dispersion_", start_col, "_", end_col, ".png")
  png(figfn_3_donors, width=600, height=500, res=120)
  print(fig2)
  invisible(dev.off())
  
  ## Standard error of dispersion (filtered)
  
  fig3 <- ggplot(group_donors %>% filter(se.phi<6.45), aes(x=se.phi))+
    geom_histogram(fill="grey70", color="grey40", bins=50)+
    facet_wrap(~donor, nrow=2, scales="free") +
    ylab(NULL) +
    xlab("SE Dispersion (filt)") +
    theme_bw()
  
  figfn_4 <- paste0(figfn, "SE_dispersion_filtered_distr/")
  if(!file.exists(figfn_4)){dir.create(figfn_4, showWarnings=F)}
  figfn_4_donors <- paste0(figfn_4, "1.2.Density_SE_dispersion_filt_", start_col, "_", end_col, ".png")
  png(figfn_4_donors, width=600, height=500, res=120)
  print(fig3)
  invisible(dev.off())

  ### 3.3 Scatter plots showing relations between mean variance and mean dispersion
               
  # 3.3.1 Scatter plot of correlation between mean vs variance
  
  print("3.3.1 Scatter plot for mean vs variance for each donor")
    
  dd3 <- group_donors %>% mutate(x=log10(mu), y=log10(va))  
  
  anno_df1 <- dd3 %>%
              group_by(donor) %>%
              nest() %>%
              mutate(corr=map(data, ~cor.test((.x)$x, (.x)$y, method="pearson")),
                     corr2=map(corr,feq),
                     r=map_dbl(corr, feq2))
  
  x <- anno_df1 %>% dplyr::select(-data, -corr, -corr2)
  
  fig1 <- ggplot(dd3, aes(x=x,y=y))+
          geom_point() +
          geom_text(data=anno_df1, aes(x=2, y=6, label=corr2), parse=T, size=2)+  ##CPM 2,9
          scale_fill_viridis_c()+
          xlab(bquote(log[10]~"("~mu~")"))+
          ylab(bquote(log[10]~"(va)"))+          
          facet_wrap(~donor)+
          geom_smooth(method="lm",formula=y~x, size=0.5)+
          theme_bw()+
          theme(legend.title=element_blank(),
                legend.text=element_text(size=7),
                legend.key.size=grid::unit(0.5,"lines"),
                strip.text=element_text(size=10))
  
  figfn_5 <- paste0(figfn, "Correlation_mean_var/")
  if(!file.exists(figfn_5)){dir.create(figfn_5, showWarnings=F)}
  figfn_5_donors <- paste0(figfn_5, "2.1.Mean_var_scatter_", start_col, "_", end_col, ".png")
  png(figfn_5_donors, width=900, height=800, res=150)
  print(fig1)
  invisible(dev.off())

  # 3.3.2 Scatter plot of the correlation between mean and dispersion
  
  print("3.3.2 Scatter plot for mean vs dispersion for each donor")
  
  dd3 <- group_donors %>% mutate(x=log10(mu), y=log10(phi))  
  
  anno_df2 <- dd3 %>%
              group_by(donor) %>%
              nest( )%>%
              mutate(corr=map(data, ~cor.test((.x)$x, (.x)$y, method="pearson")),
                     corr2=map(corr,feq),
                     r=map_dbl(corr, feq2))

  x <- anno_df2 %>% dplyr::select(-corr,-corr2)
  
  fig2 <- ggplot(dd3, aes(x=x,y=y))+
    geom_point() +
    geom_text(data=anno_df2, aes(x=2, y=2, label=corr2), parse=T, size=2)+
    scale_fill_viridis_c()+
    xlab(bquote(log[10]~"("~mu~")"))+
    scale_y_continuous(bquote(log[10]~"("~phi~")"), expand=expansion(mult=0.2))+
    facet_wrap(~donor)+
    geom_smooth(method="lm",formula=y~x, size=0.5)+
    theme_bw()+
    theme(legend.title=element_blank(),
          legend.text=element_text(size=7),
          legend.key.size=grid::unit(0.5,"lines"),
          strip.text=element_text(size=10))
  
  figfn_6 <- paste0(figfn, "Correlation_mean_dispersion/")
  if(!file.exists(figfn_6)){dir.create(figfn_6, showWarnings=F)}
  figfn_6_donors <- paste0(figfn_6, "2.2.Mean_dispersion_scatter_", start_col, "_", end_col, ".png")
  png(figfn_6_donors, width=900, height=800, res=150)
  print(fig2)
  invisible(dev.off())
  
  # 3.3.3 Scatter plots of the correlation between mean vs new dispersion
  
  print("3.3.3 Scatter plot for mean vs new dispersion for each donor")
  
  dd3 <- group_donors %>% mutate(x=log10(mu),y=log10(phiNew))
  
  anno_df3 <- dd3 %>%
              group_by(donor) %>%
              nest() %>%
              mutate(corr=map(data, ~cor.test((.x)$x, (.x)$y, method="pearson")),
                     corr2=map(corr, feq),
                     r=map_dbl(corr, feq2))
                     
  x <- anno_df3 %>% dplyr::select(-corr,-corr2)
    
  fig3 <- ggplot(dd3, aes(x=x,y=y))+
          geom_point() +
          geom_text(data=anno_df3, aes(x=1, y=2, label=corr2), parse=T, size=2)+ 
          scale_fill_viridis_c()+
          xlab(bquote(log[10]~"("~mu~")"))+
          scale_y_continuous(bquote(log[10]~"("~phi~")"), expand=expansion(mult=0.2))+
          facet_wrap(~donor)+          
          geom_smooth(method="lm",formula=y~x, size=0.5)+
          theme_bw()+
          theme(legend.title=element_blank(),
                legend.text=element_text(size=7),
                legend.key.size=grid::unit(0.5,"lines"),
                strip.text=element_text(size=10))
  
  figfn_7 <- paste0(figfn, "Correlation_mean_newdispersion/")
  if(!file.exists(figfn_7)){dir.create(figfn_7, showWarnings=F)}
  figfn_7_donors <- paste0(figfn_7, "2.3.Mean_newdispersion_scatter_", start_col, "_", end_col, ".png")
  png(figfn_7_donors, width=900, height=800, res=150)
  print(fig3)
  invisible(dev.off())

})