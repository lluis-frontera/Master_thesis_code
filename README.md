#  Variability Analysis Pipeline on Single-cell RNA-seq data

This repository contains scripts and workflows for analyzing variability in large transcriptomics datasets such as the OneK1K cohort. The main components of the analysis are:

1. **Scran Analysis**
3. **Pseudobulk Analysis**
4. **Batch Analysis**
5. **Expression/Dispersion Boxplots**
   
## Repository structure

- `Master_thesis_code/age_thresholds_plots.R`
- `Master_thesis_code/batch_analysis/boxplots_dispersion_genes.R`
- `Master_thesis_code/batch_analysis/hexbin_dispersion_allgenes.R`
- `Master_thesis_code/expression_dispersion_boxplots.R`
- `Master_thesis_code/pseudobulk/1.dispersion_preprocess.R`
- `Master_thesis_code/pseudobulk/2.dispersion_analysis.R`
- `Master_thesis_code/pseudobulk/2.2.split_batches.R`
- `Master_thesis_code/pseudobulk/pseudobulk_plots.R`
- `Master_thesis_code/pseudobulk/permutations/1.dispersion_preprocess_permutations.R`
- `Master_thesis_code/pseudobulk/permutations/2.dispersion_analysis_permutations.R`
- `Master_thesis_code/pseudobulk/permutations/pseudobulk_permutations_labels.R`
- `Master_thesis_code/scran/scran_main_analysis.R`
- `Master_thesis_code/scran/scran_merging_celltypes.R`
- `Master_thesis_code/scran/scran_permutations.R`
- `Master_thesis_code/scran/scran_permutations_labels.R`
- `Master_thesis_code/scran/scran_plots.R`

### root directory

  - `Master_thesis_code/age_thresholds_plots.R`

     This script generates the barplot related to the number of differntially variable genes (DVG) detected in CD4 Naive using different age thresholds.

  - `Master_thesis_code/expression_dispersion_boxplots.R`

     This script is used to create boxplots showing the expression and dispersion data for top DVGs per cell type using both single-cell and pseudobulk analyses.

### batch_analysis directory

  - `Master_thesis_code/batch_analysis/boxplots_dispersion_genes.R`

     This script generates boxplots to show the variability of gene expression data across different batches compared with the total one for top DVGs separating by expression bins (low, mid and high).

  - `Master_thesis_code/batch_analysis/hexbin_dispersion_allgenes.R`

     This script creates hexbin plots to visualize the difference in total vs batches dispersion for all genes according to their expression levels.

### pseudobulk directory

  - `Master_thesis_code/pseudobulk/1.dispersion_preprocess.R`

     This script is used to obtain a raw variability value per gene per donor for each cell type after aggregating single-cell RNA-seq to the donor level. It is the first step in the pseudobulk analysis.
    
  - `Master_thesis_code/pseudobulk/2.dispersion_analysis.R`

     This script is used to detect the DVGs by applying linear models for each cell type. It represents the main dispersion analysis on the pseudobulk pipeline.

  -  `Master_thesis_code/pseudobulk/2.2.split_batches.R`

     This script handles the splitting of data into batches for further analysis for each cell type.
  
  -  `Master_thesis_code/pseudobulk/pseudobulk_plots.R`

     This script generates various plots for visualizing pseudobulk data. It includes the plotting for the permutation analysis, the barplot for the number of DVGs per cell type and the variability decomposition. 

#### permutations subdirectory

   - `Master_thesis_code/pseudobulk/permutations/1.dispersion_preprocess_permutations.R`

     This script is used for the preprocessing step for permutation analyses in pseudobulk data.

   - `Master_thesis_code/pseudobulk/permutations/2.dispersion_analysis_permutations.R`

     This script performs the main variability analysis using permutation tests.

   - `Master_thesis_code/pseudobulk/permutations/pseudobulk_permutations_labels.R`

     This script perform label shuffling for the permutation tests.

### scran directory

   - `Master_thesis_code/scran/scran_main_analysis.R`

     This script performs the main analysis using the _scran_ package, which is commonly used for single-cell RNA-seq data analysis. The analysis could include normalization, variance modeling, and identification of DVGs.

   - `Master_thesis_code/scran/scran_merging_celltypes.R`

     This script is used to merge _scran_ results from different cell types. It involves aggregating results across cell types for a comprehensive analysis.

   - `Master_thesis_code/scran/scran_permutations.R`

     This script probably performs permutation tests using the scran package. It might be used to assess the significance of observed patterns in single-cell RNA-seq data by comparing them to random permutations.

   - `Master_thesis_code/scran/scran_permutations_labels.R`

     This script likely handles labeling or organizing the results of permutation tests conducted with scran. It might involve annotating permutation test results for easier interpretation and visualization.

   - `Master_thesis_code/scran/scran_plots.R`

     This script likely generates various plots using the scran package. These plots might include visualizations of gene expression, dispersion, cell type distributions, and other relevant metrics from single-cell RNA-seq data.
