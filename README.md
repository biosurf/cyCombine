
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cyCombine

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
<!-- badges: end -->

<!-- ## Clone github repository -->

<!-- ``` {sh, eval = FALSE} -->

<!-- # Run in terminal -->

<!-- git clone git@github.com:shdam/cyCombine.git -->

<!-- ``` -->

<!-- ## Restore renv library -->

<!-- ``` {r, eval = FALSE} -->

<!-- # Open project in Rstudio -->

<!-- # Install renv and restore library -->

<!-- install.packages("renv") -->

<!-- library(renv) -->

<!-- renv::restore() -->

<!-- ``` -->

## Install as package

### 1\. Create virtual environment with renv

Initialize a local R environment:

``` r
# Open project in Rstudio
# Install and initialize renv 
install.packages("renv")
library(renv)
renv::init()
```

### 2\. Install package from github

``` r
# To ensure Rstudio looks up BioConductor packages run:
setRepositories(ind = c(1:6, 8))
# Then install package with
devtools::install_github("shdam/cyCombine")
```

## Install as cloned repository

### 1\. Clone repository

``` sh
# In terminal at desired directory
git clone git@github.com:shdam/cyCombine.git
```

### 2\. Install dependecies

``` r
# Open cyCombine.Rproj in Rstudio
# Install renv
install.packages("renv")
library(renv)
# Restore from lock file
renv::restore()
```

### 3\. Load package

``` r
install.packages("devtools")
library(devtools)
devtools::load_all("~/Rprojects/cyCombine")
```

## Usage

### From a directory of .fcs files

``` r
library(cyCombine)
library(magrittr)
# Direcory containing .fcs files
data_dir <- "data/raw"
# Markers of interest
markers <- c("CD20", "CD3", "CD27", "CD45RA", "CD279", "CD5", "CD19", "CD14", "CD45RO", "GranzymeA", "GranzymeK", "FCRL6", "CD355", "CD152", "CD69", "CD33", "CD4", "CD337", "CD8", "CD197", "LAG3", "CD56", "CD137", "CD161", "FoxP3", "CD80", "CD270", "CD275", "CD134", "CD278", "CD127", "KLRG1", "CD25", "HLADR", "TBet", "XCL1")

# Compile fcs files, down-sample, and preprocess
fcs_preprocessed <- preprocess(data_dir = data_dir,
                         meta_filename = "CyTOF samples cohort.xlsx",
                         markers = markers,
                         down_sample = TRUE,
                         sample_size = 300000,
                         seed = 473,
                         cofactor = 5) 
  
# Run batch correction
fcs_corrected <- fcs_preprocessed %>%
  batch_correct(seed = 473)
```

### From a flowset

``` r
library(cyCombine)
library(magrittr)
# Load data
# Should contain the flowset, sample_ids, batch_ids, and markers of interest
load("data/flowset.Rdata")

# Convert flowset to workable datafram and transform data
fcs_preprocessed <- flowset %>%
  convert_flowset(batch_ids = batch_ids,
                  sample_ids = sample_ids,
                  down_sample = TRUE,
                  sample_size = 100000,
                  seed = 473) %>% 
  transform_asinh(markers = markers)
  
# Run batch correction
fcs_corrected <- fcs_preprocessed %>%
  batch_correct(seed = 473)
```

## Plotting

``` r
density_plots(uncorrected = fcs_preprocessed,
                corrected = fcs_corrected,
                markers = markers,
                filename = 'figs/densities_withcovar.png')

# PCA plot uncorrected
pca1 <- preprocessed_data %>%
  dimred_plot('uncorrected', type = 'pca')
  
# PCA plot corrected
pca2 <- corrected_data %>%
  dimred_plot('corrected', type = 'pca')
save_two_plots(pca1, pca2, filename = 'figs/pca.png')

# UMAP
# UMAP plot uncorrected
umap1 <- preprocessed_data %>%
  dimred_plot('uncorrected', type = 'umap')

# UMAP plot corrected
umap2 <- dimred_plot(corrected_data, 'corrected', type = 'umap')
save_two_plots(umap1, umap2, filename = 'figs/umap.png')
```
