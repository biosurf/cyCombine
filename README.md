
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

### 1\. Create virtual environment with renv (optional)

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

### 2\. Install dependencies

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
devtools::load_all("path/to/cyCombine")
```

## Usage

### From a directory of preprocessed .fcs files

``` r
library(cyCombine)
library(magrittr)
# Directory containing .fcs files
data_dir <- "data/raw"
# Markers of interest
markers <- c("CD20", "CD3", "CD27", "CD45RA", "CD279", "CD5", "CD19", "CD14", "CD45RO", "GranzymeA", "GranzymeK", "FCRL6", "CD355", "CD152", "CD69", "CD33", "CD4", "CD337", "CD8", "CD197", "LAG3", "CD56", "CD137", "CD161", "FoxP3", "CD80", "CD270", "CD275", "CD134", "CD278", "CD127", "KLRG1", "CD25", "HLADR", "TBet", "XCL1")

# Compile fcs files, down-sample, and preprocess
preprocessed <- prepare_data(data_dir = data_dir,
                             markers = markers,
                             metadata = paste0(data_dir, "/metadata.xlsx"), # Can also be .csv file or data.frame object
                             sample_ids = NULL,
                             batch_ids = "Batch",
                             filename_col = "FCS_name",
                             condition = "Set",
                             down_sample = TRUE,
                             sample_size = 500000,
                             seed = 473,
                             cofactor = 5) 
saveRDS(preprocessed, file = "_data/cycombine_raw_preprocessed.RDS")

# Run batch correction
corrected <- preprocessed %>%
  batch_correct(markers = markers,
                norm_method = "scale",
                som_type = "fsom", # Anything else to use Kohonen clustering method
                covar = "condition")
saveRDS(corrected, file = "_data/cycombine_raw_corrected.RDS")
```

### The modular alternative

``` r
library(cyCombine)
library(magrittr)
# Direcory containing .fcs files
data_dir <- "data/raw"
# Markers of interest
markers <- c("CD20", "CD3", "CD27", "CD45RA", "CD279", "CD5", "CD19", "CD14", "CD45RO", "GranzymeA", "GranzymeK", "FCRL6", "CD355", "CD152", "CD69", "CD33", "CD4", "CD337", "CD8", "CD197", "LAG3", "CD56", "CD137", "CD161", "FoxP3", "CD80", "CD270", "CD275", "CD134", "CD278", "CD127", "KLRG1", "CD25", "HLADR", "TBet", "XCL1")

# Compile fcs files, down-sample, and preprocess
flowset <- compile_fcs(data_dir = data_dir,
                   pattern = "\\.fcs")

df <- convert_flowset(metadata = paste0(data_dir, "/metadata.xlsx"),
                      sample_ids = NULL,
                      batch_ids = "Batch",
                      filename_col = "FCS_name",
                      condition = "Set",
                      down_sample = TRUE,
                      sample_size = 500000,
                      seed = 473)

preprocessed <- df %>% 
  transform_asinh(markers = markers)

saveRDS(preprocessed, file = "_data/cycombine_raw_preprocessed.RDS")

# Run batch correction
labels <- preprocessed %>%
  normalize(markers = markers,
            norm_method = "scale") %>%
  create_fsom(markers = markers) # Alternatively, use create_som() for Kohonen clustering method

corrected <- preprocessed %>%
  correct_data(label = labels,
               covar = "condition")
saveRDS(corrected, file = "_data/cycombine_raw_corrected.RDS")
```

<!-- ### From a flowset -->

<!-- ```{r, eval = FALSE} -->

<!-- library(cyCombine) -->

<!-- library(magrittr) -->

<!-- # Load data -->

<!-- # Should contain the flowset, sample_ids, batch_ids, and markers of interest -->

<!-- load("data/flowset.Rdata") -->

<!-- # Convert flowset to workable datafram and transform data -->

<!-- preprocessed <- flowset %>% -->

<!--   convert_flowset(batch_ids = batch_ids, -->

<!--                   sample_ids = sample_ids, -->

<!--                   down_sample = TRUE, -->

<!--                   sample_size = 100000, -->

<!--                   seed = 473) %>%  -->

<!--   transform_asinh(markers = markers) -->

<!-- # Run batch correction -->

<!-- corrected <- preprocessed %>% -->

<!--   batch_correct(seed = 473) -->

<!-- ``` -->

## Plotting

``` r
# Full analysis can be run with - type ?run_analysis to see how you can modify the analysis
run_analysis(tool = "cycombine", data = "raw", data_dir = "_data")

# Otherwise, plots can be made like so:
density_plots(uncorrected = preprocessed,
                corrected = corrected,
                markers = markers,
                filename = 'figs/densities_withcovar.png')

# PCA plot uncorrected
pca1 <- preprocessed %>%
  plot_dimred('uncorrected', type = 'pca')
  
# PCA plot corrected
pca2 <- corrected %>%
  plot_dimred('corrected', type = 'pca')
save_two_plots(pca1, pca2, filename = 'figs/pca.png')

# UMAP
# UMAP plot uncorrected
set.seed(473)
sample <- sample(1:nrow(uncorrected), 20000)
plot1 <- plot_dimred(uncorrected[sample,], type = 'umap', name = 'Uncorrected')
plot2 <- plot_dimred(corrected[sample,], type = 'umap', name = 'Corrected')
save_two_plots(plot1, plot2, filename = 'figs/umap.png')
```
