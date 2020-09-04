
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

## Example

``` r
library(cyCombine)
library(magrittr)
# Load input data
load("data/raw_data.Rdata")
# Preprocess
preprocessed_data <- raw_data %>% 
  preprocess()

# Batch correct
corrected_data <- preprocessed_data %>% 
  batch_correct()
```
