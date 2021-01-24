## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval=FALSE--------------------------------------------------------------
#  # Preprocessing step
#  library(cyCombine)
#  library(magrittr)
#  load("../data/raw/DFCI_panel1_data.Rdata")
#  # Run preprocessing
#  preprocessed_data <- panel1_data %>%
#    preprocess()
#  # Store result
#  save(preprocessed_data, file = "../data/01_preprocessed.Rdata")
#  
#  

## ---- eval=FALSE--------------------------------------------------------------
#  # Batch correction
#  load("../data/01_preprocessed.Rdata")
#  corrected_data <- preprocessed_data$data %>%
#    batch_correct(markers = preprocessed_data$markers)
#  

