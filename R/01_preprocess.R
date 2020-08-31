#### Data pre-processing ####

# Load data - the loaded contains objects "fcs_raw", "sample_ids", 
# "batch_ids", "all_markers", which are the raw flowSet, sample ids per row, batch 
# ids per row and all measured markers in the panel
preprocess <- function(fcs_raw){
  # De-randomize and transform data using asinh
  panel_fcs <- fcs_raw[[1]] %>% 
    flowCore::parameters() %>% 
    Biobase::pData()
  
  panel_fcs$desc <- gsub(' ', '', gsub('-', '', gsub("\\d+[A-Za-z]+_", "", panel_fcs$desc)))
  
  fcs <- flowCore::fsApply(fcs_raw, function(x, cofactor = 5){
    colnames(x) <- panel_fcs$desc
    expr <- Biobase::exprs(x)
    expr <- asinh(ceiling(expr[, all_markers]) / cofactor)
    exprs(x) <- expr
    x
  })
  return(fcs)
}
if(FALSE){
  # Load raw data
  load("data/raw/DFCI_panel1_data.Rdata")
  # Extract objects
  fcs_raw <- panel1_data$fcs_raw
  all_markers <- panel1_data$all_markers
  sample_ids <- panel1_data$sample_ids
  batch_ids <- panel1_data$batch_ids
  
  
  # Extract expression data
  combined_expr <- preprocess(fcs_raw) %>% 
    flowCore::fsApply(exprs)
}

create_sample <- function(combined_expr,
                          batch_ids,
                          sample_ids,
                          sample_size = 100000,
                          seed = 473){
  set.seed(seed)
  sample <- sample(1:nrow(combined_expr), 100000)
  combined_expr <- combined_expr[sample,] %>% 
    tibble::as_tibble() %>% 
    mutate(batch_ids = batch_ids[sample],
           sample_ids = sample_ids[sample])
  return(combined_expr)
}
# Plot of distributions - downsampled to avoid problems with memory
# set.seed(473)
# sample <- sample(1:nrow(combined_expr), 100000)
# combined_expr <- combined_expr[sample,] %>% 
#   tibble::as_tibble() %>% 
#   mutate(batch_ids = batch_ids[sample],
#          sample_ids = sample_ids[sample])
# batch_ids <- batch_ids[sample]
# sample_ids <- sample_ids[sample]
if(FALSE){
combined_expr <- create_sample(combined_expr,
                               batch_ids,
                               sample_ids)
save(combined_expr, all_markers, file = "data/01_preprocessed_data.Rdata")
}
