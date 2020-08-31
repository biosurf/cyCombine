#### Data pre-processing ####

# Load data - the loaded contains objects "fcs_raw", "sample_ids",
# "batch_ids", "all_markers", which are the raw flowSet, sample ids per row, batch
# ids per row and all measured markers in the panel


#' Transform data using asinh
#'
#'
#'
#'
#' @importFrom flowCore parameters fsApply
#' @importFrom Biobase pData exprs
#' @import magrittr
#' @family preprocess
transform_asinh <- function(fcs_raw){
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

#' Create subsample of combined expression matrix
#'
#'
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate
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

#'
#'
#' @importFrom flowCore fsApply
preprocess <- function(input,
                       sample_size = 100000,
                       seed = 473){
  # Extract objects
  fcs_raw <- input$fcs_raw
  all_markers <- input$all_markers
  sample_ids <- input$sample_ids
  batch_ids <- input$batch_ids

  #
  combined_expr <- fcs_raw %>%
    transform_asinh() %>%
    flowCore::fsApply(exprs) %>%            # Extract expression data
    create_sample(batch_ids = batch_ids,
                  sample_ids = sample_ids,
                  sample_size = sample_size,
                  seed = seed)

  return(list("data" = combined_expr,
              "markers" = all_markers))
}



