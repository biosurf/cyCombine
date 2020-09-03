#### Data pre-processing ####

# Load data - the loaded contains objects "fcs_raw", "sample_ids",
# "batch_ids", "all_markers", which are the raw flowSet, sample ids per row, batch
# ids per row and all measured markers in the panel


#' Transform data using asinh
#'
#'
#'
#'
#' @import flowCore
#' @importFrom Biobase pData exprs
#' @import magrittr
#' @importFrom dplyr select mutate_at
#' @family preprocess
#' @export
transform_asinh <- function(input, markers, cofactor = 5, panel_fcs){
  print(paste("Transforming data using asinh with a cofactor of", cofactor))
  colnames(input) <- c(gsub('[ -]', '', gsub("\\d+[A-Za-z]+_", "", panel_fcs$desc)), "Batch", "Sample")
  input <- input %>%
    # rename_at(.vars = all_of(panel_fcs$name), .funs = function(x)gsub(' ', '', gsub('-', '', gsub("\\d+[A-Za-z]+_", "", x)))) %>%
    # rename_at(.vars = all_of(panel_fcs$name), .funs = list(str_replace(., "[ -]", ""),
    # str_replace(., "-", ""),
    # str_replace(., "\\d+[A-Za-z]+_", ""))) %>% View()
    select(markers, Batch, Sample) %>%
    mutate_at(.vars = all_of(markers),
              .funs = function(x) asinh(ceiling(x)/cofactor))

  return(input)
}
# transform_asinh <- function(fcs_raw, markers){
#   # De-randomize and transform data using asinh
#   panel_fcs <- fcs_raw[[1]] %>%
#     flowCore::parameters() %>%
#     Biobase::pData()
#
#   panel_fcs$desc <- gsub(' ', '', gsub('-', '', gsub("\\d+[A-Za-z]+_", "", panel_fcs$desc)))
#
#   fcs <- flowCore::fsApply(fcs_raw, function(x, cofactor = 5){
#     colnames(x) <- panel_fcs$desc
#     expr <- Biobase::exprs(x)
#     expr <- asinh(ceiling(expr[, markers]) / cofactor)
#     exprs(x) <- expr
#     x
#   })
#   return(fcs)
# }

#' Create subsample of combined expression matrix
#'
#'
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate
#' @export
create_sample <- function(combined_expr,
                          batch_ids,
                          sample_ids,
                          sample_size = 100000,
                          seed = 473){
  print(paste("Down-sampling to", sample_size, "samples"))
  set.seed(seed)
  sample <- sample(1:nrow(combined_expr), 100000)
  combined_expr <- combined_expr[sample,] %>%
    tibble::as_tibble() %>%
    mutate(Batch = batch_ids[sample],
           Sample = sample_ids[sample])
  return(combined_expr)
}

#' Preprocess FlowSet data
#'
#' @importFrom flowCore fsApply
#' @export
preprocess <- function(input,
                       sample_size = 100000,
                       seed = 473){
  print("Extracting objects")
  # Extract objects
  fcs_raw <- input$fcs_raw
  all_markers <- input$all_markers
  sample_ids <- input$sample_ids
  batch_ids <- input$batch_ids

  panel_fcs <- fcs_raw[[1]] %>%
    flowCore::parameters() %>%
    Biobase::pData()

  print("Extracting expression data")
  combined_expr <- fcs_raw %>%
    flowCore::fsApply(exprs) %>%            # Extract expression data
    create_sample(batch_ids = batch_ids,
                  sample_ids = sample_ids,
                  sample_size = sample_size,
                  seed = seed) %>%
    transform_asinh(markers = all_markers,
                    panel_fcs = panel_fcs)
  print("Done")
  return(list("data" = combined_expr,
              "markers" = all_markers))
}



