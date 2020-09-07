#### Compile FCS files ----


#' Compile all .fcs files in a directory
#'
#' @importFrom Biobase exprs
#' @importFrom readxl read_xlsx
#' @import dplyr
#' @import stringr
#' @importFrom tibble as_tibble
#' @importFrom purrr when
#' @importFrom flowCore read.flowSet parameters fsApply
#' @param data_dir Directory containing the .fcs files
#' @export
compile_fcs <- function(data_dir, down_sample = FALSE, sample_size = 100000, seed = 473){
  # Specifying files to use
  files <- list.files(data_dir,
                      pattern="\\.fcs",
                      recursive = FALSE,
                      full.names = TRUE)
  print(str_c("Read", length(files), "file names to process",
              sep = " "))

  # Get metadata
  meta_data <- str_c(data_dir, "/CyTOF samples cohort.xlsx",
                     sep = "") %>%
    readxl::read_xlsx()

  # Read all the data files
  fcs_raw <- files %>%
    flowCore::read.flowSet(transformation = FALSE,
                           truncate_max_range = FALSE,
                           emptyValue = FALSE)

  # Clean column names
  colnames(fcs_raw) <- fcs_raw[[1]] %>%
    flowCore::parameters() %>%
    Biobase::pData() %>%
    pull(desc) %>%
    stringr::str_remove_all("[ -]") %>%
    stringr::str_remove_all("\\d+[A-Za-z]+_")


  # Get sample names
  sample_ids <- basename(files) %>%
    stringr::str_remove(".fcs") %>%
    rep(flowCore::fsApply(fcs_raw, nrow))

  print("Extracting expression data and adding sample and batch labels")
  fcs_data <- fcs_raw %>%
    flowCore::fsApply(Biobase::exprs) %>%
    as_tibble() %>%
    mutate(Batch = meta_data$Batch[match(sample_ids, meta_data$FCS_name)] %>%
             as.factor(),
           Sample = sample_ids) %>%
    # {if(down_sample) {set.seed(seed); sample_n(., sample_size)} else .}
    purrr::when(down_sample ~ slice_sample(., n = sample_size),
                ~ .)

  print("Done")
  return(fcs_data)
}

#### Prepare FlowSet ----
#' Prepare and downsample flowset
#'
#' @export
prepare_flowset <- function(flowset,
                    batch_ids,
                    sample_ids,
                    sample_size = 100000,
                    seed = 473){
  print("Extracting expression matrix")
  fcs_data <- flowset %>%
    flowCore::fsApply(Biobase::exprs) %>%
    create_sample(batch_ids,
                  sample_ids,
                  sample_size = sample_size,
                  seed = seed)
  return(fcs_data)
}

#' Create subsample of combined expression matrix
#'
#'
#' @importFrom tibble as_tibble
create_sample <- function(fcs_data, batch_ids, sample_ids,
                          sample_size = 100000,
                          seed = 473){
  print(paste("Down-sampling to", sample_size, "samples"))
  set.seed(seed)
  sample <- sample(1:nrow(fcs_data), sample_size)
  sample_set <- fcs_data[sample, ] %>%
    as_tibble() %>%
    mutate(Batch = batch_ids[sample],
           Sample = sample_ids[sample])
  return(sample_set)
}


#### Data pre-processing ----

# Load data - the loaded contains objects "fcs_raw", "sample_ids",
# "batch_ids", "all_markers", which are the raw flowSet, sample ids per row, batch
# ids per row and all measured markers in the panel


#' Transform data using asinh
#'
#'
#'
#'
#' @importFrom Biobase pData exprs
#' @import magrittr
#' @family preprocess
#' @export
transform_asinh <- function(input, markers, cofactor = 5){
  print(paste("Transforming data using asinh with a cofactor of", cofactor))
  transformed <- input %>%
    # rename_at(.vars = all_of(panel_fcs$name), .funs = function(x)gsub(' ', '', gsub('-', '', gsub("\\d+[A-Za-z]+_", "", x)))) %>%
    # rename_at(.vars = all_of(panel_fcs$name), .funs = list(str_replace(., "[ -]", ""),
    # str_replace(., "-", ""),
    # str_replace(., "\\d+[A-Za-z]+_", ""))) %>% View()
    select(all_of(c(markers, "Batch", "Sample"))) %>%
    mutate_at(.vars = all_of(markers),
              .funs = function(x) asinh(ceiling(x)/cofactor))
  return(transformed)
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



#' Preprocess FlowSet data
#'
#' @importFrom Biobase exprs
#' @export
preprocess <- function(compiled_fcs,
                       markers,
                       sample_size = 100000,
                       seed = 473){
  if("data.frame" %!in% class(compiled_fcs) || is.null(compiled_fcs$Batch)){
    stop(paste("Please make sure the input is in the right format",
               "Use the prepare() function to convert a fcs expression matrix, batch_ids, and sample_ids into a data.frame",
               "Or compile your fcs files using the compile_fcs() function.",
               sep = "\n"))
  }
  # print("Extracting objects")
  # Extract objects
  # fcs_raw <- input$fcs_raw
  # all_markers <- input$all_markers
  # sample_ids <- input$sample_ids
  # batch_ids <- input$batch_ids
  #
  # panel_fcs <- fcs_raw[[1]] %>%
  #   flowCore::parameters() %>%
  #   Biobase::pData()
  #
  # print("Extracting expression data")
  # combined_expr <- fcs_raw %>%
  #   flowCore::fsApply(Biobase::exprs) %>%            # Extract expression data
  print(paste("Down-sampling to", sample_size, "samples"))
  set.seed(seed)
  preprocessed_data <- compiled_fcs %>%
    # create_sample(sample_size = sample_size,
    #               seed = seed) %>%
    transform_asinh(markers = markers)
  print("Done")
  return(preprocessed_data)
}



