#### Compile FCS files ----


#' Compile all .fcs files in a directory to a flowset
#'
#' @importFrom readxl read_xlsx
#' @importFrom readr read_csv
#' @import dplyr
#' @import magrittr
#' @importFrom stringr str_remove
#' @importFrom flowCore read.flowSet fsApply
#' @param data_dir Directory containing the .fcs files
#' @export
compile_fcs <- function(data_dir, meta_filename){
  # Specifying files to use
  files <- list.files(data_dir,
                      pattern="\\.fcs",
                      recursive = FALSE,
                      full.names = TRUE)
  cat("Read", length(files), "file names to process", "\n",
              sep = " ")

  # Get metadata
  if(endsWith(meta_filename, suffix = ".xlsx")){
    meta_data <- file.path(data_dir, meta_filename) %>%
      readxl::read_xlsx()
  } else if(endsWith(meta_filename, suffix = ".csv")){
    meta_data <- file.path(data_dir, meta_filename) %>%
      readr::read_csv()
  } else {
    stop(stringr::str_c("Sorry, file", meta_filename, "is not in a supported format. Please use a .xlsx or .csv file.",
                        sep = " "))
  }


  # Read all the data files
  fcs_raw <- files %>%
    flowCore::read.flowSet(transformation = FALSE,
                           truncate_max_range = FALSE,
                           emptyValue = FALSE)


  # Get sample names
  sample_ids <- basename(files) %>%
    stringr::str_remove(".fcs") %>%
    rep(flowCore::fsApply(fcs_raw, nrow))
  batch_ids <- meta_data$Batch[match(sample_ids, meta_data$FCS_name)] %>%
    as.factor()
  return(list("fcs_raw" = fcs_raw,
              "sample_ids" = sample_ids,
              "batch_ids" = batch_ids))
}

# Load data - the loaded contains objects "fcs_raw", "sample_ids",
# "batch_ids", "all_markers", which are the raw flowSet, sample ids per row, batch
# ids per row and all measured markers in the panel

#' Convert a flowset into a workable dataframe
#'
#' @importFrom tibble as_tibble
#' @importFrom purrr when
#' @importFrom stringr str_remove_all
#' @importFrom Biobase exprs pData
#' @importFrom flowCore parameters
#' @export
convert_flowset <- function(flowset,
                            sample_ids,
                            batch_ids,
                            down_sample = TRUE,
                            sample_size = 100000,
                            seed = 473){
  # Down sampling setup
  if(down_sample){
    cat("Down sampling to", sample_size, "samples\n",
                         sep = " ")
    set.seed(seed)
    sample <- sample(1:length(sample_ids), sample_size) %>%
      # Sorting here enables major resource savings when down-sampling
      sort()
    sample_ids <- sample_ids[sample]
    batch_ids <- batch_ids[sample]
    # To down sample within fsApply
    nrows <- flowCore::fsApply(flowset, nrow)
  }

  cat("Extracting expression data and adding sample and batch labels", "\n")
  fcs_data <- flowset %>%
    purrr::when(down_sample ~ flowCore::fsApply(., fcs_sample,
                                                sample = sample,
                                                nrows = nrows),
                ~ flowCore::fsApply(., Biobase::exprs)) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(batch = batch_ids,
                  sample = sample_ids,
                  id = 1:nrow(.))
  # {if(down_sample) {set.seed(seed); sample_n(., sample_size)} else .}


  # Clean column names
  col_names <- flowset[[1]] %>%
    flowCore::parameters() %>%
    Biobase::pData() %>%
    dplyr::pull(desc) %>%
    stringr::str_remove_all("[ -]") %>%
    stringr::str_remove_all("\\d+[A-Za-z]+_")

  colnames(fcs_data) <- c(col_names, "batch", "sample", "id")
  cat("Your flowset is now converted into a dataframe.\n")
  return(fcs_data)
}

#' Extract from a flowset given a sample of indices
#'
#' @importFrom purrr accumulate
#' @export
fcs_sample <- function(flowframe, sample, nrows, seed = 473){
  nrows_acc <- c(0, nrows %>%
    purrr::accumulate(`+`))
  ff_number <- stringr::str_c("^", nrow(flowframe), "$") %>%
    grep(nrows)

  ff_sample <- sample - nrows_acc[ff_number]
  ff_sample <- ff_sample[ff_sample > 0]
  ff_sample <- ff_sample[ff_sample <= nrows[ff_number]]

  ff <- flowframe %>%
    Biobase::exprs()
  ff <- ff[ff_sample, ]

  return(ff)
}




#### Data transformation ----

#' Transform data using asinh
#'
#' @family preprocess
#' @export
transform_asinh <- function(input, markers, cofactor = 5){
  cat("Transforming data using asinh with a cofactor of", cofactor, "\n")
  transformed <- input %>%
    # Select markers of interest
    dplyr::select(dplyr::all_of(c(markers, "batch", "sample", "id"))) %>%
    # Transform all data on those markers
    dplyr::mutate_at(.vars = all_of(markers),
              .funs = function(x) asinh(ceiling(x)/cofactor)) %>%
    arrange(id)
  return(transformed)
}

#### Workflow function ----



#' Preprocess a directory of .fcs files
#'
#' @export
preprocess <- function(data_dir,
                       meta_filename,
                       markers,
                       down_sample = TRUE,
                       sample_size = 100000,
                       seed = 473,
                       cofactor = 5){
  if(class(data_dir) != "character"){
    stop(paste("This function only works with a directory of .fcs files.",
               "If you have already loaded a flowset, please see the convert_flowset() function.",
               sep = "\n"))
  }
  # Compile directory to flowset
  raw_flowset <- data_dir %>%
    compile_fcs(meta_filename = meta_filename)
  # Convert flowset to dataframe
  fcs_data <- raw_flowset$fcs_raw %>%
    convert_flowset(batch_ids = raw_flowset$batch_ids,
                    sample_ids = raw_flowset$sample_ids,
                    down_sample = down_sample,
                    sample_size = sample_size,
                    seed = seed) %>%
    # Transform dataset with asinh
    transform_asinh(markers = markers,
                    cofactor = cofactor)
  cat("Done!\n")
  return(fcs_data)
}



