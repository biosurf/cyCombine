#### Compile FCS files ----


#' Compile all .fcs files in a directory to a flowset
#'
#'
#'
#' @importFrom readxl read_xlsx
#' @importFrom readr read_csv
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom stringr str_remove
#' @importFrom flowCore read.flowSet fsApply
#' @param data_dir Directory containing the .fcs files
#' @param metadata Can be either a filename or data.frame of the metadata file If NULL, sample and batch ids are not predicted nor returned
#' @param sample_col The column in the metadata filename containing the sample ids. If NULL, sample ids will be the file names
#' @param batch_col The column in the metadata filename containing the sample ids
#' @param pattern The pattern to use to find the files in the folder
#'
#' @examples
#' fcs <- compile_fcs(data_dir = "_data/raw", metadata = "meta/metadata.csv", sample_col = "sample", batch_col = "batch", filename_col = "fcs_filenames")
#' @export
compile_fcs <- function(data_dir,
                        metadata = NULL,
                        sample_col = NULL,
                        batch_col = "Batch",
                        filename_col = "FCS_name",
                        pattern = "\\.fcs"){
  # Specifying files to use
  files <- list.files(data_dir,
                      pattern = pattern,
                      recursive = FALSE,
                      full.names = TRUE)
  message(paste("Read", length(files), "file names to process"))

    # Read all the data files
  fcs_raw <- files %>%
    flowCore::read.flowSet(transformation = FALSE,
                           truncate_max_range = FALSE,
                           emptyValue = FALSE)
  if(is.null(metadata)) return(fcs_raw)

  # Get metadata
  if(class(metadata) == "character"){
    if(endsWith(metadata, suffix = ".xlsx")){
      metadata <- file.path(data_dir, metadata) %>%
        readxl::read_xlsx()
    } else if(endsWith(meta_filename, suffix = ".csv")){
      metadata <- file.path(data_dir, metadata) %>%
        readr::read_csv()
    } else {
      stop(stringr::str_c("Sorry, file", metadata, "is not in a supported format. Please use a .xlsx or .csv file.",
                          sep = " "))
    }
  }


  if(!endsWith(metadata[[filename_col]][1], ".fcs")){
    metadata[[filename_col]] <- paste0(metadata[[filename_col]], ".fcs")
  }

  # Get sample and batch ids
  if (is.null(sample_col)){
    sample_ids <- basename(files) %>%
      stringr::str_remove(".fcs") %>%
      rep(flowCore::fsApply(fcs_raw, nrow))
    batch_ids <- metadata[[batch_col]][match(sample_ids, stringr::str_remove(metadata[[filename_col]], ".fcs"))] %>%
      as.factor() %>%
      rep(flowCore::fsApply(fcs_raw, nrow))
  } else{
    sample_ids <- metadata[[sample_col]][match(basename(files), metadata[[filename_col]])] %>%
      rep(flowCore::fsApply(fcs_raw, nrow)) %>%
      stringr::str_remove(".fcs")
    batch_ids <- metadata[[batch_col]][match(basename(files), metadata[[filename_col]])] %>%
      as.factor() %>%
      rep(flowCore::fsApply(fcs_raw, nrow))
  }

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
#'
#' @param flowset The flowset to convert
#' @param sample_ids Vector containing the sample ids
#' @param batch_ids Vector containing the batch ids
#' @param down_sample If TRUE, the output will be down-sampled to size sample_size
#' @param sample_size The size to down-sample to
#' @param seed The seed to use for down-sampling
#' @param panel Panel as a data.frame. Is used to predict colnames
#' @param panel_channel Only used if panel is given. It is the column name in the panel data.frame that contains the channel names
#' @param panel_antigen Only used if panel is given. It is the column name in the panel data.frame that contains the antigen names
#'
#' @examples
#' df <- convert_flowset(flowset = fcs$fcs_raw, sample_ids = fcs$sample_ids, batch_ids = fcs$batch_ids, down_sample = FALSE)
#' @export
convert_flowset <- function(flowset,
                            sample_ids = NULL,
                            batch_ids = NULL,
                            down_sample = TRUE,
                            sample_size = 100000,
                            seed = 473,
                            panel = NULL,
                            panel_channel = "fcs_colname",
                            panel_antigen = "antigen"){
  # Down sampling setup
  if(down_sample){
    # To down sample within fsApply
    nrows <- flowCore::fsApply(flowset, nrow)
    tot_nrows <- sum(nrows)
    message(paste("Down sampling to", sample_size, "samples"))
    set.seed(seed)
    sample <- sample(1:tot_nrows, sample_size) %>%
      # Sorting here enables major resource savings when down-sampling
      sort()
    if(!is.null(sample_ids) & !is.null(batch_ids)){
      sample_ids <- sample_ids[sample]
      batch_ids <- batch_ids[sample]
    }


  }

  message("Extracting expression data")
  fcs_data <- flowset %>%
    purrr::when(down_sample ~ flowCore::fsApply(., fcs_sample,
                                                sample = sample,
                                                nrows = nrows),
                ~ flowCore::fsApply(., Biobase::exprs)) %>%
    tibble::as_tibble()

  # Clean column names
  if (!is.null(panel)){
    cols <- match(colnames(fcs_data), panel$fcs_colname) %>%
      .[!is.na(.)]
    col_names <- panel$antigen[cols] %>%
      stringr::str_remove_all("[ -]") %>%
      stringr::str_remove_all("\\d+[A-Za-z]+_")

    fcs_data <- fcs_data %>%
      dplyr::select(dplyr::all_of(panel$fcs_colname))
  }else{
    col_names <- flowset[[1]] %>%
      flowCore::parameters() %>%
      Biobase::pData() %>%
      dplyr::pull(desc) %>%
      stringr::str_remove_all("[ -]") %>%
      stringr::str_remove_all("\\d+[A-Za-z]+_")
  }

  if(!is.null(sample_ids) & !is.null(batch_ids)){
    fcs_data <- fcs_data %>%
      dplyr::mutate(batch = batch_ids,
                    sample = sample_ids,
                    id = 1:nrow(.))
    colnames(fcs_data) <- c(col_names, "batch", "sample", "id")
  } else{
    fcs_data <- fcs_data %>%
      dplyr::mutate(id = 1:nrow(.))
    colnames(fcs_data) <- c(col_names, "id")
  }


  message("Your flowset is now converted into a dataframe.")
  return(fcs_data)
}




#' Extract from a flowset given a sample of indices
#'
#' @importFrom purrr accumulate
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
#' @param df The dataframe to transform
#' @param markers The markers to transform on
#' @param cofactor The cofactor to use when transforming
#' @family preprocess
#' @examples
#' preprocessed <- df %>%
#'   transform_asinh(markers = markers)
#' @export
transform_asinh <- function(df, markers, cofactor = 5){
  message(paste("Transforming data using asinh with a cofactor of", cofactor))
  transformed <- df %>%
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
#' This is a wrapper function that takes you from a directory of .fcs files to a transformed dataframe.
#'
#'
#' @inheritParams compile_fcs
#' @inheritParams convert_flowset
#' @inheritParams transform_asinh
#' @export
preprocess <- function(data_dir,
                       markers,
                       meta_filename = NULL,
                       sample_col = NULL,
                       batch_col = "Batch",
                       filename_col = "FCS_name",
                       pattern = "\\.fcs",
                       down_sample = TRUE,
                       sample_size = 500000,
                       seed = 473,
                       cofactor = 5){
  if(class(data_dir) != "character"){
    stop(paste("This function only works with a directory of .fcs files.",
               "If you have already loaded a flowset, please see the convert_flowset() function.",
               sep = "\n"))
  }
  # Compile directory to flowset
  fcs <- data_dir %>%
    compile_fcs(meta_filename = meta_filename,
                sample_col = sample_col,
                batch_col = batch_col,
                filename_col = filename_col,
                pattern = pattern)
  # Convert flowset to dataframe
  fcs_data <- fcs$fcs_raw %>%
    convert_flowset(batch_ids = fcs$batch_ids,
                    sample_ids = fcs$sample_ids,
                    down_sample = down_sample,
                    sample_size = sample_size,
                    seed = seed) %>%
    # Transform dataset with asinh
    transform_asinh(markers = markers,
                    cofactor = cofactor)
  message("Done!")
  return(fcs_data)
}



