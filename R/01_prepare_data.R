#### Compile FCS files ----


#' Compile all .fcs files in a directory to a flowset
#'
#' A simple function to compile a directory of FCS files into a flowSet object using flowCore::read.flowSet().
#'  Use the pattern argument to select only a subset of FCS files.
#'
#' @importFrom readxl read_xlsx
#' @importFrom readr read_csv
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom flowCore read.flowSet fsApply
#' @param data_dir Directory containing the .fcs files
#' @param pattern The pattern to use to find the files in the folder
#' @inheritParams flowCore::read.FCS
#' @family dataprep
#' @examples
#' fcs <- compile_fcs(data_dir = "_data/raw", pattern = "\\.fcs")
#' @export
compile_fcs <- function(data_dir,
                        pattern = "\\.fcs",
                        column.pattern = NULL,
                        invert.pattern = FALSE){

  # Error checking
  if(data_dir %>% endsWith("/")) data_dir <- data_dir %>% stringr::str_sub(end = -2)


  # Specifying files to use
  files <- list.files(data_dir,
                      pattern = pattern,
                      recursive = FALSE,
                      full.names = TRUE) %>%
    sort()
  if(length(files) == 0) stop("No files found in folder \"", data_dir, "\"")

  # Read the data files
  message(paste("Reading", length(files), "files to a flowSet.."))
  fcs_raw <- files %>%
    flowCore::read.flowSet(transformation = FALSE,
                           truncate_max_range = FALSE,
                           emptyValue = FALSE,
                           column.pattern = column.pattern,
                           invert.pattern = invert.pattern)

  return(fcs_raw)
}






#' Convert a flowSet into a tibble
#'
#' Use this function to convert a flowSet into the tibble object that the remaining
#'  functions in cyCombine relies on.
#'  A tibble is a Tidyverse implementation of a data.frame and can be treated as a such.
#'  The majority of arguments revolves adding relevant info from the metadata file/object.
#'  The panel argument is included to adjust the output column names using a panel with channel and antigen columns.
#'  Bear in mind the column names will be altered with the following:
#'  \code{stringr::str_remove_all("^\\d+\[A-Za-z\]+_") %>% stringr::str_remove_all("\[ _-\]")}.
#'
#'
#' @importFrom tibble as_tibble
#' @importFrom stringr str_remove
#' @importFrom purrr when
#' @importFrom stringr str_remove_all
#' @importFrom Biobase exprs pData
#' @importFrom flowCore parameters
#'
#' @param flowset The flowset to convert
#' @param metadata Optional: Can be either a filename or data.frame of the metadata file. Please give the full path from working directory to metadata file
#' @param sample_ids Optional: If a character, it should be the sample column in the metadata. If its a vector, it should have the same length as the total flowset. If NULL, sample ids will be the file names. If a single value, all rows will be assigned this value.
#' @param batch_ids Optional: If a character, it should be the column in the metadata containing the batch ids. If its a vector, it should have the same length as the total flowset. If a single value, all rows will be assigned this value.
#' @param filename_col Optional: The column in the metadata containing the fcs filenames. Needed if metadata is given, but sample_ids is not
#' @param condition Optional: The column in the metadata containing the condition. Will be used as the covariate in ComBat, but can be specified later. You may use this to add a different column of choice, in case you want to use a custom column in the ComBat model matrix.
#' @param anchor Optional: The column in the metadata referencing the anchor samples (control references). Will be used as a covariate in ComBat, if specified. Please be aware that this column may be confounded with the condition column. You may use this to add a different column of choice, in case you want to use a custom column in the ComBat model matrix. You may use a custom column name, but it is good practice to add the name to the 'non_markers' object exported by cyCombine, to reduce the risk of unexpected errors.
#' @param down_sample If TRUE, the output will be down-sampled to size sample_size
#' @param sample_size The size to down-sample to
#' @param sampling_type The type of down-sampling to use. "random" to randomly select cells across the entire dataset, "batch_ids" to sample evenly (sample_size) from each batch, or "sample_ids" sample evenly (sample_size) from each sample.
#' @param seed The seed to use for down-sampling
#' @param panel Optional: Panel as a filename or data.frame. Is used to define colnames from the panel_antigen column
#' @param panel_channel Optional: Only used if panel is given. It is the column name in the panel data.frame that contains the channel names
#' @param panel_antigen Optional: Only used if panel is given. It is the column name in the panel data.frame that contains the antigen names
#' @family dataprep
#' @examples
#' df <- convert_flowset(flowset = flowset,
#'  metadata = file.path(data_dir, "metadata.csv"),
#'  filename_col = "FCS_files",
#'  sample_ids = "sample_id",
#'  batch_ids = "batch_ids",
#'  down_sample = FALSE)
#' @export
convert_flowset <- function(flowset,
                            metadata = NULL,
                            filename_col = "filename",
                            sample_ids = NULL,
                            batch_ids = NULL,
                            condition = NULL,
                            anchor = NULL,
                            down_sample = TRUE,
                            sample_size = 500000,
                            sampling_type = "random",
                            seed = 473,
                            panel = NULL,
                            panel_channel = "fcs_colname",
                            panel_antigen = "antigen"){
  # Extract information necessary both with and without included metadata
  ## FlowSet row numbers
  nrows <- flowCore::fsApply(flowset, nrow)
  ## To down sample within fsApply
  tot_nrows <- sum(nrows)
  ids <- 1:tot_nrows
  ## File names from flowset
  files <- flowset@phenoData %>%
    rownames() %>%
    basename()
  # Add metadata information (long)
  if(!is.null(metadata)){
    # Error handling
    if(is.null(filename_col)){
      stop("Please specify a filename_col.")
    }
    if("character" %in% class(metadata)){
      if(!file.exists(file.path(metadata))){
        stop("File \"", file.path(metadata), "\" was not found")
      }

      # Get metadata
      if(endsWith(metadata, suffix = ".xlsx")){
        metadata <- suppressMessages(file.path(metadata) %>%
                                       readxl::read_xlsx())
      } else if(endsWith(metadata, suffix = ".csv")){
        metadata <- suppressMessages(file.path(metadata) %>%
                                       readr::read_csv())
      } else {
        stop(stringr::str_c("Sorry, file", metadata, "is not in a supported format. Please use a .xlsx or .csv file.\n",
                            "Alternatively, a data.frame of the metadata can be used.",
                            sep = " "))
      }
    }

    # Check for errors in metadata columns
    md_cols <- colnames(metadata)
    check_colname(md_cols, filename_col)


    # Extract info from metadata
    if(!endsWith(metadata[[filename_col]][1], ".fcs")){
      metadata[[filename_col]] <- paste0(metadata[[filename_col]], ".fcs")
    }
    # Remove files from metadata
    metadata <- metadata[match(files, metadata[[filename_col]]),]
    # Check that all files are represented in metadata
    if(any(files %!in% metadata[[filename_col]])){
      missing_files <- files[files %!in% metadata[[filename_col]]]
      warning(stringr::str_c("The samples ", missing_files, " were not found in the metadata file and the sample_ids will be set to NA."))
      }

    # Get sample ids
    if (is.null(sample_ids)){
      sample_ids <- metadata[[filename_col]] %>%
        stringr::str_remove(".fcs") %>%
        rep(nrows)
    } else if (length(sample_ids) == 1){
      check_colname(md_cols, sample_ids)
      sample_ids <- metadata[[sample_ids]][match(files, metadata[[filename_col]])] %>%
        stringr::str_remove(".fcs") %>%
        rep(nrows)
    }

    # Get batch ids
    if (!is.null(batch_ids)){
      if(length(batch_ids) == 1){
        check_colname(md_cols, batch_ids)
        batch_ids <- metadata[[batch_ids]][match(files, metadata[[filename_col]])] %>%
          as.factor() %>%
          rep(nrows)
      }
    }
    # Get condition
    if(!is.null(condition)){
      if(length(condition == 1)){
        check_colname(md_cols, condition)
        condition <- metadata[[condition]][match(basename(files), metadata[[filename_col]])] %>%
          as.factor() %>%
          rep(nrows)
      }
    }
    # Get anchor
    if(!is.null(anchor)){
      if(length(anchor == 1)){
        check_colname(md_cols, anchor)
        anchor <- metadata[[anchor]][match(basename(files), metadata[[filename_col]])] %>%
          as.factor() %>%
          rep(nrows)
      }
    }


  } else{ # If no metadata given
    if(is.null(sample_ids)){# Get sample ids from filenames
      sample_ids <- files %>%
        stringr::str_remove(".fcs") %>%
        rep(nrows)
    }
  }

  # Down sampling setup
  if(down_sample){
    message(paste("Down sampling to", sample_size, "samples"))
    set.seed(seed)
    if(sampling_type == "random"){
      sample <- sample(ids, sample_size) %>%
        # Sorting here enables major resource savings when down-sampling
        sort()
    } else if(sampling_type == "batch_ids" & !is.null(batch_ids)){ # even down-sampling from batches
      sample <- tibble::tibble(batch_ids, ids) %>%
        dplyr::group_by(batch_ids) %>%
        dplyr::slice_sample(n = sample_size) %>%
        dplyr::pull(ids) %>%
        sort()
    } else{ # Even down-sampling from samples
      sample <- tibble::tibble(sample_ids, ids) %>%
        dplyr::group_by(sample_ids) %>%
        dplyr::slice_sample(n = sample_size) %>%
        dplyr::pull(ids) %>%
        sort()
    }

    # Down-sample metadata columns
    ids <- ids[sample]
    if(!is.null(sample_ids) & length(sample_ids) > 1) sample_ids <- sample_ids[sample]
    if(!is.null(batch_ids) & length(batch_ids) > 1) batch_ids <- batch_ids[sample]
    if(!is.null(condition) & length(condition) > 1) condition <- condition[sample]
    if(!is.null(anchor) & length(anchor) > 1) anchor <- anchor[sample]
  }

  message("Extracting expression data..")
  fcs_data <- flowset %>%
    purrr::when(down_sample ~ flowCore::fsApply(., fcs_sample,
                                                sample = sample,
                                                nrows = nrows),
                ~ flowCore::fsApply(., Biobase::exprs)) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(id = ids) %>%
    dplyr::select(id, everything())

  # Clean column names
  if (!is.null(panel)){
    if("character" %in% class(panel)){
      if(endsWith(panel, suffix = ".xlsx")){
        panel <- suppressMessages(file.path(panel) %>%
                                       readxl::read_xlsx())
      } else if(endsWith(panel, suffix = ".csv")){
        panel <- suppressMessages(file.path(panel) %>%
                                       readr::read_csv())
      } else {
        stop(paste("Sorry, file", panel, "is not in a supported format. Please use a .xlsx or .csv file."))
      }
    }
    cols <- match(colnames(fcs_data), panel[[panel_channel]]) %>%
      .[!is.na(.)]
    col_names <- panel[[panel_antigen]][cols] %>%
      stringr::str_remove_all("^\\d+[A-Za-z]+_") %>%
      stringr::str_remove_all("[ _-]")

    fcs_data <- fcs_data %>%
      dplyr::select(id, dplyr::all_of(panel[[panel_channel]][cols]))
  }else{
    col_names <- flowset[[1]] %>%
      flowCore::parameters() %>%
      Biobase::pData() %>%
      dplyr::pull(desc) %>%
      stringr::str_remove_all("^\\d+[A-Za-z]+_") %>%
      stringr::str_remove_all("[ _-]")
  }
  colnames(fcs_data) <- c("id", col_names)

  # Add optional columns
  if(!is.null(sample_ids)) fcs_data$sample <- sample_ids
  if(!is.null(batch_ids)) fcs_data$batch <- batch_ids
  if(!is.null(condition)) fcs_data$condition <- condition
  if(!is.null(anchor)) fcs_data$anchor <- anchor

  message("Your flowset is now converted into a dataframe.")
  return(fcs_data)
}




#' Extract from a flowset given a sample of indices
#' @noRd
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
#'
#' Inverse sine transformation of a tibble.
#'  This function can also de-randomize data.
#'
#' @param df The tibble to transform
#' @param markers The markers to transform on
#' @param cofactor The cofactor to use when transforming
#' @param derand Derandomize. Should be TRUE for CyTOF data, otherwise FALSE.
#' @param .keep Keep all channels. If FALSE, channels that are not transformed are removed
#' @importFrom knitr combine_words
#' @family dataprep
#' @examples
#' uncorrected <- df %>%
#'   transform_asinh(markers = markers)
#' @export
transform_asinh <- function(df,
                            markers = NULL,
                            cofactor = 5,
                            derand = TRUE,
                            .keep = FALSE){
  if(is.null(markers)){
    markers <- df %>%
      cyCombine::get_markers()
  }
  if(any(markers %!in% colnames(df))){
    mes <- stringr::str_c("Not all given markers are in the data.\nCheck if the markers contain a _ or -:",
                 knitr::combine_words(markers),
                 "Columns:",
                 knitr::combine_words(colnames(df)),
                 sep = "\n"
    )
    stop(mes)
  }
  message(paste0("Transforming data using asinh with a cofactor of ", cofactor, ".."))
  transformed <- df %>%
    purrr::when(.keep ~ .,
                ~ dplyr::select_if(., colnames(.) %in% c(markers, cyCombine::non_markers))) %>%
    # Transform all data on those markers
    dplyr::mutate(dplyr::across(dplyr::all_of(markers),
                     .fns = function(x){
                       if(derand) asinh(ceiling(x)/cofactor) else asinh(x/cofactor)
                     }))
  return(transformed)
}

#### Wrapper function ----



#' Prepare a directory of .fcs files
#'
#' This is a wrapper function that takes you from a directory of .fcs files or a flowset to a transformed tibble.
#'
#'
#' @param flowset Optional: Prepare a flowset instead of a directory of fcs files
#' @inheritParams compile_fcs
#' @inheritParams convert_flowset
#' @inheritParams transform_asinh
#' @param transform If TRUE, the data will be transformed; if FALSE, it will not.
#' @family dataprep
#' @examples
#' uncorrected <- data_dir %>%
#'   prepare_data(metadata = "metadata.csv",
#'   markers = markers,
#'   filename_col = "FCS_name",
#'   batch_ids = "Batch",
#'   condition = "condition",
#'   down_sample = FALSE)
#' @export
prepare_data <- function(data_dir = NULL,
                         flowset = NULL,
                         markers = NULL,
                         pattern = "\\.fcs",
                         metadata = NULL,
                         filename_col = "filename",
                         sample_ids = NULL,
                         batch_ids = NULL,
                         condition = NULL,
                         anchor = NULL,
                         down_sample = TRUE,
                         sample_size = 500000,
                         sampling_type = "random",
                         seed = 473,
                         panel = NULL,
                         panel_channel = "fcs_colname",
                         panel_antigen = "antigen",
                         transform = TRUE,
                         cofactor = 5,
                         derand = TRUE,
                         .keep = FALSE){

  # Stop if no data is given
  if(is.null(data_dir) & is.null(flowset)) stop("No data given.")

  if(!is.null(data_dir)){
    # Remove slash at end of data_dir
    if(data_dir %>% endsWith("/")) data_dir <- data_dir %>% stringr::str_sub(end = -2)

    # Compile directory to flowset
    flowset <- data_dir %>%
      cyCombine::compile_fcs(pattern = pattern)

    # Look for metadata in data_dir
    if(!is.null(metadata)){
      if(!"data.frame" %in% class(metadata)){
        if(!file.exists(file.path(metadata)) & file.exists(file.path(data_dir, metadata))) metadata <- file.path(data_dir, metadata)
      }
    }
  }
  # Convert flowset to dataframe
  fcs_data <- flowset %>%
    cyCombine::convert_flowset(metadata = metadata,
                               filename_col = filename_col,
                               sample_ids = sample_ids,
                               batch_ids = batch_ids,
                               condition = condition,
                               anchor = anchor,
                               down_sample = down_sample,
                               sample_size = sample_size,
                               sampling_type = sampling_type,
                               seed = seed,
                               panel = panel,
                               panel_channel = panel_channel,
                               panel_antigen = panel_antigen) %>%
    # Transform dataset with asinh
    purrr::when(transform ~ cyCombine::transform_asinh(., markers = markers,
                                                       cofactor = cofactor,
                                                       derand = derand,
                                                       .keep = .keep),
                ~ .)

  message("Done!")
  return(fcs_data)
}





