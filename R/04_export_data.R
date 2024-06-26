#' Exporting a dataframe to SingleCellObject
#'
#' Conversion of dataframe into separate data and metadata objects for subsequent transformation.
#'  Rename variable names to fit the requirements of SCE-based tools
#'
#' @param df Tibble with expression values and metadata
#' @param markers Markers to include in exprs and counts object of SCE. If NULL, markers will be found using the \code{\link{get_markers}} function.
#' @param non_markers Non-markers to include as colData in SCE. If NULL, non_markers will be based on cyCombine::non_markers.
#' @param scatter The scatter columns to add to the SCE for FCS export.
#' @param sample_col It is the column name in the df that contains the sample names. Defaults to 'sample'.
#' @param panel Optional: Panel as a data.frame. Should have colnames Channel, Marker, Type unless otherwise specified in the panel_ args. Should be included if you want to store FCS files
#' @param panel_channel Optional: Only used if panel is given. It is the column name in the panel data that contains the channel names
#' @param panel_antigen Optional: Only used if panel is given. It is the column name in the panel data that contains the antigen names
#' @param panel_type Optional: Only used if panel is given. It is the column name in the panel data that contains the antigen types (none, state, type).
#'  "none" will be excluded from SCE. Set to NULL to disregard.
#' @param transform_cofactor The cofactor to use when reverse-transforming to raw counts
#' @param clean_names Cleans marker names.
#' @importFrom methods is
#' @family export
#' @examples
#' \dontrun{
#' sce <- df %>%
#'   df2SCE(markers = markers, non_markers = NULL, panel = panel)
#'   }
#' @export
df2SCE <- function(
    df,
    markers = NULL, non_markers = NULL,
    scatter = NULL,
    clean_names = FALSE,
    sample_col = "sample",
    panel = NULL, panel_channel = "Channel",
    panel_antigen = "Marker", panel_type = NULL,
    transform_cofactor = 5) {

  # Check for packages
  cyCombine:::missing_package("SingleCellExperiment", "Bioc")

  message("Converting dataframe to SingleCellExperiment object...")

  # Get the non markers and prepare column data
  if (is.null(non_markers)) {
    # Get non markers
    if (sample_col == "sample") {
      non_markers <- cyCombine::non_markers
    } else {
      non_markers <- c(cyCombine::non_markers, dplyr::all_of(sample_col))
    }
  }

  # Get the column data (from non markers)
  if (sample_col %in% colnames(df)) {
    colData <- df %>%
      dplyr::select(dplyr::any_of(non_markers))

  if (sample_col != 'sample_id') {
    colData <- colData %>%
      dplyr::rename(sample_id = sample_col)
  }

  } else {
    stop("Error, none of the non_markers/sample_col are available in the dataframe. You cannot make an SCE without sample names.")
  }

  # Get markers and check
  if (is.null(markers)) {
    # Get markers
    markers <- df %>%
      cyCombine::get_markers()
  }

  sapply(markers, function(x) {
    cyCombine:::check_colname(colnames(df), x, "df")})

  # Extract expression data and transpose to fit SCE format
  exprs <- t(dplyr::select(df, dplyr::all_of(c(scatter, markers))))


  # Prepare the experiment info table - first identify true meta data columns
  colCounting <- colData %>%
    dplyr::group_split(sample_id) %>%
    lapply(function(x) {(apply(x, 2, function(y) {length(table(y)) == 1}))})

  col_stability <- lapply(colCounting, function(z) {names(which(z))}) %>% unlist() %>% table()
  stable_cols <- names(which(col_stability == length(colCounting)))

  # Swap ordering
  stable_cols <- c("sample_id", stable_cols[stable_cols != "sample_id"])

  experiment_info <- colData %>%
    dplyr::select(dplyr::all_of(stable_cols)) %>%
    dplyr::group_by(sample_id) %>%
    dplyr::mutate(n_cells = dplyr::n()) %>%
    dplyr::distinct(sample_id, .keep_all = TRUE) %>%
    as.data.frame()



  # Prepare row data if available
  if (!is.null(panel) && is.data.frame(panel)) {
    # Check presence of panel's data names
    sapply(c(panel_channel, panel_antigen, panel_type), function(x) {
      cyCombine:::check_colname(colnames(panel), x, "panel")})

    rowData <- panel
    rm(panel)
    # Exclude none's
    if (!is(panel_type, "NULL")) {
      if ("none" %in% rowData[[panel_type]]){
        rowData <- rowData %>%
          dplyr::filter(.data[[panel_type]] != "none")
      }
    }


    # Change colnames to fit SCE standard
    if (panel_channel != "channel_name") {
      rowData <- rowData %>%
        dplyr::rename(channel_name = panel_channel)
    }
    if (panel_antigen != "marker_name") {
      rowData <- rowData %>%
        dplyr::rename(marker_name = panel_antigen)
    }
    if (panel_type != "marker_class" && !is(panel_type, "NULL")) {
      rowData <- rowData %>%
        dplyr::rename(marker_class = panel_type)
    }

    # Change marker names to exclude spaces and dashes
    if (clean_names) {
      rowData[["marker_name"]] <- rowData[["marker_name"]] %>%
        stringr::str_remove_all("^\\d+[A-Za-z]+_") %>%
        stringr::str_remove_all("[ _-]")
    }


    # Subset rowData to rows in exprs
    rowData <- rowData %>%
      dplyr::filter(.data[["marker_name"]] %in% rownames(exprs))

  } else {
    rowData <- NULL
    warning("To store as FCS files later, you should include panel information at this step.")
  }


  # Creating the SCE
  sce <- SingleCellExperiment::SingleCellExperiment(
    list(exprs = exprs,
         counts = exprs %>%
           t() %>%
           tibble::as_tibble() %>%
           transform_asinh(reverse = TRUE,
                           cofactor = transform_cofactor,
                           derand = FALSE,
                           .keep = TRUE,
                           markers = markers) %>%
           t()),
    colData = colData,
    rowData = rowData,
    metadata = list("experiment_info" = experiment_info)
    )
  message("Your SingleCellExperiment object is now created. The 'counts' assay contains reverse transformed expression values and 'exprs' contains expression values.")

  return(sce)
}

#' Convert SingelCellExperiment into flowSet and store as FCS files
#'
#' Wrapper function that makes it easier to go from a SCE to flowSet and written FCS files.
#'   The function uses CATALYST::sce2fcs and flowCore::write.flowSet to store the FCS files.
#'
#' @inheritParams CATALYST::sce2fcs
#' @inheritParams flowCore::write.flowSet
#' @param sce SingleCellExperiment to write to FCS files
#' @param outdir If given, the flowSet will be stored in FCS files
#' @param randomize (Default: FALSE) Logical determining whether counts are randomized for plotting purposes or not. Only works when assay = "counts" (default).
#' @family export
#' @examples
#' \dontrun{
#'  sce2FCS(sce, outdir = "fcs_files")
#'   }
#' @export
sce2FCS <- function(sce,
                    outdir = NULL,
                    split_by = "sample_id",
                    assay = "counts",
                    keep_dr = TRUE,
                    keep_cd = TRUE,
                    randomize = FALSE) {


  # Check CATALYST is installed
  cyCombine:::missing_package("CATALYST", repo = "Bioc")

  # If no panel information is available
  stopifnot("Your SingleCellExperiment should contain channel information to be stored.\n
            Consider rerunning df2sce with a panel to continue." = !is.null(CATALYST::channels(sce)))

  # Randomize counts
  if (randomize && assay == "counts") {
    SummarizedExperiment::assay(sce, "counts") <- SummarizedExperiment::assay(sce, "counts") %>%
      cyCombine:::randomize_matrix()
  } else if (randomize && assay != "counts") {
    warning("Please only use randomization with count values.")
  }

  # Convert to flowset
  fcs <- CATALYST::sce2fcs(sce,
                           keep_dr = keep_dr,
                           keep_cd = keep_cd,
                           split_by = split_by,
                           assay = assay)


  # Write FCS files
  if (!is.null(outdir)) {
    # Create output directory
    cyCombine:::check_make_dir(outdir)

    message("Writing FCS files to ", outdir)
    flowCore::write.flowSet(fcs, outdir = outdir)
  }

  return(fcs)
}
