#' Exporting a dataframe to SingleCellObject
#'
#' Conversion of dataframe into separate data and metadata objects for subsequent transformation.
#'  Rename variable names to fit the requirements of SCE-based tools
#'
#' @param df Tibble with expression values and metadata
#' @param markers Markers to include in exprs object of SCE. If NULL, markers will be found using the \code{\link{get_markers}} function.
#' @param non_markers Non-markers to include as colData in SCE. If NULL, non_markers will be based on cyCombine::non_markers.
#' @param panel_channel It is the column name in the df that contains the sample names. Defaults to 'sample'.
#' @param panel Optional: Panel as a data.frame. Should have colnames Channel, Marker, Type unless otherwise specified in the panel_ args.
#' @param panel_channel Optional: Only used if panel is given. It is the column name in the panel data that contains the channel names
#' @param panel_antigen Optional: Only used if panel is given. It is the column name in the panel data that contains the antigen names
#' @param panel_type Optional: Only used if panel is given. It is the column name in the panel data that contains the antigen types (none, state, type).
#'  "none" will be excluded from SCE.
#' @param sample_col The name of the column containing sample IDs
#' @family export
#' @examples
#' \dontrun{
#' sce <- df %>%
#'   df2SCE(markers = markers, non_markers = NULL, panel = panel)
#'   }
#' @export
df2SCE <- function(df, markers = NULL, non_markers = NULL, sample_col = 'sample', panel = NULL,
                   panel_channel = 'Channel', panel_antigen = 'Marker', panel_type = 'Type') {

  # Check for packages
  cyCombine:::missing_package("SingleCellExperiment", "Bioc")

  message("Converting dataframe to SingleCellExperiment object...")

  # Get markers and check
  if (is.null(markers)){
    # Get markers
    markers <- df %>%
      cyCombine::get_markers()
  }

  sapply(markers, function(x) {cyCombine:::check_colname(colnames(df), x, "df")})

  # Extract expression data and transpose to fit SCE format
  exprs <- t(dplyr::select(df, dplyr::all_of(markers)))


  # Get the non markers and prepare column data
  if (is.null(non_markers)){
    # Get non markers
    if (sample_col == 'sample') {
      non_markers <- cyCombine::non_markers
    } else {
      non_markers <- c(cyCombine::non_markers, sample_col)
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
    stop("Error, non of the non_markers/sample_col are available in the dataframe. You cannot make an SCE without sample names.")
  }


  # Prepare the experiment info table - first identify true meta data columns
  colCounting <- colData %>%
    dplyr::group_split(sample_id) %>%
    lapply(function(x) {(apply(x, 2, function(y) {length(table(y)) == 1}))})

  col_stability <- lapply(colCounting, function(z) {names(which(z))}) %>% unlist() %>% table()
  stable_cols <- names(which(col_stability == length(colCounting)))

  # Swap ordering
  stable_cols <- c('sample_id', stable_cols[stable_cols != 'sample_id'])

  experiment_info <- colData %>%
    dplyr::select(dplyr::all_of(stable_cols)) %>%
    dplyr::group_by(sample_id) %>%
    dplyr::mutate(n_cells = dplyr::n()) %>%
    dplyr::distinct(sample_id, .keep_all = TRUE) %>%
    as.data.frame()



  # Prepare row data if available
  if (!is.null(panel) & is.data.frame(panel)) {
    # Check presence of panel's data names
    sapply(c(panel_channel, panel_antigen, panel_type), function(x) {cyCombine:::check_colname(colnames(panel), x, "panel")})

    # Exclude none's
    rowData <- panel %>%
      dplyr::filter(.data[[panel_type]] != "none")


    # Change colnames to fit SCE standard
    if (panel_channel != 'channel_name') {
      rowData <- rowData %>%
        dplyr::rename(channel_name = panel_channel)
    }
    if (panel_antigen != 'marker_name') {
      rowData <- rowData %>%
        dplyr::rename(marker_name = panel_antigen)
    }
    if (panel_type != 'marker_class') {
      rowData <- rowData %>%
        dplyr::rename(marker_class = panel_type)
    }

    # Change marker names to exclude spaces and dashes
    rowData[['marker_name']] <- rowData[['marker_name']] %>%
      stringr::str_remove_all("^\\d+[A-Za-z]+_") %>%
      stringr::str_remove_all("[ _-]")


  } else if (!is.null(rowData) & !is.data.frame(rowData)) {
    stop("The provided panel is not a data frame.")

  } else {
    rowData = NULL

  }

  # Creating the SCE
  sce <- SingleCellExperiment::SingleCellExperiment(list(exprs = exprs),
                                                    colData = colData,
                                                    rowData = rowData,
                                                    metadata = list('experiment_info' = experiment_info))

  return(sce)
}
