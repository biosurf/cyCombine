# EMD ----

#' Compute EMD
#'
#' Given a dataframe, the Earth Movers Distance is computed, using the emdist package.
#'   A matrix is returned containing the EMD for each marker for each group of cells.
#'
#' @param df Dataframe to compute the EMD of
#' @param binSize The size of bins to use when binning data
#' @param cell_col Column name of df that contains cell population labels (or clusters)
#' @param batch_col Column name of df that contains batch numbers
#' @param markers Vector of the markers to calculate EMD for. If NULL, \code{\link{get_markers}} will be used to find markers
#' @family emd
#' @examples
#' emd <- compute_emd(df, markers = markers)
#' @export
compute_emd <- function(df,
                        binSize = 0.1,
                        cell_col = "label",
                        batch_col = "batch",
                        markers = NULL){

  # Check for package
  missing_package("emdist", "CRAN")
  missing_package("graphics", "CRAN")

  # Check colnames
  check_colname(colnames(df), cell_col, "df")
  check_colname(colnames(df), batch_col, "df")

  # Define cell columns as characters to avoid problems with factors
  df[[cell_col]] <- df[[cell_col]] %>%
    as.character()

  # Get markers if not given
  if(is.null(markers)){
    markers <- df %>%
      cyCombine::get_markers()
  }

  # Extract batches
  batches <- df %>%
    dplyr::pull(batch_col) %>%
    unique() %>%
    sort()
  # Extract cell types
  cellTypes <- df %>%
    dplyr::pull(cell_col) %>%
    unique() %>%
    sort()

  # Define limits for binning
  lower <- floor(min(df[,markers])) - binSize
  upper <- ceiling(max(df[,markers])) + binSize
  binLims <- c(lower, upper)

  # Create list of distribution matrices
  distr <- list()
  for (b in batches) {
    distr[[b]] <- list()
    for (cellType in cellTypes) {
      # Filter data on batch and cell type
      distr[[b]][[cellType]] <- df %>%
        dplyr::filter(.data[[cell_col]] == cellType,
                      .data[[batch_col]] == b) %>%
        dplyr::select(dplyr::all_of(markers)) %>%
        apply(2, function(x) {
          # Bin data
          bins <- seq(binLims[1], binLims[2], by = binSize)
          # bins <- seq(-10, 30, by = binSize)
          if (length(x) == 0) {
            rep(0, times = length(bins) - 1)
          }else{
            graphics::hist(x, breaks = bins,
                           plot = FALSE)$density
          }
        })
    }
  }
  # Compute emd from binned distributions
  distances <- list()
  for (cellType in cellTypes) {
    distances[[cellType]] <- list()
    for (marker in markers) {
      # Initiate distance matrix of a given cell type and marker
      distances[[cellType]][[marker]] <- matrix(NA, nrow = length(batches),
                                                ncol = length(batches), dimnames = list(batches,
                                                                                        batches))
      for (i in seq_along(batches)[-length(batches)]) {
        batch1 <- batches[i]
        for (j in seq(i + 1, length(batches))) {
          batch2 <- batches[j]
          # Get matrix for the two batches to be compared
          A <- matrix(distr[[batch1]][[cellType]][, marker])
          B <- matrix(distr[[batch2]][[cellType]][, marker])
          distances[[cellType]][[marker]][batch1, batch2] <- emdist::emd2d(A, B)
        }
      }
    }
  }

  return(distances)
}

#' Evaluate Earth Mover's Distance
#'
#' The function computes the Earth Mover's Distance of the two given datasets.
#'   Then the reduction is calculated as the relative change in total EMD.
#'
#' @importFrom tidyr pivot_longer
#' @inheritParams compute_emd
#' @param uncorrected Dataframe of uncorrected data
#' @param corrected Dataframe of corrected data
#' @param binSize The size of bins to use when binning data
#' @param plots If TRUE, a violin and scatter plot of the emds will be returned
#' @param filter_limit Limit for EMD removal (Removing EMDs that are below filter_limit in both before and after correction)
#' @family emd
#' @export
evaluate_emd <- function(uncorrected,
                         corrected,
                         binSize = 0.1,
                         cell_col = "label",
                         batch_col = "batch",
                         markers = NULL,
                         plots = TRUE,
                         filter_limit = 2){

  # Check for package
  missing_package("emdist", "CRAN")
  if(!plots) missing_package("plyr", "CRAN")

  # Get markers if not given
  if(is.null(markers)){
    markers <- uncorrected %>%
      cyCombine::get_markers()
  }

  # Check colnames
  check_colname(colnames(corrected), cell_col, "corrected set")
  check_colname(colnames(uncorrected), cell_col, "uncorrected set")
  check_colname(colnames(corrected), batch_col, "corrected set")
  check_colname(colnames(uncorrected), batch_col, "uncorrected set")
  check_colname(colnames(corrected), 'id', "corrected set")
  check_colname(colnames(uncorrected), 'id', "uncorrected set")

  # Define cell columns as characters to avoid problems with factors
  corrected[[cell_col]] <- corrected[[cell_col]] %>%
    as.character()
  uncorrected[[cell_col]] <- uncorrected[[cell_col]] %>%
    as.character()

  message("Computing EMD for corrected data..")
  emd_corrected <- corrected %>%
    dplyr::arrange(id) %>%
    cyCombine::compute_emd(binSize = binSize,
                           cell_col = cell_col,
                           batch_col = batch_col,
                           markers = markers)

  message("Computing EMD for uncorrected data..")
  emd_uncorrected <- uncorrected %>%
    dplyr::arrange(id) %>%
    cyCombine::compute_emd(binSize = binSize,
                           cell_col = cell_col,
                           batch_col = batch_col,
                           markers = markers)


  # Extracting EMD values
  unlist_cor <- emd_corrected %>%
    unlist()
  unlist_uncor <- emd_uncorrected %>%
    unlist()

  # Create a tibble based on the computed EMDs
  emds <- tibble::tibble(
    "Name" = names(unlist_cor),
    "Corrected" = unlist_cor,
    "Uncorrected" = unlist_uncor,
    "Reduction" = unlist_uncor - unlist_cor
  ) %>%
    dplyr::filter(!is.na(Reduction))

  # Apply filter
  message(paste("Removing EMDs below", filter_limit, "both before and after correction"))
  emds_filtered <- emds %>%
    dplyr::filter(!(Corrected < filter_limit & Uncorrected < filter_limit))


  # Calculate total reduction
  reduction <- (sum(emds_filtered$Reduction) / sum(emds_filtered$Uncorrected)) %>%
    round(2)
  message("The reduction is: ", reduction)



  if(!plots){
    return(list("reduction" = reduction,
                "emd" = emds))
  }

  message("Creating plots..")

  # Define the plotting limit
  limit <- emds_filtered$Uncorrected %>%
    max() %>%
    plyr::round_any(5, f = ceiling) + 1

  # Create violin plot
  violin <- emds_filtered %>%
    tidyr::pivot_longer(cols = dplyr::ends_with("orrected"), names_to = "corrected") %>%
    ggplot2::ggplot(ggplot2::aes(x = corrected, y = value)) +
    ggplot2::geom_violin() +
    ggplot2::geom_boxplot(width = 0.1) +
    ggplot2::labs(x = "",
         y = "Earth Mover's Distance",
         title = "Comparison of EMD before and after correction",
         subtitle = paste("Reduction: ", reduction))

  # Create scatterplot
  scatterplot <- emds %>%
    ggplot2::ggplot(ggplot2::aes(x = Corrected, y = Uncorrected)) +
    ggplot2::geom_point() +
    ggplot2::annotate("rect", xmin = 0, xmax = filter_limit, ymin = 0, ymax = filter_limit, alpha = .5) +
    ggplot2::coord_cartesian(xlim = c(-2, limit), ylim = c(-2, limit)) +
    ggplot2::geom_abline(slope = 1, intercept = 0) +
    ggplot2::labs(x = "EMD - Corrected",
         y = "EMD - Uncorrected",
         title = "",
         subtitle = paste("Reduction:", reduction))


  message("Evaluation complete.")
  return(list("violin" = violin,
              "scatterplot" = scatterplot,
              "reduction" = reduction,
              "emd" = emds))
}



# MAD ----

#' Compute MAD
#'
#' Given a dataframe, the Median Absolute Deviation (MAD) is calculated per-marker, per-batch
#'
#'
#' @param df Dataframe to compute the MADs of
#' @param cell_col Column name of df that contains cell population labels (or clusters)
#' @param batch_col Column name of df that contains batch numbers
#' @param markers Vector of the markers to calculate EMD for. If NULL, \code{\link{get_markers}} will be used to find markers
#' @family mad
#' @export
compute_mad <- function(df,
                        cell_col = "label",
                        batch_col = "batch",
                        markers = NULL){

  # Check colnames
  check_colname(colnames(df), cell_col, "df")
  check_colname(colnames(df), batch_col, "df")

  # Define cell columns as characters to avoid problems with factors
  df[[cell_col]] <- df[[cell_col]] %>%
    as.character()

  # Get markers if not given
  if(is.null(markers)){
    markers <- df %>%
      cyCombine::get_markers()
  }

  # Extract batches
  batches <- df %>%
    dplyr::pull(batch_col) %>%
    unique() %>%
    sort()

  # Extract cell types
  cellTypes <- df %>%
    dplyr::pull(cell_col) %>%
    unique() %>%
    sort()


  # Compute MADs in each batch per-marker, per-cluster
  mads <- list()
  for (cellType in cellTypes) {
    mads[[cellType]] <- list()

    for (b in batches) {
      # Calculate MAD per marker
      MAD <- df %>%
        dplyr::filter(.data[[cell_col]] == cellType,
                      .data[[batch_col]] == b) %>%
        dplyr::select(dplyr::all_of(markers)) %>%
        apply(2, stats::mad)

      mads[[cellType]][[b]] <- MAD
    }
  }

  return(mads)
}

#' Evaluate MAD
#'
#' The function computes the Median Absolute Deviation of the two given datasets.
#'   Then the reduction is calculated as the relative change in total MAD.
#'
#' @inheritParams compute_mad
#' @param uncorrected Dataframe of uncorrected data
#' @param corrected Dataframe of corrected data
#' @param filter_limit Limit for MAD removal (Removing MADs that are below or equal to filter_limit in both before and after correction). Default filters none of the values.
#' @family mad
#' @export
evaluate_mad <- function(uncorrected,
                         corrected,
                         cell_col = "label",
                         batch_col = "batch",
                         markers = NULL,
                         filter_limit = NULL){

  # Get markers if not given
  if(is.null(markers)){
    markers <- uncorrected %>%
      cyCombine::get_markers()
  }

  # Check colnames
  check_colname(colnames(corrected), cell_col, "corrected set")
  check_colname(colnames(uncorrected), cell_col, "uncorrected set")
  check_colname(colnames(corrected), batch_col, "corrected set")
  check_colname(colnames(uncorrected), batch_col, "uncorrected set")
  check_colname(colnames(corrected), 'id', "corrected set")
  check_colname(colnames(uncorrected), 'id', "uncorrected set")

  # Define cell columns as characters to avoid problems with factors
  corrected[[cell_col]] <- corrected[[cell_col]] %>%
    as.character()
  uncorrected[[cell_col]] <- uncorrected[[cell_col]] %>%
    as.character()


  message("Computing MAD for corrected data..")
  mad_corrected <- corrected %>%
    dplyr::arrange(id) %>%
    cyCombine::compute_mad(markers = markers,
                           batch_col = batch_col,
                           cell_col = cell_col)

  message("Computing MAD for uncorrected data..")
  mad_uncorrected <- uncorrected %>%
    dplyr::arrange(id) %>%
    cyCombine::compute_mad(markers = markers,
                           batch_col = batch_col,
                           cell_col = cell_col)


  # Extracting MAD values
  unlist_cor <- mad_corrected %>%
    unlist()
  unlist_uncor <- mad_uncorrected %>%
    unlist()

  # Create a tibble based on the computed MADs
  mads <- tibble::tibble(
    "Name" = names(unlist_cor),
    "Corrected" = unlist_cor,
    "Uncorrected" = unlist_uncor,
    "Difference" = abs(unlist_uncor - unlist_cor)
  )

  # Apply filter
  if (!is.null(filter_limit)) {
    message(paste("Removing MADs below or equal to", filter_limit, "both before and after correction"))
    mads_filtered <- mads %>%
      dplyr::filter(!(Corrected <= filter_limit & Uncorrected <= filter_limit))
  } else {
    mads_filtered <- mads
  }

  # Calculate combined MAD score (median of all aboslute differences between uncor/cor)
  score <- median(mads_filtered$Difference, na.rm = T) %>% round(2)

  message("The MAD score is: ", score)

  return(list("score" = score,
              "mad" = mads))

}
