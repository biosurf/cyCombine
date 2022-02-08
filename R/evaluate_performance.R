# EMD ----

#' Compute EMD
#'
#' Given a dataframe, the Earth Movers Distance is computed, using the emdist package.
#'   A comparison matrix is returned containing the maximum distance for each marker for each group of cells.
#'
#' @param df Dataframe to compute the EMD of
#' @param binSize The size of bins to use when binning data
#' @param cell_col Column name of df that contains cell population labels (or clusters)
#' @param batch_col Column name of df that contains batch numbers
#' @param markers Vector of the markers to calculate EMD for. If NULL, \code{\link{get_markers}} will be used to find markers
#' @importFrom graphics hist
#' @export
compute_emd <- function(df,
                        binSize = 0.1,
                        cell_col = "label",
                        batch_col = "batch",
                        markers = NULL){

  # Check for package
  missing_package("emdist", "CRAN")
  missing_package("graphics", "CRAN")

  # Get markers if not given
  if(is.null(markers)){
    markers <- df %>%
      get_markers()
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

  # Create list of distribution matrices
  distr <- list()
  for (b in batches) {
    distr[[b]] <- list()
    for (cellType in cellTypes) {
      # Filter data on batch and cell type
      distr[[b]][[cellType]] <- df %>%
        dplyr::filter(.data[[cell_col]] == cellType,
                      batch == b) %>%
        dplyr::select(all_of(markers)) %>%
        apply(2, function(x) {
          # Bin data
          bins <- seq(-10, 30, by = binSize)
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

#' Evaluate EMD
#'
#' The function computes the Earth Mover's Distance of the two given datasets.
#'   Then the reduction is calculated as the relative change in total EMD.
#' @importFrom tidyr pivot_longer
#' @inheritParams compute_emd
#' @param preprocessed Dataframe of uncorrected data
#' @param corrected Dataframe of corrected data
#' @param plots If TRUE, a boxplot and scatterplot of the emds will be returned
#' @export
evaluate_emd <- function(preprocessed,
                         corrected,
                         binSize = 0.1,
                         cell_col = "label",
                         batch_col = "batch",
                         markers = NULL,
                         plots = TRUE,
                         filter_limit = 2){

  # Check for package
  missing_package("emdist", "CRAN")
  missing_package("viridis", "CRAN")
  missing_package("plyr", "CRAN")

  check_colname(colnames(corrected), cell_col, "corrected set")
  check_colname(colnames(preprocessed), cell_col, "uncorrected set")
  # Define cell columns as characters to avoid problems with factors
  corrected[[cell_col]] <- corrected[[cell_col]] %>%
    as.character()
  preprocessed[[cell_col]] <- preprocessed[[cell_col]] %>%
    as.character()

  message("Computing EMD for corrected data..")
  emd_corrected <- corrected %>%
    dplyr::arrange(id) %>%
    cyCombine::compute_emd(binSize = binSize,
                           cell_col = cell_col,
                           markers = markers)

  message("Computing EMD for uncorrected data..")
  emd_uncorrected <- preprocessed %>%
    dplyr::arrange(id) %>%
    cyCombine::compute_emd(binSize = binSize,
                           cell_col = cell_col,
                           markers = markers)


  # Extractin EMD values
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

  # Calculate total reduction
  message(paste("Removing EMDs below", filter_limit, "both before and after correction"))
  emds_filtered <- emds %>%
    dplyr::filter(!(Corrected < filter_limit & Uncorrected < filter_limit))

  # Compute reduction
  reduction <- (sum(emds_filtered$Reduction) / sum(emds_filtered$Uncorrected)) %>%
    round(2)
  message("The reduction is: ", reduction)

  if(plots == FALSE){
    return(list("reduction" = reduction,
                "emd" = emds))
  }

  message("Creating plots..")

  # Extract cell types
  cellTypes <- corrected %>%
    dplyr::pull(cell_col) %>%
    unique() %>%
    sort()

  # Get markers
  if(is.null(markers)){
    markers <- corrected %>%
      cyCombine::get_markers()
  }


  # Define the plotting limit
  limit <- emds_filtered$Uncorrected %>%
    max() %>%
    plyr::round_any(5, f = ceiling) + 1

  # Create violin plots
  violin <- emds_filtered %>%
    tidyr::pivot_longer(cols = ends_with("orrected"), names_to = "corrected") %>%
    ggplot(aes(x = corrected, y = value)) +
    geom_violin() +
    geom_boxplot(width = 0.1) +
    labs(x = "",
         y = "Earth Mover's Distance",
         title = "Comparison of EMD before and after correction",
         subtitle = paste("Reduction: ", reduction))


  scatterplot <- emds %>%
    ggplot(aes(x = Corrected, y = Uncorrected)) +
    geom_point() +
    annotate("rect", xmin = 0, xmax = filter_limit, ymin = 0, ymax = filter_limit, alpha = .5) +
    coord_cartesian(xlim = c(-2, limit), ylim = c(-2, limit)) +
    geom_abline(slope = 1, intercept = 0) +
    labs(x = "EMD - Corrected",
         y = "EMD - Uncorrected",
         title = "",
         subtitle = paste("Reduction:", reduction))


  message("Evaluation complete.")
  return(list("violin" = violin,
              "scatterplot" = scatterplot,
              "reduction" = reduction,
              "emd" = emds))
}



