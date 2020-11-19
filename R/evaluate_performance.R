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
#' @export
compute_emd <- function(df,
                        binSize = 0.1,
                        cell_col = "label",
                        batch_col = "batch",
                        markers = NULL){

  # Check for package
  missing_package("emdist", "CRAN")

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
        dplyr::filter(label == cellType,
                      batch == b) %>%
        dplyr::select(all_of(markers)) %>%
        apply(2, function(x) {
          # Bin data
          bins <- seq(-5, 30, by = binSize)
          if (length(x) == 0) {
            rep(0, times = length(bins) - 1)
          }else{
            graphics::hist(x, breaks = bins,
                           plot = FALSE)$counts
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
          # if(sum(A) < 300 | sum(B) < 300){
          #   distances[[cellType]][[marker]][batch1, batch2] <- NA
          # } else{
          distances[[cellType]][[marker]][batch1, batch2] <- emdist::emd2d(A, B)
          # }

        }
      }
    }
  }

  # The comparison matrix contains the maximum distance between the batches for each cell type and marker
  comparison <- matrix(NA, nrow = length(cellTypes), ncol = length(markers),
                       dimnames = list(cellTypes, markers))
  for (cellType in cellTypes) {
    for (marker in markers) {
      comparison[cellType, marker] <- max(distances[[cellType]][[marker]],
                                          na.rm = TRUE)
    }
  }
  return(comparison)
}

#' Evaluate EMD
#'
#' The function computes the Eartch Movers Distance of the two given datasets.
#'   Then the reduction is calculated for every celltype
#' @importFrom tidyr pivot_longer
#' @inheritParams compute_emd
#' @param preprocessed Dataframe of uncorrected data
#' @param corrected Dataframe of corrected data
#' @export
evaluate_emd <- function(preprocessed,
                         corrected,
                         cell_col = "label",
                         batch_col = "batch",
                         markers = NULL){

  # Check for package
  missing_package("emdist", "CRAN")
  missing_package("viridis", "CRAN")
  missing_package("plyr", "CRAN")

  # Define cell columns as characters to avoid problems with factors
  corrected[[cell_col]] <- corrected[[cell_col]] %>%
    as.character()
  preprocessed[[cell_col]] <- preprocessed[[cell_col]] %>%
    as.character()

  message("Computing emd for corrected data")
  emd_corrected <- corrected %>%
    dplyr::arrange(id) %>%
    cyCombine::compute_emd()

  message("Computing emd for uncorrected data")
  emd_uncorrected <- preprocessed %>%
    dplyr::arrange(id) %>%
    cyCombine::compute_emd()

  # Extract cell types
  cellTypes <- corrected %>%
    dplyr::pull(cell_col) %>%
    unique() %>%
    sort()

  # Get markers
  if(is.null(markers)){
    markers <- corrected %>%
      get_markers()
  }
  message("Computing reduction in emd")
  reduction <- matrix(NA, nrow = length(cellTypes), ncol = length(markers),
                      dimnames = list(cellTypes, markers))
  for (cellType in cellTypes){
    for (marker in markers){
      emd_ori <- emd_uncorrected[cellType, ][marker]
      emd_cor <- emd_corrected[cellType, ][marker]
      if(emd_ori > 2){
        # Only compute reduction if there is a significant distance to be reduced (This avoids deviding by 0)
        reduction[cellType, marker] <- (emd_ori - emd_cor) / emd_ori
      }else{
        reduction[cellType, marker] <- NA
      }

    }
  }
  # Mean reducion (Perhaps a better aggregate function can be used)
  red <- mean(reduction, na.rm = TRUE) %>%
    round(2)

  message("Creating plots")
  # Dataframes containing the comparisons returned by compute_emd()
  scat_ori <- emd_uncorrected %>%
    tibble::as_tibble() %>%
    tidyr::pivot_longer(cols = all_of(colnames(.))) %>%
    dplyr::rename(scat_ori = value,
                  Marker = name)

  scat_cor <- emd_corrected %>%
    tibble::as_tibble() %>%
    tidyr::pivot_longer(cols = all_of(colnames(.))) %>%
    dplyr::rename(scat_cor = value,
                  Marker = name)


  # Combine dataframes
  scat <- scat_ori %>%
    dplyr::select(scat_ori) %>%
    dplyr::bind_cols(scat_cor, .name_repair = "minimal")

  # Define the plotting limit
  limit <- scat_ori$scat_ori %>%
    max() %>%
    plyr::round_any(5, f = ceiling) +1
    # ceiling() * 10

  plt <- scat %>%
    ggplot(aes(x = scat_cor, y = scat_ori, color = Marker, label = Marker)) +
    geom_text() +
    labs(x = "EMD - Corrected",
         y = "EMD - Uncorrected",
         title = "With cell population",
         subtitle = paste("Reduction:", red)) +
    coord_cartesian(xlim = c(0, limit), ylim = c(0, limit)) +
    viridis::scale_color_viridis(discrete = TRUE) +
    geom_abline(slope = 1, intercept = 0)


  message("Evaluation complete")
  return(list("plot" = plt,
              "reduction" = red,
              "emd_cor" = emd_corrected,
              "emd_uncor" = emd_uncorrected))
}

