

#' Evaluate LISI
#' @importFrom lisi compute_lisi
#' @export
evaluate_lisi <- function(df, batch_col = "batch", cell_col = "label", perplexity = 40){

  # Run PCA
  pca_df <- df %>%
    run_pca()

  # Compute LISI
  lisi_df <- lisi::compute_lisi(pca_df,
                                meta_data = tibble::as_tibble(df[, c(cell_col, batch_col)]),
                                label_colnames= c(batch_col, cell_col),
                                perplexity = perplexity)


  # Evaluate results
  iLisi <- lisi_df %>%
    pull(cell_col) %>%
    psych::harmonic.mean()
  cLisi <- lisi_df %>%
    pull(batch_col) %>%
    psych::harmonic.mean()

  score <- 2*(iLisi-1)*(cLisi)/(1-iLisi+cLisi)



  return(list("lisi" = lisi_df, "score" = score))
}

#' Compute EMD
#' @importFrom emdist emd2d
#' @export
compute_emd <- function(df, binSize = 0.1, non_markers, cell_col = "label", batch_col = "batch"){
  markers <- df %>%
    select_if(colnames(.) %!in% non_markers) %>%
    colnames()
  batches <- df %>%
    dplyr::pull(batch_col) %>%
    unique() %>%
    sort()
  cellTypes <- df %>%
    dplyr::pull(cell_col) %>%
    unique() %>%
    sort()

  distr <- list()
  for (b in batches) {
    # cat("Batch:", b, "\n")
    distr[[b]] <- list()
    for (cellType in cellTypes) {

      distr[[b]][[cellType]] <- df %>%
        filter(label == cellType,
               batch == b) %>%
        select(all_of(markers)) %>%
        apply(2, function(x) {
          bins <- seq(-1, 100, by = binSize)
          if (length(x) == 0) {
            rep(0, times = length(bins) - 1)
          }else{
            graphics::hist(x, breaks = bins,
                           plot = FALSE)$counts
          }

        })
    }
  }

  # bin_corrected <- corrected %>%
  #   group_by(batch, label) %>%
  #   tidyr::nest() %>%
  #   mutate(bin = purrr::map(data, ~function(x){
  #     binned <- x %>%
  #       select_if(colnames(.) %!in% non_markers) %>%
  #       apply(2, function(x) {
  #         graphics::hist(x, breaks = seq(-1, 50, by = binSize),
  #                        plot = FALSE)$counts
  #       })
  #     return(binned)
  #   }))
  #   select(all_of(markers)) %>%
  #   apply(2, function(x) {
  #     graphics::hist(x, breaks = seq(-1, 50, by = binSize),
  #     plot = FALSE)$counts
  #     })

  distances <- list()
  for (cellType in cellTypes) {
    distances[[cellType]] <- list()
    for (marker in markers) {
      distances[[cellType]][[marker]] <- matrix(NA, nrow = length(batches),
                                                ncol = length(batches), dimnames = list(batches,
                                                                                        batches))
      for (i in seq_along(batches)[-length(batches)]) {
        batch1 <- batches[i]
        for (j in seq(i + 1, length(batches))) {
          batch2 <- batches[j]
          A <- matrix(distr[[batch1]][[cellType]][,marker])
          B <- matrix(distr[[batch2]][[cellType]][,marker])
          distances[[cellType]][[marker]][batch1, batch2] <- emdist::emd2d(A, B)
        }
      }
    }
  }

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
#' @importFrom tidyr pivot_longer
#' @export
evaluate_emd <- function(preprocessed, corrected, cell_col = "label", batch_col = "batch", non_markers = c("batch", "sample", "covar", "som", "label", "id")){

  cat("Computing emd for corrected data\n")
  emd_corrected <- corrected %>%
    dplyr::arrange(id) %>%
    cyCombine::compute_emd(non_markers = non_markers)

  cat("Computing emd for uncorrected data\n")
  emd_uncorrected <- preprocessed %>%
    dplyr::arrange(id) %>%
    cyCombine::compute_emd(non_markers = non_markers)

  cellTypes <- corrected %>%
    dplyr::pull(cell_col) %>%
    unique() %>%
    sort()
  markers <- corrected %>%
    dplyr::select_if(colnames(.) %!in% non_markers) %>%
    colnames()
  cat("Computing reduction in emd\n")
  reduction <- matrix(NA, nrow = length(cellTypes), ncol = length(markers),
                      dimnames = list(cellTypes, markers))
  for (cellType in cellTypes){
    for (marker in markers){
      emd_ori <- emd_uncorrected[cellType, ][marker]
      emd_cor <- emd_corrected[cellType, ][marker]
      if(emd_ori > 2){
        reduction[cellType, marker] <- (emd_ori - emd_cor) / emd_ori
      }else{
        reduction[cellType, marker] <- NA
      }

    }
  }
  cat("Creating plots\n")
  scat_ori <- emd_uncorrected %>%
    tibble::as_tibble() %>%
    tidyr::pivot_longer(cols = all_of(colnames(.))) %>%
    dplyr::select(value)
  scat_cor <- emd_corrected %>%
    tibble::as_tibble() %>%
    tidyr::pivot_longer(cols = all_of(colnames(.))) %>%
    dplyr::select(value)

  scat <- scat_ori %>%
    dplyr::bind_cols(scat_cor, .name_repair = "minimal")
  colnames(scat) <- c("scat_ori", "scat_cor")

  plt <- scat %>%
    ggplot(aes(x = scat_cor, y = scat_ori)) +
    geom_point() +
    labs(x = "EMD - Corrected",
         y = "EMD - Uncorrected") +
    geom_abline(slope = 1, intercept = 0)

  red <- mean(reduction, na.rm = TRUE)
  cat("Evaluation complete\n")
  return(list("plot" = plt, "reduction" = red))

}
