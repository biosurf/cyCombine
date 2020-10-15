

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
compute_emd <- function(df, binSize = 0.1, cell_col = "label", batch_col = "batch"){
  markers <- df %>%
    dplyr::select_if(colnames(.) %!in% non_markers) %>%
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
        dplyr::filter(label == cellType,
                      batch == b) %>%
        dplyr::select(all_of(markers)) %>%
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
evaluate_emd <- function(preprocessed, corrected, cell_col = "label", batch_col = "batch"){

  cat("Computing emd for corrected data\n")
  emd_corrected <- corrected %>%
    dplyr::arrange(id) %>%
    cyCombine::compute_emd()

  cat("Computing emd for uncorrected data\n")
  emd_uncorrected <- preprocessed %>%
    dplyr::arrange(id) %>%
    cyCombine::compute_emd()

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


## EMD Without population ----

#' Compute EMD without population
#' @importFrom emdist emd2d
#' @export
compute_emd2 <- function(df, binSize = 0.1, batch_col = "batch"){
  # Define markers in dataframe
  markers <- df %>%
    dplyr::select_if(colnames(.) %!in% non_markers) %>%
    colnames()
  # Extract batches
  batches <- df %>%
    dplyr::pull(batch_col) %>%
    unique() %>%
    sort()

  # Create list of distribution matrices
  distr <- list()
  for (b in batches) {
    # Filter data on batch and bin
    distr[[b]] <- df %>%
      dplyr::filter(batch == b) %>%
      dplyr::select(dplyr::all_of(markers)) %>%
      apply(2, function(x) {
        bins <- seq(-1, 100, by = binSize)
        if (length(x) == 0) {
          # If no cells, fill with zeros
          rep(0, times = length(bins) - 1)
        }else{
          graphics::hist(x, breaks = bins,
                         plot = FALSE)$counts
        }
      })
  }

  # Compute emd from binned distributions
  distances <- list()
  for (marker in markers) {
    distances[[marker]] <- matrix(NA,
                                  nrow = length(batches),
                                  ncol = length(batches),
                                  dimnames = list(batches, batches))
    for (i in seq_along(batches)[-length(batches)]) {
      batch1 <- batches[i]
      for (j in seq(i + 1, length(batches))) {
        batch2 <- batches[j]
        A <- matrix(distr[[batch1]][, marker])
        B <- matrix(distr[[batch2]][, marker])
        distances[[marker]][batch1, batch2] <- emdist::emd2d(A, B)
      }
    }
  }
  # Extract max from each marker
  comparison <- matrix(NA, nrow = 1, ncol = length(markers),
                       dimnames = list("1", markers))
  for (marker in markers) {
    comparison[1, marker] <- max(distances[[marker]],
                                 na.rm = TRUE)
  }
  return(comparison)
}

#' Evaluate EMD without population
#' @importFrom tidyr pivot_longer
#' @export
evaluate_emd2 <- function(preprocessed, corrected, batch_col = "batch"){

  cat("Computing emd for corrected data\n")
  emd_corrected <- corrected %>%
    dplyr::arrange(id) %>%
    cyCombine::compute_emd2()

  cat("Computing emd for uncorrected data\n")
  emd_uncorrected <- preprocessed %>%
    dplyr::arrange(id) %>%
    cyCombine::compute_emd2()

  markers <- corrected %>%
    dplyr::select_if(colnames(.) %!in% non_markers) %>%
    colnames()
  cat("Computing reduction in emd\n")
  reduction <- matrix(NA, nrow = 1, ncol = length(markers),
                      dimnames = list("1", markers))
  for (marker in markers){
    emd_ori <- emd_uncorrected[1, marker]
    emd_cor <- emd_corrected[1, marker]
    if(emd_ori > 2){
      reduction[1, marker] <- (emd_ori - emd_cor) / emd_ori
    }else{
      reduction[1, marker] <- NA
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




## Silhouette scores ----

# a_i: average Euclidean between cell i and all other cells in the same group
# b_i: minimum of average distances between cell i and the cells in other groups not containing i
# s_i: b_i - a_i / max(b_i, a_i)
if(FALSE){
load("data/02_preprocessed_700k.Rdata")
load("data/02_corrected_700k.Rdata")

corrected <- corrected %>%
  dplyr::mutate(label = run_flowsom(., k = 10))
preprocessed <- preprocessed %>%
  dplyr::mutate(label = run_flowsom(., k = 10))

cor_sample <- corrected %>%
  slice_sample(n = 10000)
pre_sample <- preprocessed %>%
  slice_sample(n = 10000)

batches <- cor_sample %>%
  dplyr::pull(batch) %>%
  unique() %>%
  sort()


cors_s <- cor_sample %>%
  group_by(batch) %>%
  group_modify(function(df, ...){
    d <- df %>%
      dplyr::select(all_of(markers)) %>%
      dist() %>%
      as.matrix() %>%
      apply(1, mean)
    df <- df %>%
      dplyr::mutate(a_i = d)
    return(df)
  })
# euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
cor_sample %>%
  arrange(id) %>%
  dplyr::mutate(a_i = compute_ai(., group = "label")) %>%
  group_by(batch) %>%
  summarise(m = mean(a_i))


compute_ai <- function(df, group = "batch"){
  df <- df %>%
    group_by(.data[[group]]) %>%
    group_modify(function(df, ...){
      # Compute distances within each group
      d <- df %>%
        dplyr::select(all_of(markers)) %>%
        dist() %>%
        as.matrix() %>%
        apply(1, mean)
      df <- df %>%
        dplyr::mutate(a_i = d)
    return(df)
    }) %>%
    ungroup() %>%
    arrange(id)
  a_i <- df %>%
    pull(a_i)
  return(a_i)
}

compute_bi <- function(df, group = "batch"){
  for(i in nrow(df)){
    cell_i <- df[i, ]
    group <- cell_i %>%
      pull(group)

  }


  return(group)
}

cor_sample %>%
  mutate(test = compute_bi(.)) %>%
  select(batch, test)

for(group in batches){
  cell_group <- cor_sample %>%
    dplyr::filter(batch == group)

  d1 <- cell_group %>%
    dplyr::select(all_of(markers)) %>%
    dist() %>%
    as.matrix() %>%
    apply(1, mean)



}

for(i in nrow(cor_sample)){
  cell_i <- cor_sample[i, ]


  a_i <- cell_group %>%
    # dplyr::select(all_of(markers)) %>%
    dplyr::mutate(d = dist(select(., all_of(markers)), select(cell_i, all_of(markers))))
}

cor_sample %>%
  mutate(a_i = function(x){
    a_i <- x %>%
      dplyr::select(markers) %>%
      dplyr::filter()
      mean()
    return(a_i)
  })
}
