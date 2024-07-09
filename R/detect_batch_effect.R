#### Detect batch effects ----


#' Quicker function for detection of batch effects
#'
#' This function can be used to check if a dataset contains batch effects.
#'   The function employs three different approaches to detect the effects:
#'   1. The Earth Mover's Distance per marker when comparing each batch-batch pair.
#'   2. Density plots per marker, per batch for visual inspection.
#'   3. A MultiDimensional Scaling plot based on the median marker expression per sample.
#'   It can apply downsampling for a quicker analysis of larger datasets.
#'
#' @inheritParams detect_batch_effect
#' @param df Tibble containing the expression data and batch information. See prepare_data.
#' @param out_dir Directory for plot output
#' @param batch_col Name of column containing batch information
#' @param downsample Number of cells to include in detection. If not specified, all cells will be used.
#' @param seed Random seed for reproducibility
#' @family detect_batch_effect
#' @examples
#' \dontrun{
#' detect_batch_effect_express(df = exprs, out_dir = '/my/cycombine/dir/')
#' detect_batch_effect_express(df = exprs, downsample = 100000, out_dir = '/my/cycombine/dir/')
#' }
#' @export
detect_batch_effect_express <- function(df,
                                        out_dir,
                                        binSize = 0.1,
                                        markers = NULL,
                                        batch_col = "batch",
                                        downsample = NULL,
                                        seed = 472) {

  cyCombine:::missing_package("Matrix")
  cyCombine:::missing_package("ggridges", "CRAN")
  cyCombine:::missing_package("ggplot2", "CRAN")
  cyCombine:::missing_package("cowplot", "CRAN")
  cyCombine:::missing_package("grDevices", "CRAN")

  message('Starting the quick(er) detection of batch effects.')
  # Create output directory if missing
  cyCombine:::check_make_dir(out_dir)

  # Check batch_col and rename if necessary
  cyCombine:::check_colname(colnames(df), batch_col, location = "df")
  if (batch_col != 'batch') {
    df$batch <- df[, batch_col]
  }
  df$batch <- as.factor(df$batch)

  # Check multiple batches
  if (length(levels(df$batch)) < 2) {
    stop('Error! Please provide a datasets with more than one batch.')
  }

  # This works without clustering the data, so we set all labels to 1
  df$label <- 1

  # Check if IDs exist, otherwise add them
  if (!('id' %in% names(df))) {
    df$id <- 1:nrow(df)
  }

  # Downsample if requested
  if (!is.null(downsample)) {
    message(paste0('Downsampling to ', downsample, ' cells.'))
    if (downsample <= nrow(df)) {
      set.seed(seed)
      df <- df %>% dplyr::slice_sample(n = downsample)
    } else {
      warning('Specified downsample parameter exceeds the number of cells in the dataset.')
    }
  } # Perhaps only useful for some steps? E.g. the plotting of distributions or MDS plots



  ### Making distribution plots for all markers in each batch - good for diagnostics
  message('Making distribution plots for all markers in each batch.')

  if (is.null(markers)){
    markers <- df %>%
      cyCombine::get_markers()
  }

  # For each marker, make the plot
  grDevices::pdf(NULL)
  p <- list()
  for (c in 1:length(all_markers)) {

    p[[c]] <- ggplot2::ggplot(df, ggplot2::aes_string(x = all_markers[c],
                                             y = "batch")) +
      ggridges::geom_density_ridges(ggplot2::aes(color = batch, fill = batch),
                                    alpha = 0.4,
                                    quantile_lines = TRUE) +
      ggplot2::theme_bw()
  }

  # Save the plots
  suppressMessages(cowplot::save_plot(paste0(out_dir, '/distributions_per_batch.png'), cowplot::plot_grid(plotlist = p, nrow = round(length(all_markers) / 6)), base_width = length(all_markers) / (4/3), base_height = length(all_markers)))
  message(paste0('Saved marker distribution plots here: ', out_dir, '/distributions_per_batch.png.'))


  ### Use EMD-based batch effect detection
  message('Applying global EMD-based batch effect detection.')

  # Perform EMD calculations
  emd <- df %>%
    dplyr::arrange(id) %>%
    cyCombine::compute_emd(binSize = binSize)

  # Get summary per marker ACROSS batches
  emd_markers <- cbind.data.frame(sapply(emd[[1]], mean, na.rm = T), sapply(emd[[1]], stats::sd, na.rm = T))

  # Get mean PER BATCH per marker
  emd <- lapply(emd[[1]], function(x) {Matrix::forceSymmetric(x, uplo='U')}) # Contains the pairwise EMDs, now made symmetrical for easy extraction of values
  batch_means <- lapply(emd, function(x) {Matrix::colMeans(as.matrix(x), na.rm = T)})



  # # Flagging markers with observed high EMDs
  # flags_per_marker <- sapply(batch_means, function(x) {sum(x>4)}) # It's hard to set a numerical threshold - should it be data-based, i.e. IQR * 1.5, or?
  # flagged_markers <- names(which(flags_per_marker != 0))
  #
  # print(paste0('It looks like there are batch effects in the following markers: ', paste(flagged_markers, collapse = ', '), '.'))
  #
  # # Flagging markers with observed VERY high EMDs
  # extreme_flags_per_marker <- sapply(batch_means, function(x) {sum(x>6)})
  # extreme_flagged_markers <- names(which(extreme_flags_per_marker != 0))
  #
  # print(paste0('The effect is very strong in the following markers: ', paste(extreme_flagged_markers, collapse = ', '), '.'))
  #
  # for (m in extreme_flagged_markers) {
  #   print(paste0('The strongest effect is observed for ', m, ' in the batch: ', names(which.max(batch_means[[m]])), '.'))
  # }


  # Maybe we should make more of a scaling regarding how strong the effects are - relative to one another
  emd_markers <- cbind.data.frame(emd_markers, rownames(emd_markers)); colnames(emd_markers) <- c('mean', 'sd', 'marker')
  emd_markers$marker <- factor(emd_markers$marker, levels = emd_markers$marker[order(emd_markers$mean, decreasing=T)])

  p <- ggplot2::ggplot(emd_markers, ggplot2::aes(x = marker, y = mean)) +
    ggplot2::geom_bar(stat="identity", ggplot2::aes(fill = mean)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggplot2::ylab('Mean EMD') + ggplot2::xlab('') +
    ggplot2::geom_errorbar(ggplot2::aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                           position=ggplot2::position_dodge(.9))
  suppressMessages(ggplot2::ggsave(p, filename = paste0(out_dir, '/emd_per_marker.png')))

  message(paste0('Saved EMD plot here: ', out_dir, '/emd_per_marker.png.\n'))


  # The magnitude of expression also matters for the size EMD ?
  # marker_means <- df %>% select(get_markers(df)) %>% as.matrix() %>% matrixStats::colMeans()


  # Flagging markers with clear outlier batches with STRONG effects
  any_outliers <- F
  for (m in emd_markers$marker[emd_markers$mean > stats::median(emd_markers$mean)]) {
    if (any(batch_means[[m]] > (stats::quantile(batch_means[[m]], 0.75) + stats::IQR(batch_means[[m]])*3))) {
      found_outliers <- which(batch_means[[m]] > (stats::quantile(batch_means[[m]], 0.75) + stats::IQR(batch_means[[m]])*3))
      message(paste0(m, ' has clear outlier batch(es): ', paste(names(found_outliers), collapse = ', ')))

      summary_non <- df %>% dplyr::filter(!(batch %in% names(found_outliers))) %>% dplyr::pull(m) %>% summary()
      summary_out <- df %>% dplyr::filter(batch %in% names(found_outliers)) %>% dplyr::pull(m) %>% summary()

      message('Summary of the distribution in the OUTLIER batch(es):')
      message(paste(names(summary_out), '=', round(summary_out,2), collapse = ', '))
      cat('\n')

      message('Summary of the distribution in the non-outlier batches:')
      message(paste(names(summary_non), '=', round(summary_non,2), collapse = ', '))

      cat('\n\n')

      any_outliers <- T
    }
  }
  if (!(any_outliers)) {message('None of the markers has very strong outlier batches, consult plots for more general deviations.')}


  ### MDS plot for median protein expression per sample
  message('Making MDS plot for median protein expression per sample.')

  median_expr <- df %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise_at(cyCombine::get_markers(df), stats::median)

  dist_mat <- as.matrix(stats::dist(median_expr[, all_markers]))   # Euclidean distance
  rownames(dist_mat) <- colnames(dist_mat) <- median_expr$sample

  mds <- as.data.frame(stats::cmdscale(dist_mat, k = 2))
  mds$sample <- as.factor(rownames(mds))
  mds$batch <- as.factor(df$batch[match(mds$sample, df$sample)])

  p <- ggplot2::ggplot(mds, ggplot2::aes(V1, V2)) +
    ggplot2::geom_point(ggplot2::aes(colour = batch), size = 2) +
    ggplot2::labs(x = "MDS1", y = "MDS2", title = "MDS plot") +
    ggplot2::theme_bw()

  suppressMessages(ggplot2::ggsave(p, filename = paste0(out_dir, '/MDS_per_batch.png')))
  message(paste0('Saved MDS plot here: ', out_dir, '/MDS_per_batch.png'))

  message('Done!')
}



#' Full function for detection of batch effects using cluster proportions
#'
#' This function is used for batch effect detection in multidimensional datasets.
#'   The function applies a SOM-based clustering to a dataset in order to compare not only marker expression differences across batches,
#'   but also the cluster percentages in each batch to detect possible populations that are over-/under-represented in a single batch.
#'   This is coupled with UMAP plots to assist the interpretation of the results.
#'   However, this is primarily meaningful for sets with 3-30 batches - in cases outside this range, only the UMAPs will be generated.
#'
#' @inheritParams compute_emd
#' @param df Tibble containing the expression data and batch information. See prepare_data.
#' @param downsample Number of cells to include in detection. If not specified all cells will be used. One should be careful with the downsampling here as too strong downsampling leads to spurious results.
#' @param norm_method Normalization methods (options = 'scale' and 'rank')
#' @param xdim Grid size in x-axis for SOM (default = 8)
#' @param ydim Grid size in y-axis for SOM (default = 8)
#' @param seed Random seed for reproducibility
#' @param markers If only some markers should be used this parameter is used to define them. If not set, all markers are used.
#' @param batch_col Name of column containing batch information
#' @param label_col If existing labels should be used, this column must be present in the data
#' @param out_dir Directory for plot output
#' @param name Name of dataset - used for plot titles
#' @family detect_batch_effect
#' @examples
#' \dontrun{
#' detect_batch_effect(df = exprs)
#' detect_batch_effect(df = exprs, xdim = 8, ydim = 8, seed = 382,
#'                     markers = c('CD3', 'CD4', 'CD8a', 'CD20', 'CD19', 'CD56', 'CD33'))
#' }
#' @export
detect_batch_effect <- function(df,
                                out_dir,
                                downsample = NULL,
                                norm_method = "scale",
                                xdim = 8,
                                ydim = 8,
                                seed = 382,
                                markers = NULL,
                                binSize = 0.1,
                                batch_col = "batch",
                                label_col = "label",
                                name = 'raw data') {

  cyCombine:::missing_package("outliers")
  cyCombine:::missing_package("Matrix")
  cyCombine:::missing_package("ggplot2", "CRAN")
  cyCombine:::missing_package("cowplot", "CRAN")

  # Create output directory if missing
  cyCombine:::check_make_dir(out_dir)

  # Check batch_col and rename if necessary
  cyCombine:::check_colname(colnames(df), batch_col, location = "df")
  if (batch_col != 'batch') {
    df$batch <- df[, batch_col]
  }
  df$batch <- as.factor(df$batch)

  # Check multiple batches
  if (length(levels(df$batch)) < 2) {
    stop('Error! Please provide a datasets with more than one batch.')
  }

  # Check out dir
  if (is.null(out_dir)) {
    stop('Error! Please speicify output directory.')
  }

  # Downsample if requested
  if (!is.null(downsample)) {
    message(paste0('Downsampling to ', downsample, ' cells.'))
    if (downsample <= nrow(df)) {
      set.seed(seed)
      df <- df %>% dplyr::slice_sample(n = downsample)
    } else {
      warning('Specified downsample parameter exceeds the number of cells in the dataset.')
    }
  }

  if (is.null(markers)){
    markers <- df %>%
      cyCombine::get_markers()
  }

  # Create SOM on scaled data if no labels exist
  if (!(label_col %in% colnames(df))) {
    message('Determining new cell type labels using SOM:')
    labels <- df %>%
      cyCombine::normalize(markers = markers,
                           norm_method = norm_method) %>%
      cyCombine::create_som(markers = markers,
                            seed = seed,
                            xdim = xdim,
                            ydim = ydim)

    # Get SOM output
    df$label <- labels

  } else {
    message('Using existig cell type labels.')
    df$label <- df[[label_col]]

  }

  # Perform EMD and cluster percentage calculations - and filter to get outliers (only for sets with 3-30 batches)
  if (length(levels(df$batch)) >= 3 & length(levels(df$batch)) <= 30) {
    emd <- df %>%
      dplyr::mutate(label = as.character(label)) %>%
      cyCombine::compute_emd(binSize = binSize)

    markers_emd <- list()
    # Looking through EMDs to find culprits - using loops
    for (m in markers) {

      # Calculation of the mean per batch across all clusters for each marker
      marker_emd <- lapply(emd, function(x) {x[[m]]})
      marker_emd <- lapply(marker_emd, function(x) {Matrix::forceSymmetric(x, uplo='U')})

      marker_emd <- do.call(rbind, marker_emd)
      marker_emd <- Matrix::colMeans(as.matrix(marker_emd), na.rm = T)

      markers_emd[[m]] <- marker_emd
    }

    marker_emd <- which(stats::p.adjust(sapply(markers_emd, function(x) {outliers::dixon.test(x, opposite=F)$p.value}), method = 'BH') < 0.05 | stats::p.adjust(sapply(markers_emd,  function(x) {outliers::dixon.test(x, opposite=T)$p.value}), method = 'BH') < 0.05)
    message(paste0('\nThere are ', length(marker_emd), ' markers that appear to be outliers in a single batch:'))
    message(paste(markers[marker_emd], collapse = ', '))

    # Look for cluster over- and under-representation
    counts <- table(df$batch, df$label)
    perc <- (counts / Matrix::rowSums(counts))*100

    # use the Dixon test to look for outliers (test whether a single low or high value is an outlier)
    out_cl <- which(stats::p.adjust(apply(perc, 2, function(x) {outliers::dixon.test(x, opposite=F)$p.value}), method = 'BH') < 0.05 | stats::p.adjust(apply(perc, 2, function(x) {outliers::dixon.test(x, opposite=T)$p.value}), method = 'BH') < 0.05)

    message(paste0('\nThere are ', length(out_cl), ' clusters, in which a single cluster is strongly over- or underrepresented.'))

    for (cl in out_cl) {
      message(paste0('The cluster percentages for each batch in cluster ', cl, ' are:'))
      message(paste(names(perc[,cl]), '=', round(perc[,cl], 2), '%', collapse = ', '))

      exp_markers <- df %>%
        dplyr::filter(label == cl) %>%
        dplyr::summarise(dplyr::across(cyCombine::get_markers(df), stats::median)) %>%
        unlist() %>%
        sort(decreasing = T)

      message(paste0('The cluster expresses ', paste(names(exp_markers[which(exp_markers > 1)]), collapse = ', ')), '.')
      cat('\n')
    }
  }


  # Making a set of UMAP plots to help the diagnostics
  message('Making UMAP plots for up to 50,000 cells.')

  # SHOULD WE DOWNSAMPLE??
  if (nrow(df) > 50000) {
    set.seed(seed)
    df <- df %>% dplyr::slice_sample(n = 50000)
  }


  cyCombine::plot_dimred_full(df,
                              name,
                              type = "umap",
                              markers = NULL,
                              seed = 473,
                              out_dir)

  message(paste0('Saved UMAP plot for batches and labels here: ', out_dir, ' as UMAP_batches_labels.png.'))
  message(paste0('Saved UMAP plot colored by each marker in directory: ', out_dir, '/UMAP_markers.\n'))

  message('Done!')
}

