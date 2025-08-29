# Normalization ----


#' Batch-wise normalization of data
#'
#' This function normalizes the data in a batch-wise manner.
#'   The purpose is to minimize the impact of batch correction when clustering the data prior to batch correction.
#'   Three normalisation methods are implemented: Z-score, Rank, and Quantile normalization.
#'   Z-score is recommended in cases where batches from a single study/experiment is merged.
#'   Rank is recommend in cases where data from different studies/experiments are merged.
#'   Quantile is not recommended.
#'
#' @param df tibble or seurat with expression values
#' @param markers Markers to normalize. If NULL, markers will be found using the \code{\link{get_markers}} function.
#' @param norm_method Normalization method. Should be either 'rank', 'scale', 'CLR', 'CLR_seu', 'CLR_med', or 'qnorm'. Default: 'scale'. CLR is recommended for ADT normalization. CLR_seu mimicks the Seurat implementation of CLR. CLR_med uses the median instead of geometric mean.
#' @param ties.method The method to handle ties, when using rank. Default: 'average'. See ?rank for other options.
#' @param mc.cores Number of cores for parallelization
#' @param pb Progress bar for parallelization
#' @param ... Arguments passed to `normalize_sce` or `normalize_seurat`
#' @family batch
#' @examples
#' \dontrun{
#' df_normed <- df %>%
#'   normalize()
#'   }
#' @export
normalize <- function(
    df,
    markers = NULL,
    norm_method = c("scale", "rank", "CLR_seu", "CLR_med", "CLR", "qnorm", "none"),
    ties.method = c("average", "first", "last", "random", "max", "min"),
    mc.cores = 1,
    pb = TRUE,
    ...) {
  if (!inherits(df, "data.frame")) {
    if (inherits(df, "SummarizedExperiment")) .normalize <- normalize_sce
    else if (inherits(df, "Seurat")) .normalize <- normalize_seurat
    return(
      .normalize(
        df,
        markers = markers,
        norm_method = norm_method,
        ties.method = ties.method,
        mc.cores = mc.cores,
        pb = pb,
        ...
      )
    )
  }

  # Remove case-sensitivity
  norm_method <- match.arg(norm_method)
  ties.method <- match.arg(ties.method)
  APPLY <- set_apply(mc.cores, pb)
  # Error check
  if (norm_method == "rank" && ties.method %!in% c("average", "first", "last", "random", "max", "min")) {
    stop("When using norm_method = 'rank', please use an available ties.method (average, first, last, random, max, or min).")
  }

  if (is.null(markers)) {
    markers <- df %>%
      cyCombine::get_markers()
  }

  # Messaging
  switch(
    norm_method,
    "rank" = message("Ranking expression data.."),
    "scale" = {
      message("Scaling expression data..")
      },
    "CLR" = {
      message("CLR normalizing expression data..")
      },
    "CLR_seu" = {
      message("CLR normalizing expression data with Seurat flavor..")
      },
    "CLR_med" = {
      message("CLR normalizing expression data using median..")
      },
    "qnorm" = {
      message("Quantile normalizing expression data..")
      df <- cyCombine:::quantile_norm(df, markers = markers)
      return(df)
      },
    "none" = {return(df)}
  )

  norm_f <- switch(
    norm_method,
    "rank" = function(values) rank(values, ties.method = ties.method) / length(values),
    "scale" = scale,
    "CLR" = clr_norm_mean,
    "CLR_seu" = clr_norm,
    "CLR_med" = clr_norm_med)

  if (!"id" %in% colnames(df)) df$id <- seq_len(nrow(df))
  # Scale or rank at marker positions individually for every batch
  df <- df %>%
    split(df$batch) %>%
    APPLY(function(df_batch) {
      df_batch[, markers] <- apply(df_batch[, markers], 2, function(values) {
        if (sum(values) == 0) return(values)
        return(norm_f(values))
      })
      return(df_batch)
    })
  df <- do.call(rbind, df)
  df <- dplyr::arrange(df, .data$id)
  return(df)
}


#' Batch-wise quantile normalization per marker
#'
#' This function quantile normalizes the data in a batch-wise manner.
#'   The purpose is to minimize the impact of batch correction when clustering the data prior to batch correction.
#'
#' @param df Dataframe with expression values
#' @param markers Markers to correct. If NULL, markers will be found using the \code{\link{get_markers}} function.
#' @family batch
#' @examples
#' \dontrun{
#' df_qnorm <- quantile_norm(preprocessed)
#'   }
quantile_norm <- function(df, markers = NULL) {
  message("Quantile normalizing expression data..")
  if (is.null(markers)) {
    # Get markers
    markers <- df %>%
      cyCombine::get_markers()
  }

  # Determine goal distributions for each marker by getting quantiles across all batches
  refq <- list()
  for (m in markers) {
    # Determine the quantiles
    refq[[m]] <- stats::quantile(unlist(df[,m]), probs=seq(0,1,length.out=5), names = F)
  }

  qnormed_expr <- df
  for (b in unique(df$batch)) {
    for (m in markers) {
      qx <- stats::quantile(unlist(df[df$batch == b, m]), probs=seq(0,1,length.out=5), names = F)
      spf <- stats::splinefun(x=qx, y=refq[[m]], method="monoH.FC", ties=min)

      # Apply the spline function to adjust quantiles
      qnormed_expr[qnormed_expr$batch == b, m] <- spf(unlist(df[df$batch==b, m]))

    }
  }

  return(qnormed_expr)
}

#' Centered Log-Ratio (CLR) Normalization
#'
#' This function was added to evaluate its performance in integrating across modalities.
#'
#' @param x A vector to normalize
#'
#' @return CLR-normalized values.
#' @noRd
clr_norm <- function(x) {
  geom_mean <- expm1(sum(log1p(x[x > 0]), na.rm = TRUE) / length(x))
  clr <- log1p(x / geom_mean)
  return(clr)
}
clr_norm_mean <- function(x) {
  geom_mean <- expm1(mean(log1p(x[x >= 0]), na.rm = TRUE))
  if (geom_mean == 0) return(x)
  clr <- log((x+1) / geom_mean)
  return(clr)
}
clr_norm_med <- function(x) {
  geom_median <- median(x, na.rm = TRUE)
  if (geom_median == 0) geom_median <- 1
  clr <- log((x+1) / geom_median)
  return(clr)
}

# Clustering ----

#' Create Self-Organizing Map
#'
#' The function uses the kohonen package to create a Self-Organizing Map.
#'  It is used to segregate the cells for the batch correction to make the correction less affected
#'  by samples with high abundances of a particular cell type.
#'
#' @inheritParams kohonen::supersom
#' @inheritParams normalize
#' @param seed The seed to use when creating the SOM.
#' @param xdim The x-dimension size of the SOM. If NULL, gridsize is computed using `FuseSOM::computeGridSize`
#' @param ydim The y-dimension size of the SOM. If NULL, gridsize is computed using `FuseSOM::computeGridSize`
#' @param rlen Number of times the data is presented to the SOM network
#' @param nClus (Usable with FlowSOM, FuseSOM, and kmeans) Number of clusters to export
#' @param cluster_method Cluster method to use. Defaults to kohonen. Options: c("kohonen", "flowsom", "fusesom", "kmeans")
#' @param distf Distance metric used for kohonen ("sumofsquares", "euclidean", "manhattan") or FlowSOM ("euclidean", "cosine", "manhattan", "chebyshev")
#' @param return_model Option to return the clustering model instead of the labels for downstream usage. See used package documentation for how to use the model. E.g., the kohonen model can be used with `predict(model, newdata = df[,markers])`
#' @param ... Arguments passed to `create_som_sce`
#' @family batch
#' @importFrom stats kmeans
#' @importFrom kohonen som somgrid
#' @examples
#' \dontrun{
#' labels <- create_som(uncorrected)
#'   }
#' @export
#' @return A vector of clustering labels
create_som <- function(
    df,
    markers = NULL,
    seed = 473,
    cluster_method = c("kohonen", "flowsom", "fusesom", "kmeans"),
    distf = c("euclidean", "sumofsquares", "cosine", "manhattan", "chebyshev"),
    rlen = 10,
    mode = c("online", "batch", "pbatch"),
    xdim = 8,
    ydim = 8,
    nClus = NULL,
    return_model = FALSE,
    ...) {
  mode <- match.arg(mode)
  distf <- match.arg(distf)
  cluster_method <- match.arg(cluster_method)

  if (!inherits(df, "data.frame") & !inherits(df, "matrix")) {
    return(
      create_som_sce(
        df,
        rlen = rlen,
        mode = mode,
        seed = seed,
        xdim = xdim,
        ydim = ydim,
        nClus = nClus,
        markers = markers,
        return_model = return_model,
        cluster_method = cluster_method,
        ...
      )
    )
  }



  if (is.null(markers)) {
    markers <- df  |>
      cyCombine::get_markers()
  }

  # SOM grid on overlapping markers, extract clustering per cell

  set.seed(seed)
  df <- as.matrix(df[, markers])

  if (is.null(xdim) | is.null(ydim)) {
    cyCombine:::check_package("FuseSOM", "Bioc")
    message("Computing gridsize")
    size <- FuseSOM::computeGridSize(df)
    if (size * size < 25) {
      size <- 5
    }
    xdim <- ydim <- size
    message("Using gridsize ", size)
  }

  labels <- switch(
    cluster_method,
    "kohonen" = {
      if (!distf %in% c("sumofsquares", "euclidean", "manhattan")) {
        warning("Distance function '", dist, "' not supported by Kohonen. Setting to 'sumofsquares'")
        distf <- "sumofsquares"
      }
      message("Creating SOM grid..")
      som_model <- kohonen::som(
        df,
        grid = kohonen::somgrid(xdim = xdim, ydim = ydim),
        rlen = rlen,
        mode = mode,
        dist.fcts = distf)
      if (return_model) return(som_model)
      som_model$unit.classif
    },

    "flowsom" = {
      cyCombine:::check_package("FlowSOM", "Bioc")
      if (distf == "sumofsquares") {
        warning("Distance function '", dist, "' not supported by FlowSOM. Setting to 'euclidean'")
        distf <- "euclidean"
      }
      distf <- switch(distf, "manhattan" = 1, "euclidean" = 2, "chebyshev" = 3, "cosine" = 4)

      if (!is.null(nClus)) {
        som_model <- FlowSOM::FlowSOM(
          df, xdim = xdim, ydim = ydim, nClus = nClus, distf = distf)
        if (return_model) return(som_model)
        FlowSOM::GetMetaclusters(som_model)
      } else {
        som_model <- FlowSOM::ReadInput(df) %>%
          FlowSOM::BuildSOM(xdim = xdim, ydim = ydim, distf = distf)
        if (return_model) return(som_model)
        FlowSOM::GetClusters(som_model)
      }
    },

    "kmeans" = {
      if (is.null(nClus)) nClus <- xdim*ydim
      model <- stats::kmeans(df, centers = nClus)
      if (return_model) return(model)
      model$cluster
    },
    "fusesom" = {
      cyCombine:::check_package("FuseSOM", "Bioc")

      som_model <- FuseSOM::generatePrototypes(df, size = round(sqrt(xdim*ydim)))
      if (return_model) return(som_model)
      if (is.null(nClus)) {
        som_model$classif
      } else {
        FuseSOM::clusterPrototypes(som_model, numClusters = nClus)
      }
    }
    )

  return(labels)
}


# Batch correction ----

#' Correct data using ComBat
#'
#' Compute the batch correction on the data using the ComBat algorithm.
#'  Define a covariate, either as a character vector or name of tibble column.
#'  The covariate should preferable be the cell condition types but can be any column that infers heterogeneity in the data.
#'  The function assumes that the batch information is in the "batch" column and the data contains a "sample" column with sample information.
#'
#' @inheritParams normalize
#' @param label The cluster or cell type label. Either as a column name or vector.
#' @param covar The covariate ComBat should use. Can be a vector or a column name in the input tibble.
#'   If NULL, no covar will be used
#' @param anchor Experimental: A column or vector specifying which samples are replicates and which are not. If specified, this column will be used as a covariate in ComBat. Be aware that it may be confounded with the condition.
#' @param parametric Default: TRUE. If TRUE, the parametric version of ComBat is used. If FALSE, the non-parametric version is used.
#' @param method Default: "ComBat". Choose "ComBat" for cytometry data and "ComBat_seq" for bulk RNAseq data.
#' @param ref.batch Optional. A string of the batch that should be used as the reference for batch adjustments.
#' @param ... Arguments passed to `correct_data_mat```
#' @family batch
#' @importFrom sva ComBat ComBat_seq
#' @examples
#' \dontrun{
#' corrected <- uncorrected %>%
#'   correct_data(label = labels, covar = "condition")
#'   }
#' @export
correct_data <- function(df,
                         label = NULL,
                         markers = NULL,
                         method = c("ComBat", "ComBat_seq"),
                         covar = NULL,
                         anchor = NULL,
                         ref.batch = NULL,
                         parametric = TRUE,
                         mc.cores = 1,
                         pb = FALSE,
                         ...) {
  if (!inherits(df, "data.frame")) {
    return(
      correct_data_mat(
        df,
        pb = pb,
        covar = covar,
        anchor = anchor,
        method = method,
        markers = markers,
        mc.cores = mc.cores,
        ref.batch = ref.batch,
        parametric = parametric,
        ...
      )
    )
  }


  method <- match.arg(method)

  APPLY <- set_apply(mc.cores, pb)
  message("Batch correcting data..")
  # Check for batch column
  cyCombine:::check_colname(colnames(df), "batch", "df")
  if (is.null(markers)) {
    # Get markers
    markers <- df |>
      cyCombine::get_markers()
  }

  # Add ID column to retain data order
  if("id" %!in% colnames(df)) df$id <- seq_len(nrow(df))

  # Add label to df
  if (length(label) == 1) {
    cyCombine:::check_colname(colnames(df), label, "df")
  } else {
    df$label <- label
    label <- "label"
  }

  # Add covar to df, if given
  if (!is.null(covar)) {
    if (length(covar) == 1) {
      cyCombine:::check_colname(colnames(df), covar, "df")
      df[[covar]] <- as.factor(df[[covar]])
    } else {
      # Covar was given as a vector
      df$covar <- as.factor(covar)
      covar <- "covar"
    }
    # Ensure there is more than 1 factor level
    if (nlevels(df[[covar]]) == 1) covar <- NULL
  }

  # Add anchor to df, if given
  if (!is.null(anchor)) {
    if (length(anchor) == 1) {
      cyCombine:::check_colname(colnames(df), anchor)
      df[[anchor]] <- as.factor(df[[anchor]])
    } else {
      # Anchor was given as a vector
      df$anchor <- as.factor(anchor)
      anchor <- "anchor"
    }
    # Ensure there is more than 1 factor level
    if (nlevels(df[[anchor]]) == 1) anchor <- NULL
  }

  corrected_data <- df |>
    split(df[[label]]) |>
    APPLY(function(df_label) {
      # Correct each label group with ComBat
      # Initiate anchor and covar counter
      num_covar <- 1
      num_anchor <- 1
      # Detect if only one batch is present in the node
      num_batches <- df_label$batch |>
        unique() |>
        length()
      lab <- df_label[[label]][1] # Current label group
      if (num_batches == 1) {
        batch <- df_label$batch[1]
        message(paste("Label group", lab, "only contains cells from batch", batch))
        return(df_label)
      }
      message(paste("Correcting Label group", lab))
      # Calculate number of covars in the node
      if (!is.null(covar)) {

        # Only use covar, if it does not confound with batch
        if (!cyCombine:::check_confound(df_label$batch, stats::model.matrix(~df_label[[covar]]))) {

          # If a node is heavily skewed to a single covar, it should be treated as having only 1 covar.
          # Get number of cells in the condition with most cells
          covar_counts <- table(df_label[[covar]])
          num_covar <- length(covar_counts)

          if (sum(covar_counts) < max(covar_counts) + num_covar*5) {
            message("The label group almost exclusively consists of cells from a single covar. Therefore, covar is ignored for this label group")
            num_covar <- 1
          }
        } else {
          message("Covar is confounded with batch. Ignoring covar in this label group")
        }
      }
      # Do a similar check on anchor
      if (!is.null(anchor)) {
        if (!cyCombine:::check_confound(df_label$batch, stats::model.matrix(~df_label[[anchor]]))) {
          # If a node is heavily skewed to a single anchor, it should be treated as having only 1 covar.
          # Get number of cells in the anchor with most cells
          anchor_counts <- table(df_label[[anchor]])
          num_anchor <- length(anchor_counts)

          if (sum(anchor_counts) < max(anchor_counts) + num_anchor*5) {
            message("The label group almost exclusively consists of cells from a single anchor group. Therefore, anchor is ignored for this label group")
            num_anchor <- 1
          }
        } else {
          message("Anchor is confounded with batch. Ignoring anchor in this label group")
        }
      }
      if (num_covar > 1 & num_anchor > 1) {
        # If neither covar nor anchor confounds with batch but they do each other, prioritise covar
        if (cyCombine:::check_confound(df_label$batch, stats::model.matrix(~df_label[[covar]] + df_label[[anchor]]))) {
          num_anchor <- 1
          message("Anchor and covar are confounded. Ignoring anchor in this label group")
        }
      }

      # Determine model
      if (num_covar > 1 & num_anchor == 1) {
        mod_matrix <- stats::model.matrix(~ df_label[[covar]])
      } else if (num_covar > 1 & num_anchor > 1) {
        mod_matrix <- stats::model.matrix(~df_label[[covar]] + df_label[[anchor]])
      } else if (num_covar == 1 & num_anchor > 1) {
        mod_matrix <- stats::model.matrix(~df_label[[anchor]])
      } else if (num_covar == 1 & num_anchor == 1) {
        mod_matrix <- NULL # No model matrix needed
      }

      # Compute ComBat correction
      ComBat_output <- df_label[, markers] %>%
        as.matrix() %>%
        t() %>%
        .combat(
          batch = df_label$batch,
          mod_matrix = mod_matrix,
          parametric = parametric,
          ref.batch = ref.batch,
          method = method) %>%
        t() %>%
        tibble::as_tibble() %>%
        dplyr::bind_cols(
          dplyr::select(df_label,
                        -dplyr::all_of(markers))) %>%
        # Cap values to range of input data
        dplyr::mutate(dplyr::across(dplyr::all_of(markers),
                                    function(x) {
                                      min <- min(df_label[[dplyr::cur_column()]])
                                      max <- max(df_label[[dplyr::cur_column()]])
                                      x <- ifelse(x < min, min, x)
                                      x <- ifelse(x > max, max, x)
                                      return(x)
                                    }))
      # Only add covar column, if it is not null
      if (!is.null(covar)) ComBat_output[[covar]] <- df_label[[covar]]
      # Only add anchor column, if it is not null
      if (!is.null(anchor)) ComBat_output[[anchor]] <- df_label[[anchor]]

      return(ComBat_output)
    })
  corrected_data <- do.call(rbind, corrected_data) |>
    dplyr::arrange(id) |>
    dplyr::relocate(id)  |>
    dplyr::mutate(batch = as.factor(batch))
  return(corrected_data)
}

#' Run combat
#' @noRd
.combat <- function(x, batch, mod_matrix, parametric, ref.batch, method, ...) {
  if (method == "ComBat") {
    x <- sva::ComBat(
      x,
      batch = as.character(batch), # The as.character is to remove factor levels not present in the SOM node
      mod = mod_matrix,
      par.prior = parametric,
      ref.batch = ref.batch,
      prior.plots = FALSE
    )
  } else if (method == "ComBat_seq") {
    x <- sva::ComBat_seq(
      x,
      batch = as.character(batch),
      covar_mod = mod_matrix,
      full_mod = TRUE
    )
  }
  return(x)
}

#' Alternate correction
#'
#' This function allows running ComBat with a custom covar mod matrix.
#'  A model could look like \code{stats::model.matrix(~df$covar+df$label)}
#'
#' @inheritParams correct_data
#' @param mod Covariate model to use in ComBat.
#' @examples
#' \dontrun{
#' corrected <- uncorrected %>%
#'   correct_data(mod = stats::model.matrix(~df$covar+df$label))
#'   }
correct_data_alt <- function(df,
                             mod,
                             markers = NULL,
                             parametric = TRUE) {
  message("Batch correcting data..")
  if (is.null(markers)){
    # Get markers
    markers <- df  |>
      cyCombine::get_markers()
  }


  corrected_data <- df[, markers] %>%
      # Compute ComBat correction
    t() %>%
    # The as.character is to remove factor levels not present in the SOM node
    sva::ComBat(batch = as.character(df$batch),
                mod = mod,
                par.prior = parametric) %>%
    t() %>%
    tibble::as_tibble() %>%
    dplyr::mutate(batch = df$batch,
                  sample = df$sample,
                  id = df$id) %>%
    # Cap values to range of input data
    dplyr::mutate(dplyr::across(dplyr::all_of(markers),
                                function(x) {
                                  min <- min(df[[dplyr::cur_column()]])
                                  max <- max(df[[dplyr::cur_column()]])
                                  x <- ifelse(x < min, min, x)
                                  x <- ifelse(x > max, max, x)
                                  return(x)
                                })) %>%
    dplyr::select(id, dplyr::everything())
  return(corrected_data)
}


# Wrapper ----

#' Run batch correction on data
#'
#' This is a wrapper function for the cyCombine batch correction workflow.
#'  To run the workflow manually, type "batch_correct" to see the source code of this wrapper and follow along or read the vignettes on the GitHub page \url{https://github.com/biosurf/cyCombine}.
#'
#' @inheritParams create_som
#' @inheritParams correct_data
#' @inheritParams normalize
#' @family batch
#' @examples
#' \dontrun{
#' corrected <- uncorrected %>%
#'   batch_correct(markers = markers,
#'   covar = "condition")
#'   }
batch_correct_df <- function(df,
                          label = NULL,
                          xdim = 8,
                          ydim = 8,
                          rlen = 10,
                          mode = c("online", "batch", "pbatch"),
                          parametric = TRUE,
                          method = c("ComBat", "ComBat_seq"),
                          cluster_method = c("kohonen", "flowsom", "fusesom", "kmeans"),
                          distf = c("euclidean", "sumofsquares", "cosine", "manhattan", "chebyshev"),
                          nClus = NULL,
                          mc.cores = 1,
                          pb = FALSE,
                          ref.batch = NULL,
                          seed = 473,
                          covar = NULL,
                          anchor = NULL,
                          markers = NULL,
                          norm_method = c("scale", "rank", "CLR_seu", "CLR_med", "CLR", "qnorm", "none"),
                          ties.method = "average",
                          ...) {

  if (!inherits(df, "data.frame")) stop("Please use 'batch_correct_generic()' on non-data.frame objects.")
  # A batch column is required
  cyCombine:::check_colname(colnames(df), "batch", "df")
  if (!is.null(markers)) lapply(markers, function(marker) cyCombine:::check_colname(colnames(df), marker, "df"))

  if (any(is.na(df$batch))) { # Check for NAs
    warning("Some batches contain NAs. These will be removed")
    df <- dplyr::filter(df, !is.na(batch))
  }
  mode <- match.arg(mode)
  method <- match.arg(method)
  norm_method <- match.arg(norm_method)
  cluster_method <- match.arg(cluster_method)

  for (i in seq_len(max(length(xdim), length(ydim)))) {
    xdim_i <- xdim[min(length(xdim), i)]
    ydim_i <- ydim[min(length(ydim), i)]
    nClus_i <- nClus[min(length(nClus), i)]

    message("Batch correcting using a SOM grid of dimensions ", xdim_i,"x", ydim_i)

    # Create SOM on scaled data
    label_i <- label
    if (is.null(label)) {
      df_norm <- cyCombine::normalize(
        df,
        markers = markers,
        norm_method = norm_method,
        ties.method = ties.method,
        mc.cores = mc.cores,
        pb = pb)
      label_i <- cyCombine::create_som(
        df_norm,
        markers = markers,
        rlen = rlen,
        mode = mode,
        seed = seed,
        xdim = xdim_i,
        ydim = ydim_i,
        distf = distf,
        cluster_method = cluster_method,
        nClus = nClus_i)
      rm(df_norm)
    }


    # Run batch correction
    df <- cyCombine::correct_data(
      df,
      label = label_i,
      covar = covar,
      anchor = anchor,
      markers = markers,
      parametric = parametric,
      mc.cores = mc.cores,
      pb = pb,
      method = method,
      ref.batch = ref.batch
      )
  }
  message("Done!")
  return(df)
}
