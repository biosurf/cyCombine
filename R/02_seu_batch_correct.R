# Normalization ----

#' Batch-wise Normalization of Data Using Seurat
#'
#' This function normalizes the data in a batch-wise manner using Seurat.
#' The purpose is to minimize the impact of batch effects when clustering the data prior to batch correction.
#' Three normalization methods are implemented: Z-score, Rank, and Quantile normalization.
#' Z-score is recommended for batches from a single study/experiment.
#' Rank is recommended for data from different studies/experiments.
#' Quantile is not recommended.
#'
#' @inheritParams batch_correct_seurat
#' @param object A Seurat object
#' @param markers A vector of marker genes to normalize. Defaults to all genes if NULL.
#' @param norm_method Normalization method: "scale" (Z-score), "rank", or "qnorm" (Quantile normalization). Defaults to "scale".
#' @param ties.method Method for handling ties in rank normalization. Options are "average", "first", "last", "random", "max", or "min". Defaults to "average".
#' @param mc.cores Number of cores for parallelization
#'
#' @return A Seurat object with normalized data.
#' @export
normalize_seurat <- function(object,
                             markers = NULL,
                             layer = "counts",
                             norm_method = "scale",
                             ties.method = "average",
                             mc.cores = parallel::detectCores() - 1) {

  # Remove case-sensitivity
  norm_method <- tolower(norm_method)
  ties.method <- tolower(ties.method)

  # Error check
  if (norm_method == "rank" && !ties.method %in% c("average", "first", "last", "random", "max", "min")) {
    stop("When using norm_method = 'rank', please use an available ties.method (average, first, last, random, max, or min).")
  }

  if (is.null(markers)) {
    warning("All rows will be normalized. This can take a while.")
    # Get markers
    markers <- rownames(object)
  }

  # Messaging
  switch(norm_method,
         "rank" = message("Ranking expression data.."),
         "scale" = {
           message("Scaling expression data..")
           object <- Seurat::ScaleData(object, features = markers, split.by = "batch")
           return(object)
           },
         "CLR" = {
           message("CLR normalizing expression data..")
           object <- Seurat::NormalizeData(object, normalization.method = "CLR", margin = 2)
           return(object)
         },
         "qnorm" = {
           message("Quantile normalizing expression data..")
           object <- quantile_norm_seurat(object, markers = markers)
           return(object)
         },
         "none" = return(object),
         stop("Please use either 'scale', 'rank', or 'qnorm' as normalization method.")
  )

  # Rank at marker positions individually for every batch
  batches <- unique(object$batch)
  ranked_data <- pbmcapply::pbmclapply(
    setNames(batches, batches), function(b) {
      batch_cells <- SeuratObject::WhichCells(object, expression = batch == b)
      data <- as.matrix(SeuratObject::LayerData(object, layer)[markers, batch_cells])
      rowzeros <- rowSums(data) == 0
      # if (any(rowzeros)) {
      #   warning("A marker is 0 for an entire batch. This marker is removed.")
      #   data <- data[!rowzeros, ]
      # }
      data <- apply(data, 2, rank, ties.method = ties.method) / ncol(data)
    },
    mc.cores = mc.cores)
  ranked_data <- do.call(cbind, ranked_data)

  ranked_data <- ranked_data[, match(colnames(object), colnames(ranked_data))]



  SeuratObject::LayerData(object, "scale.data") <- ranked_data

  return(object)
}

# Batch-wise quantile normalization per marker using Seurat

quantile_norm_seurat <- function(object, markers = NULL, mc.cores = parallel::detectCores()) {
  message("Quantile normalizing expression data..")
  if (is.null(markers)) {
    markers <- rownames(object)
  }

  # Determine goal distributions for each marker by getting quantiles across all batches
  refq <- pbmcapply::pbmclapply(setNames(markers, markers), function(m) {
    quantile(SeuratObject::LayerData(object, "data")[m, ], probs = seq(0, 1, length.out = 5), names = FALSE)
  }, mc.cores = mc.cores)

  for (batch in unique(object$batch)) {
    batch_cells <- SeuratObject::WhichCells(object, ident = batch)
    data <- SeuratObject::LayerData(object, "data")[, batch_cells]
    for (m in markers) {
      qx <- quantile(data[m, ], probs = seq(0, 1, length.out = 5), names = FALSE)
      spf <- splinefun(x = qx, y = refq[[m]], method = "monoH.FC", ties = min)

      # Apply the spline function to adjust quantiles
      data[m, ] <- spf(data[m, ])
    }
    SeuratObject::LayerData(object, "scale.data") <- data
  }

  return(object)
}

# Clustering ----

#' Create Self-Organizing Map
#'
#' The function uses the kohonen package to create a Self-Organizing Map (SOM).
#' It is used to segregate the cells for batch correction to make the correction less affected
#' by samples with high abundances of a particular cell type.
#'
#' @inheritParams kohonen::supersom
#' @param object A Seurat object
#' @param markers A vector of marker genes to use for the SOM. Defaults to all genes if NULL.
#' @param seed The seed to use when creating the SOM. Defaults to 473.
#' @param xdim The x-dimension size of the SOM. Defaults to 8.
#' @param ydim The y-dimension size of the SOM. Defaults to 8.
#' @param rlen Number of times the data is presented to the SOM network. Defaults to 10.
#'
#' @return A vector of clustering labels
#' @export
create_som_seurat <- function(
    object,
    markers = NULL,
    seed = 473,
    rlen = 10,
    mode = c("online", "batch", "pbatch"),
    cluster_method = c("kohonen", "leiden"),
    resolution = 0.8,
    xdim = 8,
    ydim = 8) {

  cluster_method <- match.arg(cluster_method)
  mode <- match.arg(mode)
  # Default to all genes if markers are not specified
  if (is.null(markers)) {
    markers <- rownames(object)
  }

  if (cluster_method == "kohonen") {
    # Extract data for the markers
    data <- SeuratObject::LayerData(object, "scale.data")[markers, ]

    # SOM grid on overlapping markers, extract clustering per cell
    message("Creating SOM grid..")
    set.seed(seed)
    som_grid <- kohonen::somgrid(xdim = xdim, ydim = ydim, topo = "rectangular")
    som_model <- kohonen::som(
      t(data), grid = som_grid, rlen = rlen, dist.fcts = "euclidean", mode = mode)

    # Add labels to metadata
    object <- SeuratObject::AddMetaData(object, metadata = som_model$unit.classif, col.name = "Labels")
  } else if (cluster_method == "leiden") {
    object <- Seurat::FindClusters(object, method = "igraph", cluster.name = "Labels", resolution = resolution, random.seed = seed)
  }

  return(object)
}


# Batch correction ----



#' Correct data using ComBat
#'
#' Compute the batch correction on the data using the ComBat algorithm.
#' Define a covariate, either as a character vector or name of Seurat metadata column.
#' The covariate should preferably be the cell condition types but can be any column that infers heterogeneity in the data.
#' The function assumes that the batch information is in the "batch" column and the data contains a "sample" column with sample information.
#'
#' @inheritParams batch_correct_seurat
#' @param object A Seurat object.
#' @param markers A vector of marker genes to use for the correction. Defaults to all genes if NULL.
#' @param method The method for batch correction. Choose "ComBat" for cytometry data and "ComBat_seq" for bulk RNAseq data. Defaults to "ComBat".
#' @param covar The covariate ComBat should use. Can be a vector or a metadata column name in the input Seurat object. If NULL, no covar will be used.
#' @param anchor A column or vector specifying which samples are replicates and which are not. If specified, this column will be used as a covariate in ComBat. Be aware that it may be confounded with the condition.
#' @param ref.batch Optional. A string of the batch that should be used as the reference for batch adjustments.
#' @param parametric Logical. If TRUE, the parametric version of ComBat is used. If FALSE, the non-parametric version is used. Defaults to TRUE.
#' @param return_seurat Logical. Whether the matrix or a Seurat object should be returned.
#'
#' @return A Seurat object with corrected data.
#' @export
correct_data_seurat <- function(
    object,
    markers = NULL,
    mc.cores = 1,
    method = c("ComBat", "ComBat_seq"),
    covar = NULL,
    anchor = NULL,
    ref.batch = NULL,
    parametric = TRUE,
    return_seurat = TRUE) {


  method <- match.arg(method)
  message("Batch correcting data..")

  metadata <- object[[]]

  # Check for batch column
  if (!"batch" %in% colnames(metadata)) {
    stop("The 'batch' column is missing in the metadata.")
  }

  if (is.null(markers)) {
    markers <- rownames(object)
  }


  # Add covar to metadata if it's a vector
  if (!is.null(covar) && length(covar) > 1) {
    object$covar <- covar
    covar <- "covar"
  } else if (length(covar) == 1) {
    stopifnot("The covar column is missing" = covar %in% colnames(metadata))
  }

  # Add anchor to metadata if it's a vector
  if (!is.null(anchor) && length(anchor) > 1) {
    object$anchor <- anchor
    anchor <- "anchor"
  } else if (length(anchor) == 1) {
    stopifnot("The anchor column is missing" = anchor %in% colnames(metadata))
  }


  message("Batch correcting..")
  labels <- unique(object$Labels)
  corrected_data <- lapply(
    setNames(labels, labels),
    function(lab) {
      label_cells <- SeuratObject::WhichCells(object, expression = Labels == lab)
      object_lab <- object[markers, label_cells]
      correct_label_seurat(
        object_lab,
        covar = covar,
        anchor = anchor,
        parametric = parametric,
        ref.batch = ref.batch,
        method = method
      )}
  )

  corrected_data <- do.call(cbind, corrected_data)

  corrected_data <- corrected_data[, match(colnames(object), colnames(corrected_data))]

  if (!return_seurat) return(corrected_data)
  # Update the Seurat object with corrected data

  object[["cyCombine"]] <- SeuratObject::CreateAssayObject(data = corrected_data, key = "corrected.data_")


  return(object)
}

# Function to perform ComBat correction
combat_seurat <- function(object_lab, mod_matrix, parametric, ref.batch, method) {

  if (method == "ComBat") {
    data <- sva::ComBat(
      dat = as.matrix(SeuratObject::LayerData(object_lab, layer = "data")),
      batch = as.character(object_lab$batch),
      mod = mod_matrix,
      par.prior = parametric,
      ref.batch = ref.batch,
      prior.plots = FALSE
    )
  } else if (method == "ComBat_seq") {
    data <- sva::ComBat_seq(
      counts = as.matrix(SeuratObject::LayerData(object_lab, "counts")),
      batch = as.character(object_lab$batch),
      covar_mod = mod_matrix,
      full_mod = TRUE
    )
  }
  return(data)
}

# Function to correct each group
correct_label_seurat <- function(object_lab, covar, anchor, parametric, ref.batch, method) {
  num_covar <- 1
  num_anchor <- 1
  num_batches <- nlevels(factor(object_lab$batch))

  if (num_batches == 1) {
    message(paste("Label group", object_lab$Labels[1], "only contains cells from batch", object_lab$batch[1]))
    return(data)
  }

  if (!is.null(covar) && !check_confound(object_lab$batch, object_lab[[covar]][,1])) {
    num_covar <- nlevels(factor(object_lab[[covar]][,1]))
    if (num_covar == 1) covar <- NULL
  } else {
    covar <- NULL
  }

  if (!is.null(anchor) && !check_confound(object_lab$batch, object_lab[[anchor]][,1])) {
    num_anchor <- nlevels(factor(object_lab[[anchor]][,1]))
    if (num_anchor == 1) anchor <- NULL
  } else {
    anchor <- NULL
  }

  if (num_covar > 1 && num_anchor > 1 && check_confound(object_lab$batch, interaction(object_lab[[covar]][,1], object_lab[[anchor]][,1]))) {
    anchor <- NULL
  }

  if (!is.null(covar) && !is.null(anchor)) {
    mod_matrix <- model.matrix(~ object_lab[[covar]][,1] + object_lab[[anchor]][,1])
  } else if (!is.null(covar)) {
    mod_matrix <- model.matrix(~ object_lab[[covar]][,1])
  } else if (!is.null(anchor)) {
    mod_matrix <- model.matrix(~ object_lab[[anchor]][,1])
  } else {
    mod_matrix <- NULL
  }

  data <- combat_seurat(object_lab, mod_matrix, parametric, ref.batch, method)

  return(data)
}

# Function to check confounding
check_confound <- function(batch, covariate) {
  covariate_levels <- unique(covariate)
  for (level in covariate_levels) {
    if (length(unique(batch[covariate == level])) == 1) {
      return(TRUE)
    }
  }
  return(FALSE)
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
#' @param object A Seurat pbject
#' @param layer Layer to use from the Seurat object
#' @family batch
#' @importFrom sva ComBat ComBat_seq
#' @import stats
#' @examples
#' \dontrun{
#' corrected <- uncorrected %>%
#'   batch_correct(markers = markers,
#'   covar = "condition")
#'   }
#' @export
batch_correct_seurat <- function(
    object,
    xdim = 8,
    ydim = 8,
    rlen = 10,
    mode = c("online", "batch", "pbatch"),
    parametric = TRUE,
    method = c("ComBat", "ComBat_seq"),
    cluster_method = c("kohonen", "leiden"),
    resolution = 0.8,
    ref.batch = NULL,
    seed = 473,
    label = NULL,
    covar = NULL,
    anchor = NULL,
    markers = NULL,
    layer = "counts",
    norm_method = "scale",
    ties.method = "average",
    return_seurat = TRUE,
    mc.cores = parallel::detectCores() - 1) {

  cyCombine:::missing_package("pbmcapply")
  cyCombine:::missing_package("Seurat")

  stopifnot(
    "No 'batch' column in data." = "batch" %in% names(object[[]]))
  mode <- match.arg(mode)
  # scale_layer <- switch(
  #   norm_method,
  #   "scale" = "scale.data",
  #   "rank" = "scale.data",
  #   "none" = ifelse("scale.data" %in% Layers(object), "scale.data", "data"))

  if (is.null(markers)) {
    markers <- rownames(object)
  }

  object <- object[markers, ]

  for (i in seq_len(max(length(xdim), length(ydim)))) {
    xdim_i <- xdim[min(length(xdim), i)]
    ydim_i <- ydim[min(length(ydim), i)]

    message("Batch correcting using a SOM grid of dimensions ", xdim_i,"x", ydim_i)

    if (is(label, "NULL")) {
      # Create SOM on normalized data
      object <- normalize_seurat(
        object,
        markers = markers,
        layer = layer,
        norm_method = norm_method,
        ties.method = ties.method,
        mc.cores = mc.cores)
      # Remove excluded markers
      markers <- markers[markers %in% rownames(object)]
      object <- create_som_seurat(
        object,
        markers = markers,
        rlen = rlen,
        mode = mode,
        cluster_method = cluster_method,
        resolution = resolution,
        seed = seed,
        xdim = xdim_i,
        ydim = ydim_i)
      labels <- unique(object$Labels)
    } else {
      labels <- unique(object[[label]])
    }



    # Run batch correction

    corrected_data <- pbmcapply::pbmclapply(
      setNames(labels, labels),
      function(lab) {
        label_cells <- SeuratObject::WhichCells(object, expression = Labels == lab)
        object_lab <- object[, label_cells]

        correct_data_seurat(
          object_lab,
          covar = covar,
          anchor = anchor,
          markers = markers,
          parametric = parametric,
          method = method,
          ref.batch = ref.batch,
          mc.cores = 1,
          return_seurat = FALSE
        )
      },
      mc.cores = mc.cores
    )
    corrected_data <- do.call(cbind, corrected_data)

    if (!return_seurat) return(corrected_data)

    # Ensure data is in the same order as the original Seurat object
    flawed_cells <- !colnames(corrected_data) %in% colnames(object)
    if (sum(flawed_cells) > 0) {
      warning(sum(flawed_cells), " cell(s) were excluded in correction and are removed.")
      object <- object[, colnames(corrected_data)]
    }

    # Ensure data is in the same order as the original Seurat object
    corrected_data <- corrected_data[, match(colnames(object), colnames(corrected_data))]

    object[["cyCombine"]] <- SeuratObject::CreateAssayObject(data = corrected_data, key = "cycombine_")

    SeuratObject::DefaultAssay(object) <- "cyCombine"
  }

  message("Done!")
  return(object)
}



