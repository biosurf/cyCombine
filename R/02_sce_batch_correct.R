#' Batch-wise Normalization of Data Using Bioconductor
#'
#' This function normalizes the data in a batch-wise manner using functions
#' compatible with SingleCellExperiment objects. The purpose is to minimize
#' the impact of batch effects when clustering the data prior to batch correction.
#' Four normalization methods are implemented: Z-score ("scale"), Rank,
#' CLR ("clr"), and Log-normalization ("lognorm").
#'
#' Scale is recommended in most cases.
#' Rank sometimes performes better than scale.
#' CLR is recommended for ADT data (applied per feature across cells).
#' Log-normalization is a standard method for RNA-seq.
#'
#' @param sce A SingleCellExperiment object.
#' @inheritParams normalize
#' @param markers A character vector specifying the markers (features) to use for
#'   normalization. If NULL, all features are used.
#' @param assay Character string specifying the assay in \code{sce} to use
#'   as input (e.g., "counts", "logcounts"). Defaults to "counts".
#' @param norm_method Normalization method: "scale" (Z-score), "rank", "clr",
#'   "lognorm", or "qnorm" (Quantile normalization). Defaults to "scale".
#' @param norm_assay Character string for the name of the new assay where
#'   normalized data will be stored. Defaults to "normdata"..
#' @return A SingleCellExperiment object with a new assay containing the
#'   normalized data (named according to \code{new_assay_name}).
#' @param ... For compatibility
#'
normalize_sce <- function(
    sce,
    markers = NULL,
    assay = "exprs",
    norm_assay = "normalized",
    norm_method = c("scale", "rank", "CLR", "lognorm", "none"),
    ties.method = c("average", "first", "last", "random", "max", "min"),
    ...) {

  # Input validation
  if (!inherits(sce, "SummarizedExperiment")) {
    stop("Input must be a SingleCellExperiment object")
  }

  if (!"batch" %in% colnames(SummarizedExperiment::colData(sce))) {
    stop("Batch column 'batch' not found in colData")
  }

  if (!assay %in% SummarizedExperiment::assayNames(sce)) {
    stop("Assay '", assay, "' not found in the SingleCellExperiment object")
  }

  if (is.null(markers)) {
    # Get markers
    markers <- rownames(sce)
  }

  norm_method <- match.arg(norm_method)

  if (norm_method == "lognorm") {
    check_package("scuttle", "Bioc")
    message("Log-normalizing expression data..")
    sce <- scuttle::logNormCounts(sce)
    assay <- "lognorm"
    norm_method <- "scale"
  }

  # Messaging
  switch(
    norm_method,
    "rank" = message("Ranking expression data.."),
    "scale" = message("Scaling expression data.."),
    "CLR" = message("CLR normalizing expression data.."),
    "none" = return(sce)
  )


  # Get data
  expr_matrix <- object2mat(sce, assay = assay)

  # Initialize result matrix
  normalized_matrix <- expr_matrix

  for (batch in unique(sce$batch)) {
    batch_mask <- sce$batch == batch

    # Extract batch data
    batch_data <- expr_matrix[markers, batch_mask, drop = FALSE]

    # Apply transformation based on method
    if (norm_method == "scale") {
      row_means <- Matrix::rowMeans(batch_data)
      batch_data <- batch_data - row_means

      # Calculate row standard deviations
      row_sds <- sqrt(Matrix::rowMeans(batch_data^2))

      # Avoid division by zero
      row_sds[row_sds == 0] <- 1
      batch_data <- batch_data / row_sds

    } else if (norm_method == "rank") {
      # Rank normalization per feature within batch

      # Apply ranking row-wise (per feature)
      batch_data <- t(apply(batch_data, 1, function(x) {
        # Rank and then scale to [-1, 1]
        ranks <- rank(x, ties.method = ties.method)
        if (length(unique(ranks)) == 1) {
          # All values are the same
          rep(0, length(ranks))
        } else {
          2 * (ranks - 1) / (length(ranks) - 1) - 1
        }
      }))

    } else if (norm_method == "CLR") {
      # Centered log-ratio transformation per feature within batch
      pseudocount <- 1e-6
      batch_data <- batch_data + pseudocount

      # Calculate geometric mean for each feature (row)
      log_batch <- log(batch_data)
      row_geom_means <- Matrix::rowMeans(log_batch)

      # CLR: log(x) - mean(log(x)) for each feature
      batch_data <- log_batch - row_geom_means
    }

    normalized_matrix[markers, batch_mask] <- batch_data
  }

  # Add to SCE object
  SummarizedExperiment::assay(sce, norm_assay) <- normalized_matrix

  return(sce)
}

#' Create Self-Organizing Map or Other Clusters for Batch Correction
#'
#' This function clusters cells using various methods (SOM via Kohonen or FlowSOM,
#' graph-based Louvain/Leiden) on a specified assay (typically normalized data).
#' The resulting clusters are used to segregate cells for more targeted batch
#' correction, making the correction less affected by compositional differences
#' between batches.
#'
#' @param object A SingleCellExperiment object.
#' @param assay Character string specifying the assay containing the
#'   data to use for clustering (e.g., "normalized" from \code{normalize_sce} or "scale.data" from \code{normalize_seurat}).
#' @param markers A character vector specifying the markers (features) to use for
#'   clustering. If NULL, all features in the specified assay are used.
#' @param cluster_method Clustering method: "kohonen", "flowsom", "louvain", or "leiden".
#' @param xdim Integer, the x-dimension of the SOM grid (for "kohonen" and "flowsom").
#' @param ydim Integer, the y-dimension of the SOM grid (for "kohonen" and "flowsom").
#' @param rlen Integer, the number of iterations for SOM training.
#' @param mode Character, SOM training mode ("online", "batch", "pbatch") for "kohonen".
#' @param distf Character, distance function for SOM ("euclidean", "sumofsquares",
#'   "cosine", "manhattan", "chebyshev"). Note: "sumofsquares" is equivalent to
#'   Euclidean for Kohonen, FlowSOM uses numerical codes (see documentation).
#'   Cosine is often good for high-dimensional data.
#' @param nClus Integer, the target number of clusters for FlowSOM metaclusteirng.
#'   If NULL (default for FlowSOM), the grid node assignments are used directly.
#' @param resolution Numeric, resolution parameter for Louvain/Leiden community
#'   detection (\code{igraph}). Higher values lead to more clusters. Ignored
#'   if \code{cluster_method} is not "louvain" or "leiden".
#' @param seed Integer, random seed for reproducibility.
#' @param ... Arguments passed to `create_som`
#'
#' @return A SingleCellExperiment object with cluster labels added to \code{colData}.
#'
create_som_sce <- function(
    object,
    assay = ifelse(inherits(object, "Seurat"), "scale.data","normalized"),
    markers = NULL,
    seed = 473,
    rlen = 10,
    mode = c("online", "batch", "pbatch"),
    cluster_method = c("kohonen", "flowsom", "fusesom", "leiden"),
    resolution = 0.8,
    distf = c("euclidean", "sumofsquares", "cosine", "manhattan", "chebyshev"),
    xdim = 6,
    ydim = 6,
    nClus = NULL,
    ...) {

  mode <- match.arg(mode)
  distf <- match.arg(distf)
  cluster_method <- match.arg(cluster_method)
  # Default to all genes if markers are not specified
  if (is.null(markers)) {
    markers <- get_markers(object)
  }

  if (cluster_method == "leiden") {
    if (inherits(object, "SummarizedExperiment")) {
      check_package("SingleCellExperiment", "Bioc")
      check_package("scater", "Bioc")
      check_package("scran", "Bioc")
      check_package("igraph")

      object <- scater::runPCA(
        object,
        exprs_values = assay,
        subset_row = markers,
        ncomponents = ifelse(length(markers) > 100, 50, length(markers)/2),
        name = "pca_norm")
      snn_graph <- scran::buildSNNGraph(
        object,
        use.dimred = "pca_norm",
        type = "jaccard")

      clust <- igraph::cluster_leiden(snn_graph, objective_function = "modularity", resolution_parameter = resolution)
      label <- igraph::membership(clust)
    } else {
      check_package("Seurat")
      # if (!"pca" %in% Seurat::Reductions(object))
      object <- Seurat::RunPCA(object, npcs = ifelse(length(markers) > 100, 50, length(markers)/2), verbose = FALSE)
      object <- Seurat::FindNeighbors(
        object,
        reduction = "pca",
        compute.SNN = TRUE)

      object <- Seurat::FindClusters(
        object,
        random.seed = seed,
        algorithm = 4,
        resolution = resolution,
        random.seed = seed)
      label <- object$seurat_clusters
    }

  } else {
    label <- create_som(
      t(object2mat(object, assay = assay)),
      markers = markers,
      xdim = xdim, ydim = ydim,
      nClus = nClus,
      cluster_method = cluster_method,
      distf = distf,
      seed = seed,
      rlen = rlen,
      ...
    )
  }

  return(label)
}


#' Correct matrix using ComBat
#'
#' Compute the batch correction on the data using the ComBat algorithm.
#' Define a covariate, either as a character vector or name of metadata column.
#' The covariate should preferably be the cell condition types but can be any column that infers heterogeneity in the data.
#' The function assumes that the batch information is in the "batch" column.
#'
#' @inheritParams batch_correct
#' @param mat Data matrix
#' @param ... Arguments passed to `correct_label_mat`
#'
#' @return A SingleCellExperiment object with a new assay containing the batch-corrected data.
#'
correct_data_mat <- function(
    mat,
    metadata,
    covar = NULL,
    anchor = NULL,
    mc.cores = 1,
    pb = TRUE,
    ...) {

  APPLY <- set_apply(mc.cores, pb)
  message("Batch correcting..")

  # Check for batch column
  if (!"batch" %in% colnames(metadata)) {
    stop("The 'batch' column is missing in the metadata.")
  }

  # Add covar to metadata if it's a vector
  if (!is.null(covar) && length(covar) > 1) {
    metadata$covar <- covar
    covar <- "covar"
  } else if (length(covar) == 1) {
    stopifnot("The covar column is missing" = covar %in% colnames(metadata))
  }

  # Add anchor to metadata if it's a vector
  if (!is.null(anchor) && length(anchor) > 1) {
    metadata$anchor <- anchor
    anchor <- "anchor"
  } else if (length(anchor) == 1) {
    stopifnot("The anchor column is missing" = anchor %in% colnames(metadata))
  }


  labels <- unique(metadata$label)
  corrected_data <- APPLY(
    setNames(labels, labels),
    function(lab) {
      label_cells <- which(metadata$label == lab)
      correct_label_mat(
        mat[, label_cells],
        metadata = metadata[label_cells, ],
        ...
      )}
  )

  corrected_data <- do.call(cbind, corrected_data)
  # corrected_data <- corrected_data[, match(colnames(mat), colnames(corrected_data))]

  return(corrected_data)
}

# Function to correct each group
correct_label_mat <- function(
    mat,
    metadata,
    covar = NULL,
    anchor = NULL,
    parametric = TRUE,
    ref.batch = NULL,
    method = "ComBat",
    ...) {
  num_covar <- 1
  num_anchor <- 1
  num_batches <- nlevels(factor(metadata$batch))

  if (num_batches == 1) {
    message(paste("Label group", metadata$label[1], "only contains cells from batch", metadata$batch[1]))
    return(data)
  }

  if (!is.null(covar) && !check_confound(metadata$batch, metadata[[covar]][,1])) {
    num_covar <- nlevels(factor(metadata[[covar]][,1]))
    if (num_covar == 1) covar <- NULL
  } else {
    covar <- NULL
  }

  if (!is.null(anchor) && !check_confound(metadata$batch, metadata[[anchor]][,1])) {
    num_anchor <- nlevels(factor(metadata[[anchor]][,1]))
    if (num_anchor == 1) anchor <- NULL
  } else {
    anchor <- NULL
  }

  if (num_covar > 1 && num_anchor > 1 && check_confound(metadata$batch, interaction(metadata[[covar]][,1], metadata[[anchor]][,1]))) {
    anchor <- NULL
  }

  if (!is.null(covar) && !is.null(anchor)) {
    mod_matrix <- model.matrix(~ metadata[[covar]][,1] + metadata[[anchor]][,1])
  } else if (!is.null(covar)) {
    mod_matrix <- model.matrix(~ metadata[[covar]][,1])
  } else if (!is.null(anchor)) {
    mod_matrix <- model.matrix(~ metadata[[anchor]][,1])
  } else {
    mod_matrix <- NULL
  }

  data <- .combat(
    mat,
    batch = metadata$batch,
    mod_matrix = mod_matrix,
    parametric = parametric,
    ref.batch = ref.batch,
    method = method)

  return(data)
}



add_mat <- function(object, mat, assay = "cyCommbine", metadata = NULL) {

  if(inherits(object, "matrix") | inherits(object, "data.frame")) {
    return(mat)
  }

  # Ensure data is in the same order as the original Seurat object
  flawed_cells <- !colnames(mat) %in% colnames(object)
  if (sum(flawed_cells) > 0) {
    warning(sum(flawed_cells), " cell(s) were excluded in correction and are removed.")
    object <- object[, colnames(mat)]
    if (!is.null(metadata)) metadata <- metadata[colnames(mat),]
  }

  if (inherits(object, "Seurat")) {
    check_package("Seurat")
    # Update the Seurat object with corrected data
    object[[assay]] <- SeuratObject::CreateAssayObject(data = mat, key = "corrected_")
    SeuratObject::DefaultAssay(object) <- "cyCombine"
    if (!is.null(metadata)) object[[]] <- metadata
  } else if (inherits(object, "SummarizedExperiment")) {
    check_package("SingleCellExperiment", "Bioc")
    SummarizedExperiment::assay(object, assay) <- mat
    if (!is.null(metadata)) SummarizedExperiment::colData(object) <- metadata
  }
  return(object)
}



object2mat <- function(object, assay = "counts") {
  if (inherits(object, "Seurat")) {
    check_package("Seurat")
    mat <- as.matrix(SeuratObject::LayerData(object, layer = assay))
  } else if (inherits(object, "SummarizedExperiment")) {
    check_package("SingleCellExperiment", "Bioc")
    mat <- as.matrix(SummarizedExperiment::assay(object, assay))
  } else if(inherits(object, "matrix") | inherits(object, "data.frame")) {
    return(object)
  }
  colnames(mat) <- colnames(object)
  return(mat)
}
