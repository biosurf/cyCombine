#' Run the cyCombine workflow
#'
#' This is an alias for the `batch_correct` function.
#' See `?batch_correct` for details.
#' @param object An object of type "data.frame", "SingleCellExperiment", or "Seurat".
#' @export
cyCombine <- function(object, ...) {
  batch_correct(object, ...)
}


#' Run cyCombine batch correction
#'
#' As of cyCombine version 0.3.0, this function will run the cyCombine workflow
#' on any object of type "data.frame", "SingleCellExperiment", and "Seurat".
#'
#' @details
#' The workflow proceeds as follows:
#' 1.  **Normalization:** Normalizes data batch-wise using the method
#'     specified by \code{norm_method}.
#' 2.  **Clustering:** Clusters cells based on the normalized data using the chosen \code{cluster_method}.
#'     This step is skipped if \code{label} is provided and already exists in the object.
#' 3.  **Batch Correction:** Applies ComBat or ComBat-seq batch correction
#'     within each cluster defined in step 2.
#'
#' The process can be iterated by providing vectors for \code{xdim} and \code{ydim},
#' performing clustering and correction iteratively. In subsequent iterations, the
#' previously corrected data is used as input for normalization/clustering.
#'
#' @param object An object of type "data.frame", "SingleCellExperiment", or "Seurat".
#' @param markers A character vector specifying the markers (features) to use for
#'   normalization, clustering, and correction. If NULL, all features are used.
#' @param batch_col Character string specifying the column name in \code{colData(sce)}
#'   that contains batch information. Defaults to "batch".
#' @param xdim Integer or vector of integers, the x-dimension(s) of the SOM grid
#'   (for "kohonen" and "flowsom"). If a vector, the workflow iterates.
#' @param ydim Integer or vector of integers, the y-dimension(s) of the SOM grid.
#'   If a vector, the workflow iterates. Length should match \code{xdim} or be 1.
#' @param rlen Integer, the number of iterations for SOM training.
#' @param mode Character, SOM training mode ("online", "batch", "pbatch") for "kohonen".
#' @param cluster_method Clustering method: "kohonen", "flowsom", "fusesom", "leiden", or "kmeans.
#' @param resolution Numeric, resolution parameter for Louvain/Leiden clustering.
#' @param k_neighbors Integer, number of neighbors for k-NN graph construction.
#' @param distf Character, distance function for SOM.
#' @param nClus Integer, target number of clusters for FlowSOM metaclustering.
#' @param seed Integer, random seed for reproducibility.
#' @param assay Name of SingleCellExperiment assay to use
#' @param layer Name of Seurat layer to use.
#' @param norm_method Normalization method ("scale", "rank", "clr", "lognorm", "qnorm", "none").
#'    Applied before clustering unless \code{labels_colname} is provided. Defaults to "scale".
#' @param ties.method Method for handling ties in rank normalization.
#' @param label Character string for the cluster labels column. If this column
#'    already exists in \code{colData(sce)}, the normalization and clustering steps are
#'    skipped, and correction proceeds using these existing labels. Defaults to "Labels".
#' @param method Batch correction method ("ComBat", "ComBat_seq").
#' @param assay_combat_input Assay to use as input for ComBat (e.g., "logcounts", "normdata").
#'    Defaults to "logcounts". Ensure this assay exists if method="ComBat".
#' @param corrected_assay_name Name for the assay where the final corrected data is stored.
#'    Defaults to "corrected".
#' @param covar Optional: Covariate column name or vector for correction.
#' @param anchor Optional: Anchor column name or vector for correction.
#' @param ref.batch Optional: Reference batch for ComBat.
#' @param parametric Logical, for ComBat's parametric adjustment.
#' @param mc.cores Number of cores for parallel processing. Defaults to 1.
#' @param pb Logical, whether to show progress bars. Defaults to FALSE.
#'
#' @return An object of the same type as the input
#'
#' @export
batch_correct <- function(
    object,
    markers = NULL,
    metadata = NULL,
    label = NULL,
    xdim = 8,
    ydim = 8,
    rlen = 10,
    mode = c("online", "batch", "pbatch"),
    cluster_method = c("kohonen", "flowsom", "fusesom", "leiden", "kmeans"),
    resolution = 0.8,
    distf = c("euclidean", "sumofsquares", "cosine", "manhattan", "chebyshev"),
    nClus = NULL,
    seed = 473,
    assay = "exprs",
    layer = "data",
    norm_method = c("scale", "rank", "clr", "lognorm", "none"),
    ties.method = c("average", "first", "last", "random", "max", "min"),
    method = c("ComBat", "ComBat_seq"),
    covar = NULL,
    anchor = NULL,
    ref.batch = NULL,
    parametric = TRUE,
    mc.cores = 1,
    pb = FALSE
) {

  mode <- match.arg(mode)
  method <- match.arg(method)
  norm_method <- match.arg(norm_method)
  ties.method <- match.arg(ties.method)
  cluster_method <- match.arg(cluster_method)
  assay <- ifelse(method == "ComBat", assay, "counts")

  if (inherits(object, "matrix")) {
    stopifnot("No metadata provided" = is.null(metadata))
    if (any("CD" %in% rownames(object))) object <- t(object)
    object <- cbind(as.data.frame(object), metadata)
  }

  if (inherits(object, "Seurat")) {
    check_package("Seurat")
    stopifnot(
      "No 'batch' column in data." = "batch" %in% names(object[[]]))
    assay <- layer
    .normalize <- normalize_seurat
    .create_som <- create_som_sce
    .correct_data <- correct_data_mat

    metadata <- object[[]]
  } else if (inherits(object, "SummarizedExperiment")) {
    check_package("SingleCellExperiment", "Bioc")
    if (!"batch" %in% colnames(SummarizedExperiment::colData(sce)))
      stop("Batch column 'batch' not found in colData(sce).")
    .normalize <- normalize_sce
    .create_som <- create_som_sce
    .correct_data <- correct_data_mat

    metadata <- SummarizedExperiment::colData(sce)
  } else if (inherits(object, "data.frame")) {

    .normalize <- normalize
    .create_som <- create_som
    .correct_data <- correct_data
    metadata <- dplyr::select(object, dplyr::any_of(non_markers))

    # corrected <- batch_correct_df(
    #   object,
    #   pb = pb,
    #   xdim = xdim,
    #   ydim = ydim,
    #   rlen = rlen,
    #   mode = mode,
    #   seed = seed,
    #   label = label,
    #   distf = distf,
    #   nClus = nClus,
    #   covar = covar,
    #   method = method,
    #   anchor = anchor,
    #   markers = markers,
    #   mc.cores = mc.cores,
    #   ref.batch = ref.batch,
    #   parametric = parametric,
    #   norm_method = norm_method,
    #   ties.method = ties.method,
    #   cluster_method = cluster_method
    # )
    # return(corrected)
  } else if (!inherits(object, "matrix")) {
    stop(
      "Objects of type '" , class(object),
      "' are not yet supported. See '?batch_correct' for more details.")
    }


  if (is.null(markers)) {
    markers <- get_markers(object)
  }


  for (i in seq_len(max(length(xdim), length(ydim)))) {
    xdim_i <- xdim[min(length(xdim), i)]
    ydim_i <- ydim[min(length(ydim), i)]
    nClus_i <- nClus[min(length(nClus), i)]

    message("Batch correcting using a SOM grid of dimensions ", xdim_i,"x", ydim_i)

    if (is(label, "NULL")) {
      # Normalize data
      normalized <- .normalize(
        object,
        assay = assay,
        markers = markers,
        norm_method = norm_method,
        ties.method = ties.method,
        mc.cores = mc.cores,
        pb = pb)

      # Create SOM on normalized data
      label_i <- .create_som(
        normalized,
        rlen = rlen,
        mode = mode,
        seed = seed,
        xdim = xdim_i,
        ydim = ydim_i,
        nClus = nClus_i,
        markers = markers,
        resolution = resolution,
        cluster_method = cluster_method)
      rm(normalized)
    } else {
      if (length(label) == 1 & is(label, "character")) {
        if (label %in% colnames(metadata)) label_i <- metadata[[label]]
        else label_i <- label
      } else {
        label_i <- label
      }
    }

    metadata$label <- label_i



    # Run batch correction
    mat <- .correct_data(
      object2mat(object, assay = assay),
      metadata = metadata,
      label = label_i,
      covar = covar,
      anchor = anchor,
      markers = markers,
      parametric = parametric,
      method = method,
      ref.batch = ref.batch,
      mc.cores = mc.cores,
      pb = pb
    )
    # message(head(colnames(mat)))
    object <- add_mat(object, mat, assay = "cyCombine")
    rm(mat)
    assay <- "cyCombine"
  }

  return(object)
}
