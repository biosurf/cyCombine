# Normalization ----

#' Batch-wise Normalization of Data Using Seurat
#'
#' This function normalizes the data in a batch-wise manner using Seurat.
#' The purpose is to minimize the impact of batch effects when clustering the data prior to batch correction.
#' Four normalization methods are implemented: Z-score, Rank, CLR, and Quantile normalization.
#' Z-score is recommended for batches from a single study/experiment.
#' Rank is recommended for data from different studies/experiments.
#' CLR is recommended for ADT data.
#' Quantile is not recommended.
#'
#' @inheritParams normalize
#' @param object A Seurat object
#' @param norm_method Normalization method: "scale" (Z-score), "CLR", "lognorm", "rank", or "qnorm" (Quantile normalization). Defaults to "scale".
#' @param layer Layer to use from the Seurat object
#'
#' @return A Seurat object with normalized data.
#' @export
normalize_seurat <- function(
    object,
    markers = NULL,
    layer = "counts",
    norm_method = c("scale", "rank", "CLR", "lognorm", "qnorm", "none"),
    ties.method = c("average", "first", "last", "random", "max", "min"),
    mc.cores = 1,
    pb = FALSE) {
  cyCombine:::check_package("Seurat")
  # Remove case-sensitivity
  norm_method <- match.arg(norm_method)
  ties.method <- match.arg(ties.method)
  APPLY <- set_apply(mc.cores, pb)

  if (is.null(markers)) {
    # Get markers
    markers <- rownames(object)
  }
  object <- object[markers,]

  # Messaging
  object <- switch(norm_method,
         "rank" = {
           message("Ranking expression data..")
           rank_norm_seurat(object, markers = markers, ties.method = ties.method, APPLY = APPLY, layer = layer)
           },
         "scale" = {
           message("Scaling expression data..")
           Seurat::ScaleData(object, features = markers, split.by = "batch")
           },
         "CLR" = {
           message("CLR normalizing expression data..")
           object <- Seurat::NormalizeData(object, normalization.method = "CLR", margin = 2)
           Seurat::ScaleData(object, features = markers, split.by = "batch")
         },
         "lognorm" = {
           message("Log-normalizing expression data..")
           object <- Seurat::NormalizeData(object)
           Seurat::ScaleData(object, features = markers, split.by = "batch")
         },
         "qnorm" = {
           message("Quantile normalizing expression data..")
           quantile_norm_seurat(object, markers = markers)
         },
         "none" = {
           message("Skipping the normalization..")
           object
         }
  )

  return(object)
}

# Batch-wise quantile normalization per marker using Seurat

rank_norm_seurat <- function(object, markers = NULL, ties.method, APPLY, layer) {
  # Rank at marker positions individually for every batch
  batches <- unique(object$batch)
  ranked_data <- APPLY(
    setNames(batches, batches), function(b) {
      batch_cells <- SeuratObject::WhichCells(object, expression = batch == b)
      data <- as.matrix(SeuratObject::LayerData(object, layer)[markers, batch_cells])
      # rowzeros <- rowSums(data) == 0
      # if (any(rowzeros)) {
      #   warning("A marker is 0 for an entire batch. This marker is removed.")
      #   data <- data[!rowzeros, ]
      # }
      data <- apply(data, 2, rank, ties.method = ties.method) / ncol(data)
    })
  ranked_data <- do.call(cbind, ranked_data)

  ranked_data <- ranked_data[, match(colnames(object), colnames(ranked_data))]

  SeuratObject::LayerData(object, "scale.data") <- ranked_data
  return(object)
}

quantile_norm_seurat <- function(object, markers = NULL, mc.cores = 1, pb = FALSE) {
  APPLY <- set_apply(mc.cores, pb)
  message("Quantile normalizing expression data..")
  if (is.null(markers)) {
    markers <- rownames(object)
  }

  # Determine goal distributions for each marker by getting quantiles across all batches
  refq <- APPLY(setNames(markers, markers), function(m) {
    quantile(SeuratObject::LayerData(object, "data")[m, ], probs = seq(0, 1, length.out = 5), names = FALSE)
  })

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
#' @inheritParams create_som
#' @inheritParams normalize_seurat
#' @param resolution Resolution parameter for lauvain/leiden
#'
#' @return A vector of clustering labels
#' @export
create_som_seurat <- function(
    object,
    layer = SeuratObject::Layers(object)[length(SeuratObject::Layers(object))],
    markers = NULL,
    seed = 473,
    rlen = 10,
    mode = c("online", "batch", "pbatch"),
    cluster_method = c("kohonen", "flowsom", "lauvain", "leiden"),
    resolution = 0.8,
    distf = c("euclidean", "sumofsquares", "cosine", "manhattan", "chebyshev"),
    xdim = 8,
    ydim = 8,
    nClus = NULL) {

  mode <- match.arg(mode)
  distf <- match.arg(distf)
  cluster_method <- match.arg(cluster_method)
  # Default to all genes if markers are not specified
  if (is.null(markers)) {
    markers <- rownames(object)
  }

  if (cluster_method == "kohonen") {
    if (!distf %in% c("sumofsquares", "euclidean", "manhattan")) {
      warning("Distance function '", dist, "' not supported by Kohonen. Setting to 'sumofsquares'")
      distf <- "sumofsquares"
    }
    # Extract data for the markers
    data <- SeuratObject::LayerData(object, layer)[markers, ]

    # SOM grid on overlapping markers, extract clustering per cell
    message("Creating SOM grid..")
    set.seed(seed)
    som_grid <- kohonen::somgrid(xdim = xdim, ydim = ydim, topo = "rectangular")
    som_model <- kohonen::som(
      t(data), grid = som_grid, rlen = rlen, dist.fcts = distf, mode = mode)

    # Add labels to metadata
    object <- SeuratObject::AddMetaData(object, metadata = som_model$unit.classif, col.name = "Labels")
  } else if (cluster_method == "flowsom") {
    if (distf == "sumofsquares") {
      warning("Distance function '", dist, "' not supported by FlowSOM. Setting to 'euclidean'")
      distf <- "euclidean"
    }
    cyCombine:::check_package("FlowSOM", "Bioc")
    distf <- switch(distf, "manhattan" = 1, "euclidean" = 2, "chebyshev" = 3, "cosine" = 4)
    # Extract data for the markers
    data <- SeuratObject::LayerData(object, layer)[markers, ]
    message("Creating SOM grid..")
    set.seed(seed)
    fsom <- t(data) |>
      FlowSOM::ReadInput() |>
      FlowSOM::BuildSOM(
      xdim = xdim,
      ydim = ydim,
      rlen = rlen,
      distf = distf
    )
    if (is(nClus, "NULL")){
      # Add labels to metadata
      object <- SeuratObject::AddMetaData(object, metadata = FlowSOM::GetClusters(fsom), col.name = "Labels")
    } else{
      fsom <- FlowSOM::BuildMST(fsom, tSNE = FALSE)

        meta <- FlowSOM::metaClustering_consensus(
          fsom$map$codes,
          k = nClus,
          seed = seed)
        # Add labels to metadata
        object <- SeuratObject::AddMetaData(object, metadata = meta[fsom$map$mapping[, 1]], col.name = "Labels")
    }
  } else if (cluster_method %in% c("leiden", "lauvain")) {
    if (!"pca" %in% Seurat::Reductions(object)) object <- Seurat::RunPCA(object)
    object <- Seurat::FindNeighbors(
      object,
      reduction = "pca",
      compute.SNN = TRUE)

    object <- Seurat::FindClusters(
      object,
      random.seed = seed,
      algorithm = switch(cluster_method, "leiden" = 4, "lauvain" = 1),
      resolution = resolution,
      random.seed = seed)
    object <- SeuratObject::AddMetaData(object, metadata = object$seurat_clusters, col.name = "Labels")
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
#' @inheritParams correct_data
#' @inheritParams create_som_seurat
#' @param return_seurat Logical. Whether the matrix or a Seurat object should be returned.
#'
#' @return A Seurat object with corrected data.
#' @export
correct_data_seurat <- function(
    object,
    markers = NULL,
    method = c("ComBat", "ComBat_seq"),
    covar = NULL,
    anchor = NULL,
    ref.batch = NULL,
    parametric = TRUE,
    return_seurat = TRUE,
    mc.cores = 1,
    pb = FALSE) {
  cyCombine:::check_package("Seurat")

  method <- match.arg(method)
  APPLY <- set_apply(mc.cores, pb)
  message("Batch correcting..")

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


  labels <- unique(object$Labels)
  corrected_data <- APPLY(
    setNames(labels, labels),
    function(lab) {
      label_cells <- SeuratObject::WhichCells(object, expression = `Labels` == lab)
      correct_label_seurat(
        object[markers, label_cells],
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

  # Ensure data is in the same order as the original Seurat object
  flawed_cells <- !colnames(corrected_data) %in% colnames(object)
  if (sum(flawed_cells) > 0) {
    warning(sum(flawed_cells), " cell(s) were excluded in correction and are removed.")
    object <- object[, colnames(corrected_data)]
  }

  # Update the Seurat object with corrected data

  object[["cyCombine"]] <- SeuratObject::CreateAssayObject(data = corrected_data, key = "corrected_")
  SeuratObject::DefaultAssay(object) <- "cyCombine"

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

  layer <- ifelse(method == "ComBat", "data", "counts")
  data <- .combat(
    as.matrix(SeuratObject::LayerData(object_lab, layer = layer)),
    batch = object_lab$batch,
    mod_matrix,
    parametric,
    ref.batch,
    method)

  return(data)
}

# Function to check confounding
# check_confound <- function(batch, covariate) {
#   covariate_levels <- unique(covariate)
#   for (level in covariate_levels) {
#     if (length(unique(batch[covariate == level])) == 1) {
#       return(TRUE)
#     }
#   }
#   return(FALSE)
# }

# Wrapper ----

#' Run batch correction on data
#'
#' This is a wrapper function for the cyCombine batch correction workflow.
#'  To run the workflow manually, type "batch_correct" to see the source code of this wrapper and follow along or read the vignettes on the GitHub page \url{https://github.com/biosurf/cyCombine}.
#'
#' @inheritParams batch_correct
#' @inheritParams create_som_seurat
#' @inheritParams correct_data_seurat
#' @inheritParams normalize_seurat
#' @family batch
#' @importFrom sva ComBat ComBat_seq
#' @importFrom methods slot slot<-
#' @import stats
#' @examples
#' \dontrun{
#' seu <- batch_correct_seurat(
#'  seu,
#'  markers = markers,
#'  covar = "condition"
#'  )
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
    cluster_method = c("kohonen", "flowsom", "lauvain", "leiden"),
    resolution = 0.8,
    ref.batch = NULL,
    seed = 473,
    label = NULL,
    covar = NULL,
    anchor = NULL,
    markers = NULL,
    layer = "counts",
    norm_method = c("scale", "rank", "CLR", "lognorm", "qnorm", "none"),
    ties.method = "average",
    mc.cores = 1,
    pb = FALSE) {

  APPLY <- set_apply(mc.cores, pb)
  cyCombine:::check_package("Seurat")

  stopifnot(
    "No 'batch' column in data." = "batch" %in% names(object[[]]))

  if (is.null(markers)) {
    markers <- rownames(object)
  }

  if (SeuratObject::DefaultAssay(object) == "SCT"){
    ## Workaround for correctly subsetting SCT assays.
    ## Published on github by longmanz
    ## https://github.com/satijalab/seurat-object/issues/208

    tmp_SCT_features_attributes <- slot(object = object[['SCT']], name = "SCTModel.list")[[1]]@feature.attributes
    tmp_SCT_features_attributes <- tmp_SCT_features_attributes[markers, ]
    slot(object = object[['SCT']], name = "SCTModel.list")[[1]]@feature.attributes <- tmp_SCT_features_attributes
    object <- object[markers, ] #subset(object, features = markers)
  } else {
    object <- object[markers, ]
  }

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
        mc.cores = mc.cores,
        pb = pb)
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
    object <- correct_data_seurat(
      object,
      covar = covar,
      anchor = anchor,
      markers = markers,
      parametric = parametric,
      method = method,
      ref.batch = ref.batch,
      mc.cores = mc.cores,
      pb = pb,
      return_seurat = TRUE
    )
    layer <- "data"
  }

  message("Done!")
  return(object)
}



