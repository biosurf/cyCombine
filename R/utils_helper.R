#' Inverted versions of in, is.null and is.na
#'
#' @noRd
#'
#' @examples
#' 1 %not_in% 1:10
#' not_null(NULL)
`%!in%` <- Negate(`%in%`)


#' Run PCA analysis
#' @importFrom stats prcomp
run_pca <- function(df, pcs = 20, non_markers = c("batch", "sample", "covar", "som", "label", "id")){
  df <- df %>%
    dplyr::select_if(names(.) %!in% non_markers)
  # Run PCA
  pca <- df %>%
    stats::prcomp(scale. = TRUE, center = TRUE)

  return(pca$x[, 1:pcs])
}

#' Cumpute flowsom clustering
#' @importFrom flowCore flowFrame colnames
#' @importFrom FlowSOM ReadInput BuildSOM BuildMST metaClustering_consensus
run_flowsom <- function(dataset, k = 7, seed = 473, non_markers = c("batch", "sample", "covar", "som", "label", "id")){
  # Create FlowFrame from data

  data_FlowSOM <- dataset %>%
    select_if(colnames(.) %!in% non_markers) %>%
    as.matrix() %>%
    flowCore::flowFrame()
  # Set seed for reproducibility
  set.seed(seed)

  # Run FlowSOM (initial steps prior to meta-clustering)
  fSOM <- data_FlowSOM %>%
    FlowSOM::ReadInput(transform = FALSE, scale = FALSE, silent = TRUE) %>%
    FlowSOM::BuildSOM(colsToUse = colnames(data_FlowSOM), silent = TRUE) %>%
    FlowSOM::BuildMST(silent = TRUE)
  # Extract cluster labels (pre meta-clustering) from output object
  labels_pre <- fSOM$map$mapping[, 1]

  # Run FlowSOM meta-clustering
  fSOM_Cluster <- FlowSOM::metaClustering_consensus(fSOM$map$codes, seed = seed, k = k)

  labels <- fSOM_Cluster[labels_pre]

  return(labels)
}