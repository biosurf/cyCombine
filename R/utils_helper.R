#' Inverted versions of in, is.null and is.na
#'
#' @noRd
#'
#' @examples
#' 1 %!in% 1:10
#' not_null(NULL)
`%!in%` <- Negate(`%in%`)


#' Wrapper for missing packages
#'
#' @noRd
missing_package <- function(package, repo = "CRAN", git_repo = ""){

  if (repo == "CRAN"){
    install_function <- "install.packages('"
  } else if (repo == "github") {
    install_function <- paste0("devtools::install_github('", git_rep, "/")
  } else if (repo == "Bioc"){
    install_function <- "BiocManager::install('"
  }

  if(package %!in% rownames(installed.packages())){
    stop(paste0("Package ", package," is not installed.\n",
         "Please run: ", install_function, package, "')"))
  }
  requireNamespace(package, quietly = TRUE)
}



#' Get markers from a dataframe
#'
#' This function uses the global variable "non_markers".
#'   If the output contains markers you did not expect, you can add to non_markers like this:
#'   \code{non_markers <- c(non_markers, "remove1", "remove2")} and rerun get_markers()
#' @param df dataframe to get the markers from
#' @importFrom stringr str_to_lower
#' @export
get_markers <- function(df){
  marker_pos <- stringr::str_to_lower(colnames(df)) %!in% non_markers
  markers <- colnames(df)[marker_pos]
  return(markers)
}

#' Check colname
#' @noRd
check_colname <- function(df_colnames, col_name, location = "metadata"){
  if(!is.null(col_name)){
    if(col_name %!in% df_colnames){
      stop("Column \"", col_name, "\" was not found in the ", location)
    }}
}

#' Run PCA analysis
#' @importFrom stats prcomp
#' @noRd
run_pca <- function(df, pcs = 20){
  missing_package("stats", "CRAN")
  pca <- df %>%
    dplyr::select_if(names(.) %!in% non_markers) %>%
    # Run PCA
    stats::prcomp(scale. = TRUE, center = TRUE)

  return(pca$x[, 1:pcs])
}

# @importFrom FlowSOM ReadInput BuildSOM BuildMST metaClustering_consensus


#' Compute flowsom clustering
#' @importFrom flowCore flowFrame colnames
#' @noRd
run_flowsom <- function(dataset, k = 7, seed = 473){

  # Check for package
  missing_package("FlowSOM", "Bioc")

  # Create FlowFrame from data

  data_FlowSOM <- dataset %>%
    dplyr::select_if(colnames(.) %!in% non_markers) %>%
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

#' Check if directory exists, if not, make it
#' @noRd
check_make_dir <- function(dir.path) {
  if (!dir.exists(dir.path)) {dir.create(dir.path)}
}

col_max <- function(df, m){
  for (ma in m){
    val <- cor %>%
      dplyr::select(dplyr::all_of(ma)) %>%
      max()
    message(ma, ": ", round(val, 2))
  }
}


col_min <- function(df, m){
  for (ma in m){
    val <- cor %>%
      dplyr::select(dplyr::all_of(ma)) %>%
      min()
    message(ma, ": ", round(val, 2))
  }
}


#' Check if the batch is confounded with the provided model
#'
#' @noRd
check_confound <- function(dat, batch, markers = NULL, mod = NULL) {

  if (is.null(markers)){
    # Get markers
    markers <- df %>%
      cyCombine::get_markers()
  }
  ### Code adapted from sva::ComBat

  ## coerce dat into a matrix
  dat <- as.matrix(dat[, markers])

  batchmod <- model.matrix(~-1+batch)


  ## A few other characteristics on the batches
  n.batch <- nlevels(batch)
  batches <- list()
  for (i in 1:n.batch) {
    batches[[i]] <- which(batch == levels(batch)[i])
  } # list of samples in each batch
  n.batches <- sapply(batches, length)

  n.array <- sum(n.batches)
  ## combine batch variable and covariates
  design <- cbind(batchmod,mod)

  ## check for intercept in covariates, and drop if present
  check <- apply(design, 2, function(x) all(x == 1))

  design <- as.matrix(design[,!check])

  ## Check if the design is confounded
  if(qr(design)$rank < ncol(design)) {
    ## if(ncol(design)<=(n.batch)){stop("Batch variables are redundant! Remove one or more of the batch variables so they are no longer confounded")}
    if(ncol(design)==(n.batch+1)) {
      stop("The covariate is confounded with batch! Remove the covariate")
    }
    if(ncol(design)>(n.batch+1)) {
      if((qr(design[,-c(1:n.batch)])$rank<ncol(design[,-c(1:n.batch)]))){
        stop('The covariates are confounded! Please remove one or more of the covariates so the design is not confounded')
      } else {
        stop("At least one covariate is confounded with batch! Please remove confounded covariates")
      }
    }
  }
}
