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
#' @export
get_markers <- function(df){
  marker_pos <- colnames(df) %!in% non_markers
  markers <- colnames(df)[marker_pos]
  return(markers)
}

#' Check colname
#' @noRd
check_colname <- function(df_colnames, col_name){
  if(!is.null(col_name)){
    if(col_name %!in% df_colnames){
      stop("Column \"", col_name, "\" was not found in the metadata.")
    }}
}

#' Run PCA analysis
#' @importFrom stats prcomp
run_pca <- function(df, pcs = 20){
  missing_package("stats", "CRAN")
  pca <- df %>%
    dplyr::select_if(names(.) %!in% non_markers) %>%
    # Run PCA
    stats::prcomp(scale. = TRUE, center = TRUE)

  return(pca$x[, 1:pcs])
}

# @importFrom FlowSOM ReadInput BuildSOM BuildMST metaClustering_consensus


#' Cumpute flowsom clustering
#' @importFrom flowCore flowFrame colnames
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
