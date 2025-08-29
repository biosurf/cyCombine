#' Inverted versions of in
#'
#' @noRd
#' @examples
#' 1 %!in% 1:10
`%!in%` <- Negate(`%in%`)

#' Utility function for NULL coalescing
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x

#'
#' Define apply function to use
#' @noRd
set_apply <- function(mc.cores, pb = TRUE) {
  if(mc.cores == 1 & !pb) {
    APPLY <- lapply
  } else {
    if(pb) {
      cyCombine:::check_package("pbmcapply")
      APPLY <- pbmcapply::pbmclapply
    } else {
      cyCombine:::check_package("parallel")
      APPLY <- parallel::mclapply
    }
    formals(APPLY)$mc.cores <- mc.cores
  }
  return(APPLY)
}

#' Randomize a matrix of values
#'
#' Subtracts a random value between 0 and 1 from all values. Then convert negative values to zeros.
#' @param mat A matrix to randomize
randomize_matrix <- function(mat) {
  # Matrix dimensions
  ncols <- ncol(mat)
  nrows <- nrow(mat)

  # Ceiling values
  randomized <- ceiling(mat)
  # Subtract random number
  randomized <- randomized - matrix(stats::runif(n = ncols*nrows, min = 0, max = 0.9999),
                                    nrow = nrows,
                                    ncol = ncols)
  # Ceiling negative values
  randomized[randomized < 0] <- 0
  return(randomized)
}


#' Get markers from a dataframe or vector of markers
#'
#' This function uses the global variable "non_markers".
#'   If the output contains markers you did not expect, you can add to non_markers like this:
#'   \code{non_markers <- c(non_markers, "remove1", "remove2")} and rerun get_markers()
#' @param x dataframe or vector to get the markers from
#' @export
get_markers <- function(x) {
  # Use global non_markers if available
  if (!is.null(.GlobalEnv$non_markers)) non_markers <- .GlobalEnv$non_markers

  if (inherits(x, "data.frame")) x <- colnames(x)
  else  x <- rownames(x)
  marker_pos <- tolower(x) %!in% tolower(non_markers)
  markers <- x[marker_pos]
  markers <- markers[which(!is.na(markers))]
  return(markers)
}


#' Run PCA analysis
#' @noRd
run_pca <- function(df, pcs = 20) {
  missing_package("stats", "CRAN")
  pca <- df %>%
    dplyr::select_if(names(.) %!in% non_markers) %>%
    # Run PCA
    stats::prcomp(scale. = TRUE, center = TRUE)

  return(pca$x[, 1:pcs])
}


#' Check colname
#' @noRd
check_colname <- function(df_colnames, col_name, location = "metadata") {
  if (!is.null(col_name)) {
    if (col_name %!in% df_colnames) {
      stop("Column \"", col_name, "\" was not found in the ", location)
    }}
}



#' Check if directory exists, if not, make it
#' @noRd
check_make_dir <- function(dir.path) {
  if (!dir.exists(dir.path)) {dir.create(dir.path)}
}


#' Check if the batch is confounded with the provided model
#'
#' Code adapted from sva::ComBat
#'
#' @noRd
check_confound <- function(batch, mod = NULL) {


  ### Code adapted from sva::ComBat

  ## Create batch model
  batch <- as.factor(batch)
  batchmod <- stats::model.matrix(~ -1 + batch)


  ## A few other characteristics on the batches
  n.batch <- nlevels(batch)
  batches <- list()
  for (i in 1:n.batch) {
    batches[[i]] <- which(batch == levels(batch)[i])
  } # list of samples in each batch
  n.batches <- sapply(batches, length)

  n.array <- sum(n.batches)
  ## combine batch variable and covariates
  design <- cbind(batchmod, mod)

  ## check for intercept in covariates, and drop if present
  check <- apply(design, 2, function(x) all(x == 1))

  design <- as.matrix(design[, !check])

  ## Check if the design is confounded
  if (qr(design)$rank < ncol(design)) {
    ## if(ncol(design)<=(n.batch)){stop("Batch variables are redundant! Remove one or more of the batch variables so they are no longer confounded")}
    if (ncol(design) == (n.batch+1)) {
      # message("The covariate is confounded with batch!") # Remove the covariate")
      return(TRUE)
    }
    if (ncol(design) > (n.batch+1)) {
      if ((qr(design[, -c(1:n.batch)])$rank<ncol(design[, -c(1:n.batch)]))) {
        # message("The covariates are confounded!") # Please remove one or more of the covariates so the design is not confounded")
        return(TRUE)
      } else {
        # message("At least one covariate is confounded with batch!") # Please remove confounded covariates")
        return(TRUE)
      }
    }
  }
  return(FALSE)
}

#' Wrapper for missing packages
#'
#' @noRd
check_package <- function(package, repo = "CRAN", git_repo = "") {

  if (repo == "CRAN") {
    install_function <- "install.packages('"
  } else if (repo == "github") {
    install_function <- paste0("devtools::install_github('", git_repo, "/")
  } else if (repo == "Bioc") {
    install_function <- "BiocManager::install('"
  }

  if (!requireNamespace(package, quietly = TRUE)) {
    stop(
      paste0("Package ", package, " is not installed.\n",
             "Please run: ", install_function, package, "')"))
  }
  requireNamespace(package, quietly = TRUE)
}

missing_package <- function(...) {
  check_package(...)
}
