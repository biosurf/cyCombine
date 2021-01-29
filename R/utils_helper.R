#' Inverted versions of in
#'
#' @noRd
#' @examples
#' 1 %!in% 1:10
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
check_colname <- function(df_colnames, col_name, location = "metadata"){
  if(!is.null(col_name)){
    if(col_name %!in% df_colnames){
      stop("Column \"", col_name, "\" was not found in the ", location)
    }}
}


#' Check if directory exists, if not, make it
#' @export
check_make_dir <- function(dir.path) {
  if (!dir.exists(dir.path)) {dir.create(dir.path)}
}
