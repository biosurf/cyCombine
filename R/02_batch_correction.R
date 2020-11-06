#### Batch correction ####


#' Batch-wise scaling of data
#'
#' This function scales the data in a batch-wise manner.
#'   The purpose is to minimize the impact of batch correction when clustering the data prior to batch correction.
#'
#' @param df Dataframe with expression values
#' @family batch
#' @export
scale_expr <- function(df){
  cat("Scaling expression data\n")
  # Get markers
  markers <- df %>%
    cyCombine::get_markers()
  # Scale at marker positions
  scaled_expr <- df %>%
    dplyr::group_by(.data$batch) %>%
    dplyr::mutate_at(dplyr::vars(markers), .funs = scale) %>%
    dplyr::ungroup()
  return(scaled_expr)
}


#' Create Self-Organizing Map
#'
#' The function uses the kohonen package to create a Self-Organizing Map.
#'   It is used to segregate the cells for the batch correction to make the correction less affected
#'  by samples with high abundances of a particular cell type.
#'
#' @importFrom kohonen som somgrid
#' @inheritParams scale_expr
#' @param seed The seed to use when creating the SOM
#' @param xdim The x-dimension size of the SOM
#' @param ydim The y-dimension size of the SOM
#' @family batch
#' @export
create_som <- function(df,
                       seed = 473,
                       xdim = 10,
                       ydim = 10){
  # 10x10 SOM grid on overlapping markers, extract clustering per cell
  cat("Creating SOM grid\n")
  set.seed(seed)
  som <- df %>%
    dplyr::select_if(colnames(.) %!in% non_markers) %>%
    as.matrix() %>%
    kohonen::som(grid = kohonen::somgrid(xdim = xdim,
                                         ydim = ydim),
                 dist.fcts = "euclidean")
  return(som)
}




#' Correct data using ComBat
#'
#' Deprecated
#' @importFrom sva ComBat
correct_data_prev <- function(df,
                              som_classes){
  # Create copy dataset
  corrected_data <- df %>%
    dplyr::mutate(covar = "")


  for (s in sort(unique(som_classes))) {
    # Extract original (non-scaled+ranked) data for cluster
    data_subset <- df[which(som_classes==s), ]
    cat("som class:", s, "\n", sep = " ")


    # ComBat batch correction using disease status as covariate
    covar <- rep('CLL', nrow(data_subset))
    covar[grep('^HD', data_subset$sample)] <- 'HD'
    magic_output <- data_subset %>%
      dplyr::select_if(colnames(.) %!in% non_markers) %>%
      t() %>%
      sva::ComBat(batch = data_subset$batch, mod = stats::model.matrix(~covar)) %>%
      t() %>%
      tibble::as_tibble() %>%
      dplyr::mutate(batch = data_subset$batch,
             sample = data_subset$sample,
             covar = covar,
             id = data_subset$id)

    # Fill copy dataset with corrected data
    corrected_data[which(som_classes==s), ] <- magic_output
  }
  # Fix values below zero
  corrected_data[corrected_data < 0] <- 0
  corrected_data <- corrected_data %>%
    arrange(id)

  return(corrected_data)
}

#' Correct data using ComBat
#'
#' This function computes the batch correction on the preprocessed data using the ComBat algorithm.
#'
#' @importFrom sva ComBat
#' @importFrom stats model.matrix
#' @param som_classes The classes as returned by the \code{\link{create_som}} function
#' @param covar The covariate ComBat uses. If NULL, MAKE MORE FLEXIBLE!
#' @inheritParams scale_expr
#' @family batch
#' @export
correct_data <- function(df,
                         som_classes,
                         covar = NULL,
                         markers = NULL){
  if (is.null(markers)){
    # Get markers
    markers <- df %>%
      cyCombine::get_markers()
  }


  if(is.null(covar)){
    df <- df %>%
      dplyr::mutate(som = som_classes,
                    # Determine covariate
                    covar = case_when(stringr::str_starts(string = sample,
                                                          pattern = "HD") ~ "HD",
                                      TRUE ~ "CLL") %>%
                      as.factor())
  }else{
    df <- df %>%
      dplyr::mutate(som = som_classes,
                    covar = as.factor(covar))
  }

  corrected_data <- df %>%
    dplyr::group_by(som) %>%
    # Run ComBat on each SOM class
    dplyr::group_modify(function(df, ...){

      ComBat_output <- df %>%
        dplyr::select_if(colnames(.) %!in% non_markers) %>%
        t() %>%
        # The as.character is to remove factor levels not present in the SOM node
        sva::ComBat(batch = as.character(df$batch), mod = stats::model.matrix(~df$covar)) %>%
        t() %>%
        tibble::as_tibble() %>%
        dplyr::mutate(batch = df$batch,
                      sample = df$sample,
                      covar = df$covar,
                      id = df$id)
      return(ComBat_output)
    }) %>%
    dplyr::ungroup() %>%
    # dplyr::select(-som) %>%
    # Reduce all negative values to zero
    dplyr::mutate_at(dplyr::vars(markers),
                     function(x) {
                       x[x < 0] <- 0
                       return(x)
                       }) %>%
    dplyr::arrange(id)
  return(corrected_data)
}



#' Run batch correction on preprocessed data
#'
#' This is a wrapper function for the cyCombine batch correction workflow.
#'   To run the workflow manually, type "batch_correct" to see the source code of this wrapper and follow along.
#'
#' @inheritParams create_som
#' @inheritParams correct_data
#' @param preprocessed The preprocessed dataframe to run batch correction on
#' @family batch
#' @export
batch_correct <- function(preprocessed,
                          xdim = 10,
                          ydim = 10,
                          seed = 473,
                          covar = NULL,
                          markers = NULL){

  # Create SOM on scaled data
  som <- preprocessed %>%
    scale_expr() %>%
    create_som(seed = seed,
               xdim = xdim,
               ydim = ydim)

  # Run batch correction
  cat("Batch correcting data\n")
  corrected <- preprocessed %>%
    correct_data(som_classes = som$unit.classif,
                 covar = covar,
                 markers = markers)
  cat("Done!\n")
  return(corrected)
}





