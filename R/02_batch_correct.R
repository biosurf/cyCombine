#### Batch correction ####


#' Batch-wise scaling of data
#'
#' This function scales the data in a batch-wise manner.
#'   The purpose is to minimize the impact of batch correction when clustering the data prior to batch correction.
#'
#' @param df Dataframe with expression values
#' @param markers Markers to scale. If NULL, markers will be found using the \code{\link{get_markers}} function.
#' @family batch
#' @examples
#' df_scaled <- preprocessed %>%
#'   scale_expr()
#' @export
scale_expr <- function(df, markers = NULL){
  message("Scaling expression data..")
  if(is.null(markers)){
    # Get markers
    markers <- df %>%
      cyCombine::get_markers()
  }

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
#' @examples
#' som_ <- preprocessed %>%
#'   create_som()
#' @export
create_som <- function(df,
                       markers = NULL,
                       seed = 473,
                       xdim = 8,
                       ydim = 8){
  if(is.null(markers)){
    # Get markers
    markers <- df %>%
      cyCombine::get_markers()
  }
  # 10x10 SOM grid on overlapping markers, extract clustering per cell
  message("Creating SOM grid.. (Depending on the size of the data set, this may take a while)")
  set.seed(seed)
  som <- df %>%
    dplyr::select(markers) %>%
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
    message(paste("SOM class:", s))


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
#' @param label The cluster or cell type label. Either as a column name or vector.
#' @param covar The covariate ComBat should use. Can be a vector or a column name in the input datafrome.
#'   If NULL, no covar will be used
#' @param parametric Default: TRUE. If TRUE, the parametric version of ComBat is used. If FALSE, the non-parametric version is used.
#' @inheritParams scale_expr
#' @family batch
#' @examples
#' corrected <- preprocesed %>%
#'   correct_data(som_classes = som_$unit.classif, covar = "condition")
#' @export
correct_data <- function(df,
                         label,
                         covar = NULL,
                         markers = NULL,
                         parametric = TRUE){
  message("Batch correcting data..")
  if (is.null(markers)){
    # Get markers
    markers <- df %>%
      cyCombine::get_markers()
  }

  # Add label to df
  if(length(label) == 1){
    check_colname(colnames(df), label)
  }else{
    df$label <- label
    label <- "label"
  }

  # Add covar to df, if given
  if(is.null(covar)){
    # No covar is given
    num_covar <- 1
  }else if(length(covar) == 1){
    check_colname(colnames(df), covar)
  } else{
    # Covar was given as a vector
    df$covar <- as.factor(covar)
    covar <- "covar"
  }

  corrected_data <- df %>%
    dplyr::group_by(.data[[label]]) %>%
    # Correct (modify) each label group with ComBat
    dplyr::group_modify(.keep = TRUE, function(df, ...){
      # Detect if only one batch is present in the node
      num_batches <- df$batch %>%
        unique() %>%
        length()
      if(num_batches == 1){
        lab <- df[[label]][1]
        batch <- df$batch[1]
        message(paste("Label group", lab, "only contains cells from batch", batch))
        df <- df %>% select(-c(label))
        return(df)
      }
      message(paste("Correcting Label group", df[[label]][1]))
      # Calculate number of covars in the node
      if(!is.null(covar)){
        num_covar <- df[[covar]] %>%
          unique() %>%
          length()

        covar_counts <- df %>%
          count(.data$batch, .data[[covar]]) %>%
          pull(.data$n)

        if(sum(covar_counts) < max(covar_counts) + num_covar*5){
          num_covar <- 1
        }
      }

      # Compute ComBat correction
      ComBat_output <- df %>%
        dplyr::select(all_of(markers)) %>%
        t() %>%
        # The as.character is to remove factor levels not present in the SOM node
        purrr::when(num_covar > 1 ~ sva::ComBat(.,
                                                batch = as.character(df$batch),
                                                mod = stats::model.matrix(~as.character(df[[covar]])),
                                                par.prior = parametric),
                    ~ sva::ComBat(.,
                                  batch = as.character(df$batch),
                                  par.prior = parametric)
                    ) %>%
        t() %>%
        tibble::as_tibble() %>%
        dplyr::mutate(batch = df$batch,
                      sample = df$sample,
                      id = df$id) %>%
        # Only add covar column, if it is not null
        purrr::when(!is.null(covar) ~ dplyr::mutate(., covar = df[[covar]]),
                    ~ .)
      return(ComBat_output)
    }) %>%
    dplyr::ungroup() %>%
    # Reduce all negative values to zero
    dplyr::mutate_at(dplyr::vars(all_of(markers)),
                     function(x) {
                       x[x < 0] <- 0
                       return(x)
                       }) %>%
    dplyr::arrange(id) %>%
    select(id, everything())
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
batch_correct <- function(df,
                          label = NULL,
                          xdim = 8,
                          ydim = 8,
                          parametric = TRUE,
                          seed = 473,
                          covar = NULL,
                          markers = NULL,
                          batch_col = "batch"){

  check_colname(colnames(df), batch_col)

  # Create SOM on scaled data
  if(is.null(label)){
    som_ <- df %>%
      scale_expr(markers = markers) %>%
      create_som(markers = markers,
                 seed = seed,
                 xdim = xdim,
                 ydim = ydim)
    label <- som_$unit.classif
  }


  # Run batch correction
  corrected <- df %>%
    correct_data(label = label,
                 covar = covar,
                 markers = markers,
                 parametric = parametric)
  message("Done!")
  return(corrected)
}





# ALT ----
#' Correct data using ComBat
#'
#' This function computes the batch correction on the preprocessed data using the ComBat algorithm.
#'
#' @importFrom sva ComBat
#' @importFrom stats model.matrix
#' @param som_classes The classes as returned by the \code{\link{create_som}} function
#' @param covar The covariate ComBat should use. Can be a vector or a column name in the input datafrome.
#'   If NULL, no covar will be used
#' @param parametric Default: TRUE. If TRUE, the parametric version of ComBat is used. If FALSE, the non-parametric version is used.
#' @inheritParams scale_expr
#' @family batch
#' @examples
#' corrected <- preprocesed %>%
#'   correct_data(som_classes = som_$unit.classif, covar = "condition")
correct_data_alt <- function(df,
                         som_classes,
                         covar = NULL,
                         markers = NULL,
                         parametric = TRUE){
  message("Batch correcting data..")
  if (is.null(markers)){
    # Get markers
    markers <- df %>%
      cyCombine::get_markers()
  }


  if(is.null(covar)){
    # No covar is given
    df$som <- som_classes
    num_covar <- 1
  }else if(class(covar) == "character" & length(covar) == 1){
    # check_colname(colnames(df), covar)
    # Covar is in the data.frame
    df <- df %>%
      dplyr::mutate(som = som_classes,
                    covar = df[[covar]] %>%
                      as.factor())
  } else{
    # Covar was given as a vector
    df <- df %>%
      dplyr::mutate(som = som_classes,
                    covar = as.factor(covar))
  }

  # Detect if only one batch is present in the node
  num_batches <- df$batch %>%
    unique() %>%
    length()
  if(num_batches == 1){
    stop("Only one batch in data.")
  }
  if(!is.null(covar)){
    num_covar <- df$covar %>%
      unique() %>%
      length()
  }


    # Compute ComBat correction
    ComBat_output <- df %>%
      dplyr::select(all_of(markers)) %>%
      t() %>%
      # The as.character is to remove factor levels not present in the SOM node
      purrr::when(num_covar > 1 ~ sva::ComBat(.,
                                              batch = as.character(df$batch),
                                              mod = stats::model.matrix(~df$covar + df$som),
                                              par.prior = parametric),
                  ~ sva::ComBat(.,
                                batch = as.character(df$batch),
                                mod = stats::model.matrix(~df$som),
                                par.prior = parametric)
      ) %>%
      t() %>%
      tibble::as_tibble() %>%
      dplyr::mutate(batch = df$batch,
                    sample = df$sample,
                    id = df$id) %>%
      # Only add covar column, if it is not null
      purrr::when(!is.null(covar) ~ dplyr::mutate(., covar = df$covar),
                  ~ .) %>%
      # Reduce all negative values to zero
      dplyr::mutate_at(dplyr::vars(all_of(markers)),
                       function(x) {
                         x[x < 0] <- 0
                         return(x)
                         }) %>%
      dplyr::arrange(id)
    return(ComBat_output)
}
