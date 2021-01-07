#### Batch correction ####


#' Batch-wise scaling of data
#'
#' This function scales the data in a batch-wise manner.
#'   The purpose is to minimize the impact of batch correction when clustering the data prior to batch correction.
#'
#' @param df Dataframe with expression values
#' @param markers Markers to normalize. If NULL, markers will be found using the \code{\link{get_markers}} function.
#' @param norm_method Normalization method. Should be either 'rank', 'scale' or 'qnorm'. Default: 'rank'
#' @param ties.method The method to handle ties, when using rank. Default: "min". See ?rank for other options.
#' @family batch
#' @examples
#' df_normed <- df %>%
#'   normalize()
#' @export
normalize <- function(df, markers = NULL, norm_method = "rank", ties.method = "min"){

  norm_method <- norm_method %>% stringr::str_to_lower()

  if(norm_method == "rank") message("Ranking expression data..")
  else if(norm_method == "scale") message("Scaling expression data..")
  else if(norm_method == "qnorm") {
    # message("Quantile normalizing expression data..")
    df_normed <- quantile_norm(df, markers = markers)
    return(df_normed)
    } else stop("Please use either 'scale', 'rank', or 'qnorm' as normalization method." )

  if(is.null(markers)){
    # Get markers
    markers <- df %>%
      cyCombine::get_markers()
  }

  # Scale at marker positions
  df_normed <- df %>%
    dplyr::group_by(.data$batch) %>%
    purrr::when(norm_method == "rank"  ~ dplyr::mutate(., dplyr::across(dplyr::all_of(markers),
                                                            .fns = ~ {rank(.x, ties.method = ties.method) / length(.x)})),
                norm_method == "scale" ~ dplyr::mutate(., dplyr::across(dplyr::all_of(markers),
                                                            .fns = scale))
    ) %>%
    dplyr::ungroup()
  return(df_normed)
}


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


#' Batch-wise quantile normalization per marker
#'
#' This function quantile normalizes the data in a batch-wise manner.
#'   The purpose is to minimize the impact of batch correction when clustering the data prior to batch correction.
#'
#' @param df Dataframe with expression values
#' @param markers Markers to correct. If NULL, markers will be found using the \code{\link{get_markers}} function.
#' @family batch
#' @examples
#' df_qnorm <- preprocessed %>%
#'   quantile_norm()
#' @export
quantile_norm <- function(df, markers = NULL){
  message("Quantile normalizing expression data..")
  if(is.null(markers)){
    # Get markers
    markers <- df %>%
      cyCombine::get_markers()
  }

  # Determine goal distributions for each marker by getting quantiles across all batches
  refq <- list()
  for (m in markers) {
    # Determine the quantiles
    refq[[m]] <- quantile(unlist(df[,m]), probs=seq(0,1,length.out=5), names = F)
  }

  qnormed_expr <- df
  for (b in unique(df$batch)) {
    for (m in markers) {
      qx <- quantile(unlist(df[df$batch == b, m]), probs=seq(0,1,length.out=5), names = F)
      spf <- splinefun(x=qx, y=refq[[m]], method="monoH.FC", ties=min)

      # Apply the spline function to adjust quantiles
      qnormed_expr[qnormed_expr$batch == b, m] <- spf(unlist(df[df$batch==b, m]))

    }
  }

  return(qnormed_expr)
}

#' Batch-wise ranking per marker
#'
#' This function ranks and scales the data in a batch-wise manner.
#'   The purpose is to minimize the impact of batch correction when clustering the data prior to batch correction.
#'
#' @param df Dataframe with expression values
#' @param markers Markers to correct. If NULL, markers will be found using the \code{\link{get_markers}} function.
#' @family batch
#' @examples
#' df_ranked <- preprocessed %>%
#'   ranking()
#' @export
ranking <- function(df, markers = NULL, ties.method = "min") {

  message("Ranking expression data..")
  if(is.null(markers)){
    # Get markers
    markers <- df %>%
      cyCombine::get_markers()
  }


  # # Ranking of each marker in each cell - reversing order
  # ranking_per_row <- t(apply(df[,markers], 1, rank, ties.method = 'average')) # Probably wont work that well due to many 0's and consequently assigning 0.1 a large rank...

  # Ranking per batch/sample...? for each marker
  adj_df <- df
  for (m in markers) {
    for (b in unique(df$batch)) {
      #adj_df[df$batch == b, m] <- scale(rank(df[df$batch == b, m], ties.method = 'average')) # We can discuss the ties.method, we have to normalize within the batch to avoid large batches getting larger max rank
      adj_df[df$batch == b, m] <- rank(df[df$batch == b, m], ties.method = ties.method) / nrow(df[df$batch == b,]) # We can discuss the ties.method, we have to normalize within the batch to avoid large batches getting larger max rank

      # scaled <- apply(df[df$batch == b, markers], 2, scale)
      # adj_df[df$batch == b, markers] <- t(apply(scaled, 1, rank))

    }
  }

  return(adj_df)
}



#' Create Self-Organizing Map
#'
#' The function uses the kohonen package to create a Self-Organizing Map.
#'   It is used to segregate the cells for the batch correction to make the correction less affected
#'  by samples with high abundances of a particular cell type.
#'
#' @importFrom kohonen som somgrid
#' @importFrom stats predict
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

  # Predict runtime
  pred <- stats::predict(model, tibble::tibble("Size" = nrow(df)))

  # 10x10 SOM grid on overlapping markers, extract clustering per cell
  message("Creating SOM grid.. (This is estimated to take ", round(pred, 2), " minutes)")
  set.seed(seed)
  som_grid <- df %>%
    dplyr::select(markers) %>%
    as.matrix() %>%
    kohonen::som(grid = kohonen::somgrid(xdim = xdim,
                                         ydim = ydim),
                 dist.fcts = "euclidean")
  return(som_grid)
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
#'   correct_data(label = som_$unit.classif, covar = "condition")
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
                       x[x > 30] <- 30
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
#' @inheritParams normalize
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
                          norm_method = 'rank',
                          ties.method = "min"){
  # A batch column is required
  check_colname(colnames(df), "batch")

  # Create SOM on scaled data
  if(is.null(label)) {
    som_ <- df %>%
      normalize(markers = markers,
                norm_method = norm_method,
                ties.method = ties.method) %>%
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
