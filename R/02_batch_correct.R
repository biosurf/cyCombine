# Normalization ----


#' Batch-wise normalization of data
#'
#' This function normalizes the data in a batch-wise manner.
#'   The purpose is to minimize the impact of batch correction when clustering the data prior to batch correction.
#'   Three normalisation methods are implemented: Z-score, Rank, and Quantile normalization.
#'   Z-score is recommended in cases where batches from a single study/experiment is merged.
#'   Rank is recommend in cases where data from different studies/experiments are merged.
#'   Quantile is not recommended.
#'
#' @param df tibble with expression values
#' @param markers Markers to normalize. If NULL, markers will be found using the \code{\link{get_markers}} function.
#' @param norm_method Normalization method. Should be either 'rank', 'scale' or 'qnorm'. Default: 'scale'
#' @param ties.method The method to handle ties, when using rank. Default: "average". See ?rank for other options.
#' @family batch
#' @examples
#' df_normed <- df %>%
#'   normalize()
#' @export
normalize <- function(df,
                      markers = NULL,
                      norm_method = "scale",
                      ties.method = "average"){

  # Remove case-sensitivity
  norm_method <- norm_method %>% stringr::str_to_lower()

  # Messaging
  if(norm_method == "rank") message("Ranking expression data..")
  else if(norm_method == "scale") message("Scaling expression data..")
  else if(norm_method == "qnorm") {
    # message("Quantile normalizing expression data..")
    # Run quantile normalization
    df_normed <- quantile_norm(df, markers = markers)
    return(df_normed)
    } else stop("Please use either 'scale', 'rank', or 'qnorm' as normalization method." )

  if(is.null(markers)){
    # Get markers
    markers <- df %>%
      cyCombine::get_markers()
  }

  # Scale or rank at marker positions individually for every batch
  df_normed <- df %>%
    dplyr::group_by(.data$batch) %>%
    purrr::when(norm_method == "rank"  ~ dplyr::mutate(., dplyr::across(dplyr::all_of(markers),
                                                                        .fns = ~{
                                                                          if(sum(.x) == 0) stop("A marker is 0 for an entire batch. Please remove this markers.")
                                                                          rank(.x, ties.method = ties.method) / length(.x)})),
                norm_method == "scale" ~ dplyr::mutate(., dplyr::across(dplyr::all_of(markers),
                                                                        .fns = ~{
                                                                          if(sum(.x) == 0) stop("A marker is 0 for an entire batch. Please remove this markers.")
                                                                          scale(.x)}))
    ) %>%
    dplyr::ungroup()
  return(df_normed)
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



# Clustering ----

#' Create Self-Organizing Map
#'
#' The function uses the FlowSOM or kohonen package to create a Self-Organizing Map.
#'  It is used to segregate the cells for the batch correction to make the correction less affected
#'  by samples with high abundances of a particular cell type.
#'  FlowSOM is recommended with Z-score normalization.
#'  kohonen is recommended with rank normalization
#'
#' @importFrom kohonen som somgrid
#' @importFrom stats predict
#' @importFrom FlowSOM SOM
#' @inheritParams normalize
#' @param seed The seed to use when creating the SOM.
#' @param xdim The x-dimension size of the SOM.
#' @param ydim The y-dimension size of the SOM.
#' @param rlen Number of times the data is presented to the SOM network
#' @family batch
#' @examples
#' labels <- uncorrected %>%
#'   create_som()
#' @export
#' @return A vector of clustering labels
create_som <- function(df,
                       markers = NULL,
                       seed = 473,
                       rlen = 10,
                       xdim = 8,
                       ydim = 8){
  if(is.null(markers)){
    # Get markers
    markers <- df %>%
      cyCombine::get_markers()
  }

  # SOM grid on overlapping markers, extract clustering per cell
  message("Creating SOM grid..")
  set.seed(seed)
  labels <- df %>%
    dplyr::select(markers) %>%
    as.matrix() %>%
    kohonen::som(grid = kohonen::somgrid(xdim = xdim, ydim = ydim),
                 rlen = rlen,
                 dist.fcts = "euclidean")

  labels <- labels$unit.classif

  return(labels)
}

#' Compute flowsom clustering (Deprecated)
#' @importFrom FlowSOM SOM
#' @export
create_fsom <- function(df,
                        markers = NULL,
                        seed = 473,
                        xdim = 8,
                        ydim = 8){

  warning("This function is deprecated. Please use 'create_som(som_type = 'fsom')' instead.")
  # Check for package
  missing_package("FlowSOM", "Bioc")

  # Get markers
  if(is.null(markers)){
    # Get markers
    markers <- df %>%
      cyCombine::get_markers()
  }


  # Create SOM grid
  set.seed(seed)
  fsom <- df %>%
    dplyr::select(dplyr::all_of(markers)) %>%
    as.matrix() %>%
    FlowSOM::SOM(xdim = xdim, ydim = ydim)

  label <- fsom$mapping[, 1]

  return(label)
}


# Batch correction ----

#' Correct data using ComBat
#'
#' Compute the batch correction on the data using the ComBat algorithm.
#'  Define a covariate, either as a character vector or name of tibble column.
#'  The covariate should preferable be the cell condition types, but can be any column that infers heteroneity in the data.
#'  The function assumes that the batch information is in the "batch" column and the data contains a "sample" column with sample information.
#'
#' @importFrom sva ComBat
#' @importFrom stats model.matrix
#' @param label The cluster or cell type label. Either as a column name or vector.
#' @param covar The covariate ComBat should use. Can be a vector or a column name in the input tibble.
#'   If NULL, no covar will be used
#' @param parametric Default: TRUE. If TRUE, the parametric version of ComBat is used. If FALSE, the non-parametric version is used.
#' @inheritParams normalize
#' @family batch
#' @examples
#' corrected <- uncorrected %>%
#'   correct_data(label = labels, covar = "condition")
#' @export
correct_data <- function(df,
                         label,
                         covar = NULL,
                         markers = NULL,
                         parametric = TRUE){
  message("Batch correcting data..")
  # Check for batch column
  check_colname(colnames(df), "batch", "df")
  check_colname(colnames(df), "sample", "df")
  if (is.null(markers)){
    # Get markers
    markers <- df %>%
      cyCombine::get_markers()
  }

  # Add ID column to retain data order
  if("id" %!in% colnames(df)) df$id <- 1:nrow(df)

  # Add label to df
  if(length(label) == 1){
    check_colname(colnames(df), label, "df")
  }else{
    df$label <- label
    label <- "label"
  }

  # Add covar to df, if given
  if(is.null(covar)){
    # No covar is given
    num_covar <- 1
  }else if(length(covar) == 1){
    check_colname(colnames(df), covar, "df")
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
      lab <- df[[label]][1] # Current label group
      if(num_batches == 1){
        batch <- df$batch[1]
        message(paste("Label group", lab, "only contains cells from batch", batch))
        df <- df %>% select(-label) # Not removed from output, but removed here to prevent bug
        return(df)
      }
      message(paste("Correcting Label group", lab))
      # Calculate number of covars in the node
      if(!is.null(covar)){
        num_covar <- df[[covar]] %>%
          unique() %>%
          length()

        # If a node is heavily skewed to a single covar, it should be treated as having only 1 covar
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
                                                mod = stats::model.matrix(~as.factor(df[[covar]])),
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
        # Cap values to range of input data
        dplyr::mutate(dplyr::across(dplyr::all_of(markers),
                                    function(x) {
                                      min <- min(df[[dplyr::cur_column()]])
                                      max <- max(df[[dplyr::cur_column()]])
                                      x <- ifelse(x < min, min, x)
                                      x <- ifelse(x > max, max, x)
                                      return(x)
                                    })) %>%
      # Only add covar column, if it is not null
      if(!is.null(covar)) ComBat_output[[covar]] <- df[[covar]]

      return(ComBat_output)
    }) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(id) %>%
    select(id, everything())
  return(corrected_data)
}

#' Alternate correction
#'
#' This function allows running ComBat with a custom covar mod matrix.
#'  A model could look like \code{stats::model.matrix(~df$covar+df$label)}
#'
#' @inheritParams correct_data
#' @param mod Covariate model to use in ComBat.
#' @examples
#' corrected <- uncorrected %>%
#'   correct_data(mod = stats::model.matrix(~df$covar+df$label))
#' @export
correct_data_alt <- function(df,
                             mod,
                             markers = NULL,
                             parametric = TRUE){
  message("Batch correcting data..")
  if (is.null(markers)){
    # Get markers
    markers <- df %>%
      cyCombine::get_markers()
  }


  corrected_data <- df %>%
      # Compute ComBat correction
    dplyr::select(all_of(markers)) %>%
    t() %>%
    # The as.character is to remove factor levels not present in the SOM node
    sva::ComBat(batch = as.character(df$batch),
                mod = mod,
                par.prior = parametric) %>%
    t() %>%
    tibble::as_tibble() %>%
    dplyr::mutate(batch = df$batch,
                  sample = df$sample,
                  id = df$id) %>%
    # Cap values to range of input data
    dplyr::mutate(dplyr::across(dplyr::all_of(markers),
                                function(x) {
                                  min <- min(df[[dplyr::cur_column()]])
                                  max <- max(df[[dplyr::cur_column()]])
                                  x <- ifelse(x < min, min, x)
                                  x <- ifelse(x > max, max, x)
                                  return(x)
                                })) %>%
    select(id, everything())
  return(corrected_data)
}


# Wrapper ----

#' Run batch correction on data
#'
#' This is a wrapper function for the cyCombine batch correction workflow.
#'  To run the workflow manually, type "batch_correct" to see the source code of this wrapper and follow along.
#'  som_type = "fsom" is recommended when merging batches from a single study/experiment.
#'  som_type = "kohonen is recommended when merging data from different studies/experiments.
#'
#' @inheritParams create_som
#' @inheritParams correct_data
#' @inheritParams normalize
#' @family batch
#' @examples
#' corrected <- uncorrected %>%
#'   batch_correct(markers = markers,
#'   covar = "condition")
#' @export
batch_correct <- function(df,
                          label = NULL,
                          xdim = 8,
                          ydim = 8,
                          rlen = 10,
                          parametric = TRUE,
                          seed = 473,
                          covar = NULL,
                          markers = NULL,
                          norm_method = 'scale',
                          ties.method = "average"){
  # A batch column is required
  check_colname(colnames(df), "batch", "df")

  # Create SOM on scaled data
  if(is.null(label)) {
    label <- df %>%
      normalize(markers = markers,
                norm_method = norm_method,
                ties.method = ties.method) %>%
      create_som(markers = markers,
                 rlen = rlen,
                 seed = seed,
                 xdim = xdim,
                 ydim = ydim)
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
