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
#' @param ties.method The method to handle ties, when using rank. Default: 'average'. See ?rank for other options.
#' @family batch
#' @examples
#' \dontrun{
#' df_normed <- df %>%
#'   normalize()
#'   }
#' @export
normalize <- function(df,
                      markers = NULL,
                      norm_method = "scale",
                      ties.method = "average") {

  # Remove case-sensitivity
  norm_method <- norm_method %>% stringr::str_to_lower()
  ties.method <- ties.method %>% stringr::str_to_lower()

  # Error check
  if (norm_method == "rank" && ties.method %!in% c("average", "first", "last", "random", "max", "min")) {
    stop("When using norm_method = 'rank', please use an available ties.method (average, first, last, random, max, or min).")
  }

  # Messaging
  if (norm_method == "rank") {message("Ranking expression data..")
  } else if (norm_method == "scale") {message("Scaling expression data..")
  } else if (norm_method == "qnorm") {
    # message("Quantile normalizing expression data..")
    # Run quantile normalization
    df_normed <- cyCombine:::quantile_norm(df, markers = markers)
    return(df_normed)
  } else stop("Please use either 'scale', 'rank', or 'qnorm' as normalization method." )

  if (is.null(markers)) {
    # Get markers
    markers <- df %>%
      cyCombine::get_markers()
  }

  # Scale or rank at marker positions individually for every batch
  df_normed <- df %>%
    dplyr::group_by(.data$batch) %>%
    purrr::when(
      norm_method == "rank"  ~ dplyr::mutate(
        ., dplyr::across(dplyr::all_of(markers),
                         .fns = ~ {
                           if(sum(.x) == 0) stop("A marker is 0 for an entire batch. Please remove this marker.")
                           rank(.x, ties.method = ties.method) / length(.x)})),
      norm_method == "scale" ~ dplyr::mutate(., dplyr::across(dplyr::all_of(markers),
                                                              .fns = ~{
                                                                if(sum(.x) == 0) stop("A marker is 0 for an entire batch. Please remove this marker.")
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
#' \dontrun{
#' df_qnorm <- preprocessed %>%
#'   quantile_norm()
#'   }
quantile_norm <- function(df, markers = NULL) {
  message("Quantile normalizing expression data..")
  if (is.null(markers)) {
    # Get markers
    markers <- df %>%
      cyCombine::get_markers()
  }

  # Determine goal distributions for each marker by getting quantiles across all batches
  refq <- list()
  for (m in markers) {
    # Determine the quantiles
    refq[[m]] <- stats::quantile(unlist(df[,m]), probs=seq(0,1,length.out=5), names = F)
  }

  qnormed_expr <- df
  for (b in unique(df$batch)) {
    for (m in markers) {
      qx <- stats::quantile(unlist(df[df$batch == b, m]), probs=seq(0,1,length.out=5), names = F)
      spf <- stats::splinefun(x=qx, y=refq[[m]], method="monoH.FC", ties=min)

      # Apply the spline function to adjust quantiles
      qnormed_expr[qnormed_expr$batch == b, m] <- spf(unlist(df[df$batch==b, m]))

    }
  }

  return(qnormed_expr)
}



# Clustering ----

#' Create Self-Organizing Map
#'
#' The function uses the kohonen package to create a Self-Organizing Map.
#'  It is used to segregate the cells for the batch correction to make the correction less affected
#'  by samples with high abundances of a particular cell type.
#'
#' @inheritParams normalize
#' @param seed The seed to use when creating the SOM.
#' @param xdim The x-dimension size of the SOM.
#' @param ydim The y-dimension size of the SOM.
#' @param rlen Number of times the data is presented to the SOM network
#' @family batch
#' @examples
#' \dontrun{
#' labels <- uncorrected %>%
#'   create_som()
#'   }
#' @export
#' @return A vector of clustering labels
create_som <- function(df,
                       markers = NULL,
                       seed = 473,
                       rlen = 10,
                       xdim = 8,
                       ydim = 8) {
  if (is.null(markers)) {
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


# Batch correction ----

#' Correct data using ComBat
#'
#' Compute the batch correction on the data using the ComBat algorithm.
#'  Define a covariate, either as a character vector or name of tibble column.
#'  The covariate should preferable be the cell condition types but can be any column that infers heterogeneity in the data.
#'  The function assumes that the batch information is in the "batch" column and the data contains a "sample" column with sample information.
#'
#' @param label The cluster or cell type label. Either as a column name or vector.
#' @param covar The covariate ComBat should use. Can be a vector or a column name in the input tibble.
#'   If NULL, no covar will be used
#' @param anchor Experimental: A column or vector specifying which samples are replicates and which are not. If specified, this column will be used as a covariate in ComBat. Be aware that it may be confounded with the condition.
#' @param parametric Default: TRUE. If TRUE, the parametric version of ComBat is used. If FALSE, the non-parametric version is used.
#' @inheritParams normalize
#' @family batch
#' @examples
#' \dontrun{
#' corrected <- uncorrected %>%
#'   correct_data(label = labels, covar = "condition")
#'   }
#' @export
correct_data <- function(df,
                         label,
                         covar = NULL,
                         anchor = NULL,
                         markers = NULL,
                         parametric = TRUE) {
  message("Batch correcting data..")
  # Check for batch column
  cyCombine:::check_colname(colnames(df), "batch", "df")
  if (is.null(markers)) {
    # Get markers
    markers <- df %>%
      cyCombine::get_markers()
  }

  # Add ID column to retain data order
  if("id" %!in% colnames(df)) df$id <- seq_len(nrow(df))

  # Add label to df
  if (length(label) == 1) {
    cyCombine:::check_colname(colnames(df), label, "df")
  } else {
    df$label <- label
    label <- "label"
  }

  # Add covar to df, if given
  if (!is.null(covar)) {
    if (length(covar) == 1) {
      cyCombine:::check_colname(colnames(df), covar, "df")
      df[[covar]] <- as.factor(df[[covar]])
    } else {
      # Covar was given as a vector
      df$covar <- as.factor(covar)
      covar <- "covar"
    }
    # Ensure there is more than 1 factor level
    if (nlevels(df[[covar]]) == 1) covar <- NULL
  }

  # Add anchor to df, if given
  if (!is.null(anchor)) {
    if (length(anchor) == 1) {
      cyCombine:::check_colname(colnames(df), anchor)
      df[[anchor]] <- as.factor(df[[anchor]])
    } else {
      # Anchor was given as a vector
      df$anchor <- as.factor(anchor)
      anchor <- "anchor"
    }
    # Ensure there is more than 1 factor level
    if (nlevels(df[[anchor]]) == 1) anchor <- NULL
  }

  corrected_data <- df %>%
    dplyr::group_by(.data[[label]]) %>%
    # Correct (modify) each label group with ComBat
    dplyr::group_modify(.keep = TRUE, function(df, ...) {
      # Initiate anchor and covar counter
      num_covar <- 1
      num_anchor <- 1
      # Detect if only one batch is present in the node
      num_batches <- df$batch %>%
        factor() %>%
        nlevels()
      lab <- df[[label]][1] # Current label group
      if (num_batches == 1) {
        batch <- df$batch[1]
        message(paste("Label group", lab, "only contains cells from batch", batch))
        df <- df %>% dplyr::select(-label) # Not removed from output, but removed here to prevent bug
        return(df)
      }
      message(paste("Correcting Label group", lab))
      # Calculate number of covars in the node
      if (!is.null(covar)) {

        # Only use covar, if it does not confound with batch
        if (!cyCombine:::check_confound(df$batch, stats::model.matrix(~df[[covar]]))) {
          num_covar <- df[[covar]] %>%
            factor() %>%
            nlevels()

          # If a node is heavily skewed to a single covar, it should be treated as having only 1 covar.
          # Get number of cells in the condition with most cells
          covar_counts <- df %>%
            dplyr::count(.data[[covar]]) %>%
            dplyr::pull(n)

          if (sum(covar_counts) < max(covar_counts) + num_covar*5) {
            message("The label group almost exclusively consists of cells from a single covar. Therefore, covar is ignored for this label group")
            num_covar <- 1
          }
        } else {
          message("Covar is confounded with batch. Ignoring covar in this label group")
        }
      }
      # Do a similar check on anchor
      if (!is.null(anchor)) {
        if (!cyCombine:::check_confound(df$batch, stats::model.matrix(~df[[anchor]]))) {
          num_anchor <- df[[anchor]] %>%
            factor() %>%
            nlevels()

          # If a node is heavily skewed to a single anchor, it should be treated as having only 1 covar.
          # Get number of cells in the anchor with most cells
          anchor_counts <- df %>%
            dplyr::count(.data[[anchor]]) %>%
            dplyr::pull(n)

          if (sum(anchor_counts) < max(anchor_counts) + num_anchor*5) {
            message("The label group almost exclusively consists of cells from a single anchor group. Therefore, anchor is ignored for this label group")
            num_anchor <- 1
          }
        } else {
          message("Anchor is confounded with batch. Ignoring anchor in this label group")
        }
      }
      if (num_covar > 1 & num_anchor > 1) {
        # If neither covar nor anchor confounds with batch but they do each other, prioritise covar
        if (cyCombine:::check_confound(df$batch, stats::model.matrix(~df[[covar]] + df[[anchor]]))) {
          num_anchor <- 1
          message("Anchor and covar are confounded. Ignoring anchor in this label group")
        }
      }
      # Compute ComBat correction
      ComBat_output <- df %>%
        dplyr::select(dplyr::all_of(markers)) %>%
        t() %>%
        # The as.character is to remove factor levels not present in the SOM node
        purrr::when(num_covar > 1 & num_anchor == 1 ~
                      sva::ComBat(.,
                                  batch = as.character(df$batch),
                                  mod = stats::model.matrix(~ df[[covar]]),
                                  par.prior = parametric),
                    num_covar > 1 & num_anchor > 1 ~
                      sva::ComBat(.,
                                  batch = as.character(df$batch),
                                  mod = stats::model.matrix(~df[[covar]] + df[[anchor]]),
                                  par.prior = parametric),
                    num_covar == 1 & num_anchor > 1 ~
                      sva::ComBat(.,
                                  batch = as.character(df$batch),
                                  mod = stats::model.matrix(~df[[anchor]]),
                                  par.prior = parametric),
                    num_covar == 1 & num_anchor == 1 ~
                      sva::ComBat(.,
                                  batch = as.character(df$batch),
                                  par.prior = parametric)
        ) %>%
        t() %>%
        tibble::as_tibble() %>%
        dplyr::bind_cols(
          dplyr::select(df,
                        -dplyr::all_of(c(markers, label)))) %>%
        # Cap values to range of input data
        dplyr::mutate(dplyr::across(dplyr::all_of(markers),
                                    function(x) {
                                      min <- min(df[[dplyr::cur_column()]])
                                      max <- max(df[[dplyr::cur_column()]])
                                      x <- ifelse(x < min, min, x)
                                      x <- ifelse(x > max, max, x)
                                      return(x)
                                    }))
      # Only add covar column, if it is not null
      if (!is.null(covar)) ComBat_output[[covar]] <- df[[covar]]
      # Only add anchor column, if it is not null
      if (!is.null(anchor)) ComBat_output[[anchor]] <- df[[anchor]]

      return(ComBat_output)
    }) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(id) %>%
    dplyr::select(id, dplyr::everything()) %>%
    dplyr::mutate(batch = as.factor(batch))
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
#' \dontrun{
#' corrected <- uncorrected %>%
#'   correct_data(mod = stats::model.matrix(~df$covar+df$label))
#'   }
correct_data_alt <- function(df,
                             mod,
                             markers = NULL,
                             parametric = TRUE) {
  message("Batch correcting data..")
  if (is.null(markers)){
    # Get markers
    markers <- df %>%
      cyCombine::get_markers()
  }


  corrected_data <- df %>%
      # Compute ComBat correction
    dplyr::select(dplyr::all_of(markers)) %>%
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
    dplyr::select(id, dplyr::everything())
  return(corrected_data)
}


# Wrapper ----

#' Run batch correction on data
#'
#' This is a wrapper function for the cyCombine batch correction workflow.
#'  To run the workflow manually, type "batch_correct" to see the source code of this wrapper and follow along or read the vignettes on the GitHub page \url{https://github.com/biosurf/cyCombine}.
#'
#' @inheritParams create_som
#' @inheritParams correct_data
#' @inheritParams normalize
#' @family batch
#' @examples
#' \dontrun{
#' corrected <- uncorrected %>%
#'   batch_correct(markers = markers,
#'   covar = "condition")
#'   }
#' @export
batch_correct <- function(df,
                          label = NULL,
                          xdim = 8,
                          ydim = 8,
                          rlen = 10,
                          parametric = TRUE,
                          seed = 473,
                          covar = NULL,
                          anchor = NULL,
                          markers = NULL,
                          norm_method = "scale",
                          ties.method = "average") {
  # A batch column is required
  cyCombine:::check_colname(colnames(df), "batch", "df")
  if (any(is.na(df$batch))) { # Check for NAs
    message("Some batches contain NAs. These will be removed")
    warning("Some batches contain NAs. These will be removed")
    df <- df %>%
      dplyr::filter(!is.na(batch))
    }

  # Create SOM on scaled data
  if (is.null(label)) {
    label <- df %>%
      cyCombine::normalize(markers = markers,
                           norm_method = norm_method,
                           ties.method = ties.method) %>%
      cyCombine::create_som(markers = markers,
                            rlen = rlen,
                            seed = seed,
                            xdim = xdim,
                            ydim = ydim)
  }


  # Run batch correction
  corrected <- df %>%
    cyCombine::correct_data(label = label,
                            covar = covar,
                            anchor = anchor,
                            markers = markers,
                            parametric = parametric)
  message("Done!")
  return(corrected)
}
