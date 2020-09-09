#### Batch correction ####


#' Batch-wise scaling of data
#'
#' @family batch
#' @export
scale_expr <- function(input){
  cat("Scaling expression data\n")
  scaled_expr <- input %>%
    dplyr::group_by(Batch) %>%
    dplyr::mutate_at(dplyr::vars(-c("Batch", "Sample")), .funs = scale) %>%
    dplyr::ungroup()
  return(scaled_expr)
}


#' Create SOM
#'
#' @importFrom kohonen som somgrid
#'
#' @family batch
#' @export
create_som <- function(scaled_expr,
                       seed = 548,
                       xdim = 10,
                       ydim = 10){
  # 10x10 SOM grid on overlapping markers, extract clustering per cell
  cat("Creating SOM grid\n")
  set.seed(seed)
  som <- scaled_expr %>%
    dplyr::select(-c("Batch", "Sample")) %>%
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
#' @family batch
correct_data_prev <- function(input,
                         som_classes){
  # Create copy dataset
  corrected_data <- input %>%
    dplyr::mutate(covar = "")


  for (s in sort(unique(som_classes))) {
    # Extract original (non-scaled+ranked) data for cluster
    data_subset <- input[which(som_classes==s), ]
    cat("som class:", s, "\n", sep = " ")


    # ComBat batch correction using disease status as covariate
    covar <- rep('CLL', nrow(data_subset))
    covar[grep('^HD', data_subset$Sample)] <- 'HD'
    magic_output <- data_subset %>%
      dplyr::select(-c("Batch", "Sample")) %>%
      t() %>%
      sva::ComBat(batch = data_subset$Batch, mod = model.matrix(~covar)) %>%
      t() %>%
      tibble::as_tibble() %>%
      dplyr::mutate(Batch = data_subset$Batch,
             Sample = data_subset$Sample,
             covar = covar)

    # Fill copy dataset with corrected data
    corrected_data[which(som_classes==s), ] <- magic_output
  }
  # Fix values below zero
  corrected_data[corrected_data < 0] <- 0
  corrected_data <- corrected_data %>%
    arrange(Batch)

  return(corrected_data)
}

#' Correct data using ComBat
#'
#' @importFrom sva ComBat
#' @family batch
#' @export
correct_data <- function(input,
                         som_classes){


  corrected_data <- input %>%
    dplyr::mutate(som = som_classes,
                  # Determine covariate
                  covar = case_when(stringr::str_starts(string = Sample,
                                                        pattern = "HD") ~ "HD",
                                    TRUE ~ "CLL")) %>%
    dplyr::group_by(som) %>%
    # Run ComBat on each SOM class
    dplyr::group_modify(function(df, ...){
      ComBat_output <- df %>%
        dplyr::select(-c("Batch", "Sample", "covar")) %>%
        t() %>%
        sva::ComBat(batch = df$Batch, mod = model.matrix(~df$covar)) %>%
        t() %>%
        tibble::as_tibble() %>%
        dplyr::mutate(Batch = df$Batch,
               Sample = df$Sample,
               covar = df$covar)
      return(ComBat_output)
    }) %>%
    dplyr::ungroup() %>%
    dplyr::select(-som) %>%
    # Reduce all negative values to zero
    dplyr::mutate_at(vars(-c("Batch", "Sample", "covar")),
                     function(x) {
                       x[x < 0] <- 0
                       return(x)
                       }) %>%
    dplyr::arrange(Batch)
  return(corrected_data)
}



#' Run batch correction on preprocessed data
#'
#'
#' @family batch
#' @export
batch_correct <- function(preprocessed,
                          xdim = 10,
                          ydim = 10,
                          seed = 473){


  # Create SOM on scaled data
  som <- preprocessed %>%
    scale_expr() %>%
    create_som(seed = seed,
               xdim = xdim,
               ydim = ydim)
  # Run batch correction
  cat("Batch correcting data\n")
  corrected_data <- preprocessed %>%
    correct_data(som_classes = som$unit.classif)
  cat("Done!")
  return(corrected_data)
}





