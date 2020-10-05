#### Batch correction ####


#' Batch-wise scaling of data
#'
#' @family batch
#' @export
scale_expr <- function(input){
  cat("Scaling expression data\n")
  scaled_expr <- input %>%
    dplyr::group_by(batch) %>%
    dplyr::mutate_at(dplyr::vars(-c("batch", "sample", "id")), .funs = scale) %>%
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
                       seed = 473,
                       xdim = 10,
                       ydim = 10){
  # 10x10 SOM grid on overlapping markers, extract clustering per cell
  cat("Creating SOM grid\n")
  set.seed(seed)
  som <- scaled_expr %>%
    dplyr::select(-c("batch", "sample", "id")) %>%
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
    covar[grep('^HD', data_subset$sample)] <- 'HD'
    magic_output <- data_subset %>%
      dplyr::select(-c("batch", "sample", "id")) %>%
      t() %>%
      sva::ComBat(batch = data_subset$batch, mod = model.matrix(~covar)) %>%
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
#' @importFrom sva ComBat
#' @family batch
#' @export
correct_data <- function(input,
                         som_classes){

  corrected_data <- input %>%
    dplyr::mutate(som = som_classes,
                  # Determine covariate
                  covar = case_when(stringr::str_starts(string = sample,
                                                        pattern = "HD") ~ "HD",
                                    TRUE ~ "CLL") %>%
                    as.factor()) %>%
    dplyr::group_by(som) %>%
    # Run ComBat on each SOM class
    dplyr::group_modify(function(df, ...){

      ComBat_output <- df %>%
        dplyr::select_if(colnames(.) %!in% c("batch", "sample", "covar", "id")) %>%
        t() %>%
        # The as.character is to remove factor levels not present in the SOM node
        sva::ComBat(batch = as.character(df$batch), mod = model.matrix(~df$covar)) %>%
        t() %>%
        tibble::as_tibble() %>%
        dplyr::mutate(batch = df$batch,
               sample = df$sample,
               covar = df$covar,
               id = df$id)
      return(ComBat_output)
    }) %>%
    dplyr::ungroup() %>%
    dplyr::select(-som) %>%
    # Reduce all negative values to zero
    dplyr::mutate_at(vars(-c("batch", "sample", "covar", "id")),#, "som")),
                     function(x) {
                       x[x < 0] <- 0
                       return(x)
                       }) %>%
    dplyr::arrange(id)
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
  corrected <- preprocessed %>%
    correct_data(som_classes = som$unit.classif)
  cat("Done!\n")
  return(corrected)
}





