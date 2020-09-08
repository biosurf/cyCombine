#### Batch correction ####


#' Batch-wise scaling of data
#'
#' @family batch
#' @export
scale_expr <- function(input){
  cat("Scaling expression data")
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
  cat("Creating SOM grid")
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
#' @importFrom sva ComBat
#' @family batch
#' @export
correct_data <- function(input,
                         som_classes){
  # sample_ids <- combined_expr$sample_ids
  # batch_ids <- combined_expr$batch_ids
  # Processing and correcting data per-cluster
  #combined_expr <- combined

  # Create empty dataset
  corrected_data <- input %>%
    dplyr::mutate(covar = "")
  # corrected_data <- tibble::tibble(.rows = nrow(input)) %>%
  #   tibble::add_column(!!!set_names(as.list(rep(0, ncol(input))), nm = markers)) %>%
  #   mutate(Batch = 0,
  #          Sample = "",
  #          covar = "")


  for (s in sort(unique(som_classes))) {
    # Extract original (non-scaled+ranked) data for cluster
    data_subset <- input[which(som_classes==s), ]
    cat("som class:", s, sep = " ")


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

    corrected_data[which(som_classes==s), ] <- magic_output
  }
  # Fix values below zero
  corrected_data[corrected_data < 0] <- 0

  return(corrected_data)
}


correct_data2 <- function(input,
                          som_classes,
                          markers){
  # sample_ids <- combined_expr$sample_ids
  # batch_ids <- combined_expr$batch_ids
  # Processing and correcting data per-cluster
  #combined_expr <- combined

  # Create empty dataset
  # corrected_data <- input %>%
  #   mutate(covar = "")

  corrected_data2 <- input %>%
    dplyr::mutate(som = som_classes,
           covar = case_when(stringr::str_starts(string = Sample,
                                        pattern = "HD") ~ "HD",
                             TRUE ~ "CLL")) %>%
    # nest_by(som) %>%
    # mutate(mean = mean(data$CD20))
    dplyr::group_by(som) %>%
    dplyr::group_modify(function(df, ...){
      magic_output <- df %>%
        dplyr::select(markers) %>%
        t() %>%
        sva::ComBat(batch = df$Batch, mod = model.matrix(~df$covar)) %>%
        t() %>%
        tibble::as_tibble() %>%
        dplyr::mutate(Batch = df$Batch,
               Sample = df$Sample,
               covar = df$covar)
      magic_output[magic_output < 0] <- 0
      return(magic_output)
    }) %>%
    dplyr::ungroup() %>%
    dplyr::select(-som)
  return(corrected_data2)
}



#' Run batch correction on preprocessed data
#'
#'
#' @family batch
#' @export
batch_correct <- function(preprocessed_data){


  # Create SOM on scaled data
  som <- preprocessed_data %>%
    scale_expr() %>%
    create_som()
  # Run batch correction
  cat("Batch correcting data")
  corrected_data <- preprocessed_data %>%
    correct_data(som_classes = som$unit.classif)
  cat("Done!")
  return(corrected_data)
}





