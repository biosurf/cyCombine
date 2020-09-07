#### Batch correction ####


#' Batch-wise scaling of data
#' @importFrom dplyr mutate_at group_by ungroup vars
#' @family batch
#' @export
scale_expr <- function(input){
  print("Scaling expression data")
  scaled_expr <- input %>%
    group_by(Batch) %>%
    mutate_at(vars(-c("Batch", "Sample")), .funs = scale) %>%
    ungroup()
  return(scaled_expr)
}


#' Create SOM
#'
#' @importFrom dplyr select all_of
#' @importFrom kohonen som somgrid
#'
#' @family batch
#' @export
create_som <- function(scaled_expr,
                       seed = 548,
                       xdim = 10,
                       ydim = 10){
  # 10x10 SOM grid on overlapping markers, extract clustering per cell
  print("Creating SOM grid")
  set.seed(seed)
  som <- scaled_expr %>%
    select(-c("Batch", "Sample")) %>%
    as.matrix() %>%
    kohonen::som(grid = kohonen::somgrid(xdim = xdim,
                                         ydim = ydim),
                 dist.fcts = "euclidean")
  return(som)
}




#' Correct data using ComBat
#'
#' @importFrom tibble tibble add_column
#' @importFrom purrr when
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
    mutate(covar = "")
  # corrected_data <- tibble::tibble(.rows = nrow(input)) %>%
  #   tibble::add_column(!!!set_names(as.list(rep(0, ncol(input))), nm = markers)) %>%
  #   mutate(Batch = 0,
  #          Sample = "",
  #          covar = "")


  for (s in sort(unique(som_classes))) {
    # Extract original (non-scaled+ranked) data for cluster
    data_subset <- input[which(som_classes==s), ]
    data_subset %>% count(Batch) %>%
      print()
    print(s)


    # ComBat batch correction using disease status as covariate
    covar <- rep('CLL', nrow(data_subset))
    covar[grep('^HD', data_subset$Sample)] <- 'HD'
    magic_output <- data_subset %>%
      select(-c("Batch", "Sample")) %>%
      t() %>%
      sva::ComBat(batch = data_subset$Batch, mod = model.matrix(~covar)) %>%
      t() %>%
      as_tibble() %>%
      mutate(Batch = data_subset$Batch,
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
    mutate(som = som_classes,
           covar = case_when(str_starts(string = Sample,
                                        pattern = "HD") ~ "HD",
                             TRUE ~ "CLL")) %>%
    # nest_by(som) %>%
    # mutate(mean = mean(data$CD20))
    group_by(som) %>%
    group_modify(function(df, ...){
      magic_output <- df %>%
        select(markers) %>%
        t() %>%
        sva::ComBat(batch = df$Batch, mod = model.matrix(~df$covar)) %>%
        t() %>%
        as_tibble() %>%
        mutate(Batch = df$Batch,
               Sample = df$Sample,
               covar = df$covar)
      magic_output[magic_output < 0] <- 0
      return(magic_output)
    }) %>%
    ungroup() %>%
    select(-som)
  return(corrected_data2)
}



#' Run batch correction on preprocessed data
#'
#'
#' @family batch
#' @export
batch_correct <- function(preprocessed_data,
                          markers){


  # Create SOM on scaled data
  som <- preprocessed_data %>%
    scale_expr() %>%
    create_som()
  # Run batch correction
  print("Batch correcting data")
  corrected_data <- preprocessed_data %>%
    correct_data(som_classes = som$unit.classif)
  print("Done")
  return(corrected_data)
}





