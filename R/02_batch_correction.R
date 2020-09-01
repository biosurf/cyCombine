#### Batch correction ####


#' Batch-wise scaling of data
#' @importFrom dplyr mutate_at group_by ungroup vars
#' @family batch
scale_expr <- function(combined_expr){
  scaled_expr <- combined_expr %>%
    group_by(batch_ids) %>%
    mutate_at(vars(-c("batch_ids", "sample_ids")), .funs = scale) %>%
    ungroup()
  return(scaled_expr)
}


#' Create SOM
#'
#' @importFrom dplyr select all_of
#' @importFrom kohonen som somgrid
#'
#' @family batch
create_som <- function(scaled_expr,
                       seed = 548,
                       xdim = 10,
                       ydim = 10){
  # 10x10 SOM grid on overlapping markers, extract clustering per cell
  set.seed(seed)
  som <- scaled_expr %>%
    select(all_of(all_markers)) %>%
    as.matrix() %>%
    kohonen::som(grid = kohonen::somgrid(xdim = xdim,
                                         ydim = ydim),
                 dist.fcts = "euclidean")
  return(som)
}




#' Correct data using ComBat
#'
#' @importFrom tibble tibble add_column
#' @importFrom sva ComBat
#' @family batch
correct_data <- function(combined_expr,
                         som_classes,
                         markers){
  # sample_ids <- combined_expr$sample_ids
  # batch_ids <- combined_expr$batch_ids
  # Processing and correcting data per-cluster
  #combined_expr <- combined

  # Create empty dataset
  corrected_data <- tibble::tibble(.rows = nrow(combined_expr)) %>%
    tibble::add_column(!!!set_names(as.list(rep(0, length(markers))), nm = markers))


  for (s in sort(unique(som_classes))) {

    # Extract original (non-scaled+ranked) data for cluster
    data <- combined_expr[which(som_classes==s), ]


    # ComBat batch correction using disease status as covariate
    covar <- rep('CLL', nrow(data))
    covar[grep('^HD', data$sample_ids)] <- 'HD'
    magic_output <- data %>%
      select(-c("batch_ids", "sample_ids")) %>%
      t() %>%
      sva::ComBat(batch = data$batch_ids, mod = model.matrix(~covar)) %>%
      t()
      # t(sva::ComBat(t(data), batch = batches, mod = model.matrix(~covar)))

    corrected_data[which(som_classes==s), ] <- magic_output
  }
  # Fix values below zero
  corrected_data[corrected_data < 0] <- 0

  return(corrected_data)
}


#' Run batch correction on preprocessed data
#'
#'
#' @family batch
batch_correct <- function(data,
                          markers){

  # Create SOM on scaled data
  som <- data %>%
    scale_expr() %>%
    create_som()
  # Run batch correction
  corrected_data <- combined_expr %>%
    correct_data(som_classes = som$unit.classif,
                 markers = markers)
  return(corrected_data)
}





### Plotting ----
if (FALSE){
  # Plot densities - uncorrected vs. corrected
  density_plots(uncorrected = select(combined_expr, -c("batch_ids","sample_ids")),
                corrected = corrected_data,
                combined_expr$batch_ids,
                filename = 'figs/02_densities_withcovar.png')

  # PCA plot uncorrected
  pca1 <- combined_expr %>% select(-batch_ids) %>%
    dimred_plot(combined_expr$batch_ids, 'uncorrected', type = 'pca')


  # UMAP plot uncorrected
  umap1 <- combined_expr %>% select(-batch_ids) %>%
    dimred_plot(combined_expr$batch_ids, 'uncorrected', type = 'umap')



  # PCA plot corrected
  pca2 <- corrected_data %>%
    dimred_plot(combined_expr$batch_ids, 'corrected', type = 'pca')
  save_two_plots(pca1, pca2, filename = 'figs/02_pca.png')


  # UMAP plot corrected
  umap2 <- dimred_plot(corrected_data, combined_expr$batch_ids, 'corrected', type = 'umap')
  save_two_plots(umap1, umap2, filename = 'figs/02_umap.png')

}
