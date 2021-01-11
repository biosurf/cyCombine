#### Panel merging ####


#' Salvaging problematic channels - one at a time
#'
#' This function imputes the expression values for a single channel based on the information in the other channels
#'   The purpose is to be able to salvage a misstained channel in single batches of an experiment
#'
#' @param df Dataframe with expression values
#' @param correct_batches Batches with misstaining (or other problem) for which to perform the imputation of 'channel'
#' @param channel Channel (marker) to impute
#' @param sample_size Number of cells to base imputation on - defaults to all cells in df not in correct_batches
#' @param exclude Channels to exclude (e.g. other misstained channels not to be used for imputation)
#' @family merging
#' @examples
#' df_fixed <- df %>%
#'   salvage_problematic(correct_batches = 1, channel = 'CD127')
#' @export
salvage_problematic <- function(df, correct_batches, channel, sample_size = NULL, exclude = NULL) {
  
  if (is.null(sample_size)) {
    sample_size = nrow(df)
  }
  
  print(paste('Started imputation for', channel, 'in batch(es)', paste(correct_batches, collapse = ', ')))
  
  # Get all other channels in samples to impute for
  impute_for <- df[df$batch %in% correct_batches, !colnames(df) == channel]
  
  # Get the observations not in problematic batch(es)
  complete_obs <- df[!(df$batch %in% correct_batches),]
  complete_obs <- complete_obs[sample(nrow(complete_obs),sample_size),]
  
  # Get 8 x 8 SOM classes for each complete event - this takes some time, which is why I downsample
  print('Calculating SOM')
  som_res <- som(as.matrix(complete_obs[,!colnames(complete_obs) %in% c(channel, "batch", "sample", "id", exclude)]), 
                 grid=somgrid(xdim = 8, ydim = 8), 
                 dist.fcts = "euclidean")
  
  som_classes <- som_res$unit.classif
  
  # For each complete cluster
  # Calculate the distance from each missing event to centroid (calculated based on mean) of each cluster
  # AND
  # Calculate the density distribution of the channel to be imputed
  clustdist <- matrix(nrow = nrow(impute_for), ncol = max(som_classes))
  densList <- list()
  for (s in sort(unique(som_classes))) {
    
    clustdist[,s] <- dista(x = as.matrix(impute_for[,!colnames(impute_for) %in% c("batch", "sample", "id", exclude)]), 
                           xnew = matrix(colMeans(complete_obs[which(som_classes==s),!colnames(complete_obs) %in% c(channel, "batch", "sample", "id", exclude)]), nrow = 1), 
                           type = "euclidean")
    
    densList[[s]] <- complete_obs[which(som_classes==s),] %>% 
      pull(channel) %>% 
      density()
    
  }
  
  # Assign cluster to each test event by finding closest centroid
  missingCluster <- sapply(1:nrow(impute_for), function(x) {which.min(clustdist[x,])})
  
  # For each missing event in each cluster, impute values based on density draws in respective cluster (here represented by random sample and addition of a number with mean 0 and bandwidth sd)
  print('Performing density draws')
  imputed <- rep(0, length(missingCluster))
  for (s in sort(unique(missingCluster))) {
    
    imputed[missingCluster==s] <- complete_obs[which(som_classes==s),] %>%
      pull(channel) %>%
      sample(size = sum(missingCluster==s), replace = T) + rnorm(sum(missingCluster==s), 0, densList[[s]]$bw)
    
    # imputed[missingCluster==s] <- sample(as.numeric(unlist(complete_obs[which(som_classes==s),channel])), length(which(missingCluster==s)), replace = T) + rnorm(length(which(missingCluster==s)), 0, densList[[s]]$bw)
  }
  # Fix values below 0
  imputed[imputed < 0] <- 0
  
  # Put the values back in the dataset
  df[df$batch %in% correct_batches, channel] <- imputed
  
  return(df)
}




#' Impute non-overlapping channels for whole data sets
#'
#' This function imputes the expression values for a whole dataset based on the overlapping markers contained in another dataset
#'   The purpose is to be able to merge multi-panel data or to merge differently run datasets which have non-overlapping markers
#'
#' @param impute_for Dataframe with expression values to impute markers for
#' @param complete_obs Dataframe with expression values to base imputation on
#' @param overlap_channels Channels (markers) that overlap between impute_for and complete_obs, which can be used to base imputation on
#' @param impute_channels Channels to impute and add to the impute_for set (must be present in complete_obs)
#' @family merging
#' @examples
#' df_imputed <- df_panel1 %>%
#'   impute_across_panels(complete_obs = df_panel2, overlap_channels = intersect(df_panel1_markers, df_panel2_markers), impute_channels = unique_df_panel2_markers)
#' @export
impute_across_panels <- function(impute_for, complete_obs, overlap_channels, impute_channels) {
  
  print(paste('Started imputation for', paste(impute_channels, collapse = ', ')))
  
  
  # Checking colnames
  if (!all(impute_channels %in% complete_obs)) {
    stop("Error: Some of your impute_channels are not found among the complete_obs column names.")
  }
  if (!all(overlap_channels %in% impute_for)) {
    stop("Error: Some of your overlap_channels are not found among the impute_for column names.")
  }
  if (!all(overlap_channels %in% complete_obs)) {
    stop("Error: Some of your overlap_channels are not found in among complete_obs column names.")
  }
  
  
  panel_imputed <- NULL
  
  # Get 8 x 8 SOM classes for complete_obs and impute_for on overlapping channels - takes some time...
  print('Calculating SOM')
  overlapping_data <- as.matrix(rbind(complete_obs[,overlap_channels], impute_for[,overlap_channels]))
  som_res <- som(overlapping_data,
                 grid=somgrid(xdim = 8, ydim = 8), 
                 dist.fcts = "euclidean")
  
  # Extract and split SOM classes
  som_classes <- som_res$unit.classif
  complete_obs_som <- som_classes[1:nrow(complete_obs)]
  impute_obs_som <- som_classes[(nrow(complete_obs)+1):length(som_classes)]
  
  
  # For each missing event in each cluster, impute values based on density draws in the same cluster
  print('Performing density draws')
  imputed <- matrix(nrow=nrow(impute_for), ncol=length(impute_channels), dimnames = list("rn"=NULL, "cn"=impute_channels))
  
  for (s in sort(unique(impute_obs_som))) {
    
    # Check if there's at least 50 cells from the 'training' data in the given cluster and if so, impute - otherwise set values to NA and throw a warning
    if (sum(complete_obs_som==s) >= 50) {
      
      # Estimate the density per marker that needs imputation
      dens <- apply(complete_obs[which(complete_obs_som==s), impute_channels], 2, function(x) {density(x)$bw})
      
      # Performing the imputation
      imputed[which(impute_obs_som==s),] <- complete_obs[which(complete_obs_som==s),] %>%
        select(all_of(impute_channels)) %>%
        sample_n(size = sum(impute_obs_som==s), replace = T) %>% 
        as.matrix() + sapply(impute_channels, function(ch) {rnorm(length(which(impute_obs_som==s)), 0, dens[ch])})
    } else {
      
      imputed[which(impute_obs_som==s),impute_channels] <- NA
      warning('Be aware that a cluster contains cells only from the dataset you wish to impute for. As a result, imputations were not made for those cells.')
    }
  }
  
  # Fix values below 0
  imputed[imputed < 0] <- 0
  
  # Add new columns to impute_for
  impute_for[,impute_channels] <- imputed
  
  impute_for <- impute_for %>%
    relocate(c(batch, sample, id), .after = all_of(impute_channels))
  
  return(impute_for)
}

