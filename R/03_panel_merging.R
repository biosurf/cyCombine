#### Panel merging ####


#' Salvaging problematic channels - one at a time
#'
#' This function imputes the expression values for a single channel based on the information in the other channels.
#'   The purpose is to be able to salvage a misstained channel in single batches of an experiment
#'
#' @param df Dataframe with expression values
#' @param correct_batches Batches with misstaining (or other problem) for which to perform the imputation of 'channel'
#' @param channel Channel (marker) to impute
#' @param label (Optional) Column with group labels to base imputations on
#' @param sample_size Number of cells to base imputation on - defaults to all cells in df not in 'correct_batches'
#' @param exclude Channels to exclude (e.g. other misstained channels not to be used for imputation)
#' @inheritParams create_som
#' @family merging
#' @examples
#' \dontrun{
#' df_fixed <- df %>%
#'   salvage_problematic(correct_batches = 1, channel = 'CD127')
#'   }
#' @export
salvage_problematic <- function(
    df,
    correct_batches,
    channel,
    markers = NULL,
    label = NULL,
    sample_size = NULL,
    exclude = NULL,
    xdim = 8,
    ydim = 8,
    seed = 482) {


  cyCombine:::check_colname(colnames(df), label, "df")

  message(paste("Started imputation for", channel, "in batch(es)", paste(correct_batches, collapse = ", ")))

  # Get all other channels in samples to impute for
  impute_for <- df %>%
    dplyr::select(-dplyr::all_of(channel)) %>%
    dplyr::filter(batch %in% correct_batches)

  # Get the observations not in problematic batch(es)
  complete_obs <- df %>%
    dplyr::filter(batch %!in% correct_batches)

  if (!is.null(sample_size)) {
    complete_obs <- complete_obs[sample(nrow(complete_obs),sample_size), ]
  }

  # Use global non_markers if available
  if(!is.null(.GlobalEnv$non_markers)) non_markers <- .GlobalEnv$non_markers
  # Combine the data to impute for and from based on the overlapping markers
  overlapping_data <- rbind(complete_obs[,!colnames(complete_obs) %in% c(channel, non_markers, exclude)],
                            impute_for[,!colnames(impute_for) %in% c(non_markers, exclude)])


  # Get SOM classes for each event - on overlapping markers
  if (is(label, "NULL")){
    som_classes <- overlapping_data %>%
      cyCombine::create_som(markers = markers,
                            xdim = xdim,
                            ydim = ydim,
                            seed = seed)
  } else{
    som_classes <- overlapping_data[[label]]
  }


  # Split SOM classes to each original dataset
  complete_obs_som <- som_classes[1:nrow(complete_obs)]
  impute_obs_som <- som_classes[(nrow(complete_obs)+1):length(som_classes)]


  # Calculate the density distribution of the channel to be imputed
  # AND
  # For each missing event in each cluster, impute values based on density draws in respective cluster
  # (here represented by random sample and addition of a number with mean 0 and bandwidth sd)

  message("Performing density draws.")

  imputed <- rep(0, length(impute_obs_som))

  for (s in sort(unique(impute_obs_som))) {

    # Check if there's at least 50 cells from the 'training' data in the given cluster and if so, impute - otherwise set values to NA and throw a warning
    if (sum(complete_obs_som==s) >= 50) {

      # Estimate the density per marker that needs imputation
      dens <- complete_obs[which(complete_obs_som==s), ] %>%
        dplyr::pull(channel) %>%
        stats::density()

      # For each missing event in each cluster, impute values based on density draws in respective cluster (here represented by random sample and addition of a number with mean 0 and bandwidth sd)
      set.seed(seed)
      imputed[which(impute_obs_som==s)] <- complete_obs[which(complete_obs_som==s), ] %>%
        dplyr::pull(channel) %>%
        sample(size = sum(impute_obs_som==s), replace = T) + stats::rnorm(sum(impute_obs_som==s), 0, dens$bw)

      # Fix values outside input range (per marker) - per SOM node
      imputed[which(impute_obs_som==s)][imputed[which(impute_obs_som==s)] < min(complete_obs[which(complete_obs_som==s),channel])] <- min(complete_obs[which(complete_obs_som==s),channel])
      imputed[which(impute_obs_som==s)][imputed[which(impute_obs_som==s)] > max(complete_obs[which(complete_obs_som==s),channel])] <- max(complete_obs[which(complete_obs_som==s),channel])

    } else {

      imputed[which(impute_obs_som==s)] <- NA
      warning("Be aware that a cluster contains cells primarily from the dataset you wish to impute for. As a result, imputations were not made for those cells.")
    }
  }

  # Put the values back in the dataset
  df[df$batch %in% correct_batches, channel] <- imputed

  # Print statement warning users about directly interpreting imputed values
  message("Caution! Analysis based on imputed values can lead to false inferences. We recommend only using imputed marker expressions for visualization and not for differential expression analysis.")

  return(df)
}




#' Impute non-overlapping channels for whole data sets
#'
#' This function imputes the expression values for a whole dataset based on the overlapping markers contained in another dataset.
#'   The purpose is to be able to merge multi-panel data or to merge differently run datasets that have non-overlapping markers
#'
#' @param dataset1 Dataframe with expression values 1
#' @param dataset2 Dataframe with expression values 2
#' @param overlap_channels Channels (markers) that overlap between impute_for and complete_obs, which can be used to base imputation on
#' @param impute_channels1 Channels to impute and add to the impute_for dataset1 (must be present in dataset2). Alternatively can be set to NULL.
#' @param impute_channels2 Channels to impute and add to the impute_for dataset2 (must be present in dataset1). Alternatively can be set to NULL.
#' @inheritParams create_som
#' @family merging
#' @examples
#' \dontrun{
#' dfs_imputed <- impute_across_panels(dataset1 = df_panel1, dataset2 = df_panel2,
#'  overlap_channels = intersect(df_panel1_markers, df_panel2_markers),
#'  impute_channels1 = unique_df_panel2_markers, , impute_channels2 = unique_df_panel1_markers)
#' }
#' @export
impute_across_panels <- function(
    dataset1,
    dataset2,
    overlap_channels,
    impute_channels1,
    impute_channels2,
    label = NULL,
    xdim = 8,
    ydim = 8,
    seed = 482) {


  # Checking colnames
  if (!all(impute_channels1 %in% colnames(dataset2))) {
    stop("Error: Some of your impute_channels1 are not found among the dataset2 column names.")
  }
  if (!all(impute_channels2 %in% colnames(dataset1))) {
    stop("Error: Some of your impute_channels2 are not found among the dataset1 column names.")
  }
  if (!all(overlap_channels %in% colnames(dataset1))) {
    stop("Error: Some of your overlap_channels are not found among the dataset1 column names.")
  }
  if (!all(overlap_channels %in% colnames(dataset2))) {
    stop("Error: Some of your overlap_channels are not found among the dataset2 column names.")
  }

  # Checking if choice is to not impute for one dataset
  if (is.null(impute_channels1) && is.null(impute_channels2)) {
    stop("Error: You have not specified any channels to impute!")
  } else if (is.null(impute_channels1) || is.null(impute_channels2)) {
    message("Be aware that one of your panels will not have any imputed markers.")
  }


  # Get SOM classes for datasets on overlapping channels
  overlapping_data <- rbind(dataset1[,overlap_channels], dataset2[,overlap_channels])


  if (is(label, "NULL")){
    som_classes <- overlapping_data %>%
      cyCombine::create_som(markers = overlap_channels,
                            xdim = xdim,
                            ydim = ydim,
                            seed = seed)
  } else{
    som_classes <- overlapping_data[[label]]
  }


  # Split SOM classes to each original dataset
  dataset1_som <- som_classes[1:nrow(dataset1)]
  dataset2_som <- som_classes[(nrow(dataset1)+1):length(som_classes)]


  # Now, we impute the impute_channels
  imputed_dfs <- list(); other_set <- c(2,1)
  for (i in 1:2) {

    # Set the parameters for the loop to run
    impute_for <- eval(parse(text=paste0("dataset", i)))
    impute_channels <- eval(parse(text=paste0("impute_channels", i)))
    impute_obs_som <- eval(parse(text=paste0("dataset", i, "_som")))


    complete_obs <- eval(parse(text=paste0("dataset", other_set[i])))
    complete_obs_som <- eval(parse(text=paste0("dataset", other_set[i], "_som")))


    # Check if impute_channels is NULL, in that case - skip imputation for the panel
    if (is.null(impute_channels)) {
      message(paste0("Skipping imputation for dataset", i, "."))

      impute_for <- impute_for %>%
        dplyr::relocate(dplyr::any_of(non_markers), .after = dplyr::last_col())

      imputed_dfs[[paste0("dataset", i)]] <- impute_for
      next
    }


    # For each missing event in each cluster, impute values based on density draws in the same cluster
    message(paste0("Performing density draws for dataset", i, "."))
    imputed <- matrix(nrow=nrow(impute_for), ncol=length(impute_channels), dimnames = list("rn"=NULL, "cn"=impute_channels))

    for (s in sort(unique(impute_obs_som))) {

      # Check if there's at least 50 cells from the 'training' data in the given cluster and if so, impute - otherwise set values to NA and throw a warning
      if (sum(complete_obs_som==s) >= 50) {

        # Estimate the density per marker that needs imputation
        dens <- apply(complete_obs[which(complete_obs_som==s), impute_channels], 2, function(x) {stats::density(x)$bw})

        # Performing the imputation
        set.seed(seed)
        imputed[which(impute_obs_som==s), ] <- complete_obs[which(complete_obs_som==s), ] %>%
          dplyr::select(dplyr::all_of(impute_channels)) %>%
          dplyr::slice_sample(n = sum(impute_obs_som==s), replace = T) %>%
          as.matrix() + sapply(impute_channels, function(ch) {stats::rnorm(sum(impute_obs_som==s), 0, dens[ch])})

        # Fix values outside input range (per marker)
        for (m in impute_channels) {
          imputed[which(impute_obs_som==s), m][imputed[which(impute_obs_som==s), m] < min(complete_obs[which(complete_obs_som==s), m])] <- min(complete_obs[which(complete_obs_som==s), m])
          imputed[which(impute_obs_som==s), m][imputed[which(impute_obs_som==s), m] > max(complete_obs[which(complete_obs_som==s), m])] <- max(complete_obs[which(complete_obs_som==s), m])
        }

      } else {

        imputed[which(impute_obs_som==s),impute_channels] <- NA
        warning("Be aware that a cluster contains cells primarily/only from the dataset you wish to impute for. As a result, imputations were not made for those cells.\n")
      }
    }



    # Add new columns to impute_for
    impute_for[,impute_channels] <- imputed

    impute_for <- impute_for %>%
      dplyr::relocate(dplyr::any_of(non_markers), .after = dplyr::all_of(impute_channels))

    imputed_dfs[[paste0("dataset", i)]] <- impute_for
  }

  # Print statement warning users about directly interpreting imputed values
  message("Caution! Analysis based on imputed values can lead to false inferences. We recommend only using imputed marker expressions for visualization and not for differential expression analysis.")

  return(imputed_dfs)
}

