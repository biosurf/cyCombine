#### Plotting functions ----

#' Density ridges for two sets
#' @import
density_plots <- function(uncorrected, corrected, batch_ids, uncorrected_order = as.factor(batch_ids), corrected_order = as.factor(batch_ids), filename) {

  # Extract data into dataframe format
  uncor_df <- cbind.data.frame(uncorrected, uncorrected_order, rep('Uncorrected', nrow(uncorrected)))
  cor_df <- cbind.data.frame(corrected, corrected_order, rep('Corrected', nrow(corrected)))
  colnames(cor_df)[(ncol(cor_df)-1):ncol(cor_df)] <- colnames(uncor_df)[(ncol(uncor_df)-1):ncol(uncor_df)] <- c('Batch', 'Type')

  df <- rbind.data.frame(uncor_df, cor_df)

  # For each marker, make the plot
  p <- list()
  for (c in 1:length(all_markers)) {

    p[[c]] <- ggplot(df, aes_string(x = all_markers[c], y = "Batch")) +
      geom_density_ridges(aes(color = Type, fill = Type), alpha = 0.4) +
      theme_bw()
  }

  # Save the plots
  save_plot(filename, plot_grid(plotlist = p, ncol = 6), base_width = 28, base_height = 40)

}



# Dimensionality reduction plot
dimred_plot <- function(data, batch_ids, name, type = 'pca', plot = 'batch', marker = NULL) {

  if (type == 'pca') {
    # Run PCA
    pca <- data %>%
      select_if(is.numeric) %>%
      prcomp(scale. = TRUE, center = TRUE)

    # Make dataframe with output
    if (plot == 'batch') {
      df <- pca$x %>%
        as_tibble() %>%
        mutate(Batch = as.factor(batch_ids))
      #cbind.data.frame(pca$x, as.factor(batch_ids)); colnames(df)[ncol(df)] <- 'Batch'
    } else {
      df <- cbind.data.frame(pca$x, as.factor(batch_ids), data[,marker]); colnames(df)[(ncol(df)-1):ncol(df)] <- c('Batch', marker)
    }

  } else if (type == 'umap') {
    # Run UMAP
    set.seed(758)
    umap <- data %>%
      select_if(is.numeric) %>%
      umap(n_neighbors = 15, min_dist = 0.2, metric = 'euclidean')

    # Make dataframe with output
    if (plot == 'batch') {
      df <- umap %>%
        as_tibble() %>%
        mutate(Batch = as.factor(batch_ids)) %>%
        rename(UMAP1 = V1,
               UMAP2 = V2)
      # cbind.data.frame(umap, as.factor(batch_ids)); colnames(df) <- c('UMAP1', 'UMAP2', 'Batch')
    } else {
      df <- cbind.data.frame(umap, as.factor(batch_ids), data[,marker]); colnames(df) <- c('UMAP1', 'UMAP2', 'Batch', marker)
    }
  }

  # Make the plot
  if (plot == 'batch') {
    plot <- df %>%
      ggplot(aes_string(x = colnames(df)[1], y = colnames(df)[2])) +
      geom_point(aes_string(color = "Batch"), alpha = 0.3, size = 0.4) +
      guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) +
      theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
      ggtitle(paste(toupper(type), '-', name))
  } else {
    plot <- ggplot(df, aes_string(x = colnames(df)[1], y = colnames(df)[2])) +
      geom_point(aes_string(color = marker), alpha = 0.3, size = 0.4) +
      #guides(colour = guide_legend(override.aes = list(alpha = 1, size = 1))) +
      theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
      ggtitle(paste(toupper(type), '-', name))
  }

  return(plot)
}



# Save two plots aligned with cowplot
save_two_plots <- function(plot1, plot2, filename) {
  plot <- plot_grid(plot1, plot2, align = 'v', scale = 0.9)
  save_plot(filename = filename, plot, base_width = 12, base_height = 6)
}

# Save four plots aligned with cowplot
save_four_plots <- function(plot1, plot2, plot3, plot4, filename) {
  plot <- plot_grid(plot1, plot2, plot3, plot4, align = 'v', scale = 0.9, nrow = 2)
  save_plot(filename = filename, plot, base_width = 12, base_height = 12)
}
