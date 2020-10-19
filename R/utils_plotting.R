#### Plotting functions ----

# @import ggplot2
# @import ggridges
# @import cowplot

# @importFrom dplyr select_if bind_rows rename
#' Density ridges for two sets
#' @import ggplot2
#' @export
plot_density <- function(uncorrected, corrected, markers = NULL, filename) {

  # Check for packages
  missing_package("ggridges", "CRAN")
  missing_package("ggplot2", "CRAN")
  missing_package("cowplot", "CRAN")



  if(is.null(markers)){
    markers <- uncorrected %>%
      get_markers()
  }

  uncorrected <- uncorrected %>%
    dplyr::select(all_of(markers)) %>%
    dplyr::mutate(Type = "Uncorrected",
                  batch = as.factor(uncorrected$batch))

  df <- corrected %>%
    dplyr::select(all_of(markers)) %>%
    dplyr::mutate(Type = "Corrected",
                  batch = as.factor(corrected$batch)) %>%
    dplyr::bind_rows(uncorrected)


  # Extract data into dataframe format
  # uncor_df <- cbind.data.frame(uncorrected, uncorrected_order, rep("Uncorrected", nrow(uncorrected)))
  # cor_df <- cbind.data.frame(corrected, corrected_order, rep("Corrected", nrow(corrected)))
  # colnames(df)[(ncol(df)-1):ncol(df)] <- c("Batch", "Type")
  #
  # df <- rbind.data.frame(uncor_df, cor_df)

  # For each marker, make the plot
  p <- list()
  for (c in 1:length(markers)) {

    p[[c]] <- df %>%
      ggplot(aes_string(x = markers[c], y = "batch")) +
      ggridges::geom_density_ridges(aes(color = .data$Type, fill = .data$Type), alpha = 0.4) +
      theme_bw()
  }

  # Save the plots
  # ggsave(filename = filename, plot = plot_grid(plotlist = p, ncol = 6),
         # device = "png", width = 28, height = 40)
  cowplot::save_plot(filename, cowplot::plot_grid(plotlist = p, ncol = 6), base_width = 28, base_height = 40)

}



# @importFrom uwot umap
#' Dimensionality reduction plot
#' @export
plot_dimred <- function(df, name, type = "pca", plot = "batch", marker = NULL) {

  missing_package("uwot", "CRAN")
  missing_package("ggplot2", "CRAN")
  missing_package("ggridges", "CRAN")


  Batch <- df$batch %>%
    as.factor()
  df <- df %>%
    dplyr::select_if(names(.) %!in% c("batch", "sample", "covar", "id"))

  if (type == "pca2") {
    # Run PCA
    pca <- df %>%
      prcomp(scale. = TRUE, center = TRUE)

    # Make dataframe with output
    if (plot == "batch") {
      df <- pca$x %>%
        tibble::as_tibble() %>%
        dplyr::mutate(Batch = as.factor(Batch))
      #cbind.data.frame(pca$x, as.factor(batch_ids)); colnames(df)[ncol(df)] <- "Batch"
    } else {
      df <- cbind.data.frame(pca$x, as.factor(Batch), data[,marker]); colnames(df)[(ncol(df)-1):ncol(df)] <- c("Batch", marker)
    }

  } else if (type == "umap") {
    # Run UMAP
    set.seed(758)
    umap <- df %>%
      uwot::umap(n_neighbors = 15, min_dist = 0.2, metric = "euclidean")

    # Make dataframe with output
    if (plot == "batch") {
      colnames(umap) <- c("UMAP1", "UMAP2")
      df <- umap %>%
        tibble::as_tibble() %>%
        dplyr::mutate(Batch = Batch) #%>%
        # dplyr::rename(UMAP1 = V1,
                      # UMAP2 = V2)
      # cbind.data.frame(umap, as.factor(batch_ids)); colnames(df) <- c("UMAP1", "UMAP2", "Batch")
    } else {
      df <- cbind.data.frame(umap, Batch, data[,marker]); colnames(df) <- c("UMAP1", "UMAP2", "Batch", marker)
    }
  }

  # Make the plot
  if (plot == "batch") {
    plot <- df %>%
      ggplot(aes_string(x = colnames(df)[1], y = colnames(df)[2])) +
      geom_point(aes_string(color = "Batch"), alpha = 0.3, size = 0.4, shape = 1) +
      guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) +
      theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
      ggtitle(paste(toupper(type), "-", name))
  } else {
    plot <- ggplot(df, aes_string(x = colnames(df)[1], y = colnames(df)[2])) +
      geom_point(aes_string(color = marker), alpha = 0.3, size = 0.4) +
      #guides(colour = guide_legend(override.aes = list(alpha = 1, size = 1))) +
      theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
      ggtitle(paste(toupper(type), "-", name))
  }

  return(plot)
}




#' Save two plots aligned with cowplot
#' @export
plot_save_two <- function(plot1, plot2, filename) {
  missing_package("cowplot", "CRAN")

  plot <- cowplot::plot_grid(plot1, plot2, align = "v", scale = 0.9)
  cowplot::save_plot(filename = filename, plot, base_width = 12, base_height = 6)
}


#' Save four plots aligned with cowplot
#' @export
plot_save_four <- function(plot1, plot2, plot3, plot4, filename) {
  missing_package("cowplot", "CRAN")

  plot <- cowplot::plot_grid(plot1, plot2, plot3, plot4, align = "v", scale = 0.9, nrow = 2)
  cowplot::save_plot(filename = filename, plot, base_width = 12, base_height = 12)
}
