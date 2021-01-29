#### Plotting functions ----

# @import ggplot2
# @import ggridges
# @import cowplot

# @importFrom dplyr select_if bind_rows rename
#' Density ridges for two sets
#'
#' Compare densities between each batch both before and after correction
#' @import ggplot2
#' @param format Plotting format (1 = 1 row per batch, 2 = all batches in same row.)
#' @param uncorrected tibble with uncorrected data
#' @param corrected tibble with corrected data
#' @param markers Character vector with the markers to plot
#' @param filename Output figure filename. If NULL plots are returned
#' @param y The column to stack densities with. Default: "batch". If set to "Type", the dataset_names will be used to stack densities.
#' @param xlim The limit of the x axis
#' @param dataset_names Change the names of the datasets from Uncorrected and Corrected to something else. Format: c("Uncorrected", "Corrected")
#' @param ncol Number of density plots in a single row of plots. Default: 6
#' @examples
#' plot_density(uncorrected, corrected, y = 'batch', filename = 'my/dir/batchcor_plot.pdf')
#' plot_density(imputed1, imputed2, y = 'Type', dataset_names = paste('Panel', 1:2), filename = 'my/dir/merging_plot.pdf')
#' @export
plot_density <- function(uncorrected,
                         corrected,
                         markers = NULL,
                         filename = NULL,
                         y = "batch",
                         xlim = 10,
                         dataset_names = NULL,
                         ncol = 6,
                         format = 1) {

  # Check for packages
  missing_package("ggridges", "CRAN")
  missing_package("ggplot2", "CRAN")
  missing_package("cowplot", "CRAN")

  # Get markers
  if (is.null(markers)) {
    markers <- uncorrected %>%
      get_markers()
  }

  # Rename datasets
  if (is.null(dataset_names)) {
    dataset_names <- c('Uncorrected', 'Corrected')
  } else if (length(dataset_names) != 2) {
    dataset_names <- c('Dataset 1', 'Dataset 2')
  }

  # Define how densities are stacked
  if (y == 'Type') {
    batch1 <- dataset_names[1]
    batch2 <- dataset_names[2]
  } else {
    batch1 <- as.factor(uncorrected[[y]])
    batch2 <- as.factor(corrected[[y]])
  }

  # Isolate markers and add Type and batch columns
  uncorrected <- uncorrected %>%
    dplyr::select(all_of(markers)) %>%
    dplyr::mutate(Type = dataset_names[1],
                  batch = batch1)

  # Combine uncorrected and corrected data
  df <- corrected %>%
    dplyr::select(all_of(markers)) %>%
    dplyr::mutate(Type = dataset_names[2],
                  batch = batch2) %>%
    dplyr::bind_rows(uncorrected)

  # 1 row per batch
  if (format == 1) {
    # For each marker, make the plot
    p <- list()
    for (c in 1:length(markers)) {

      p[[c]] <- df %>%
        ggplot(aes_string(x = markers[c], y = 'batch')) +
        ggridges::geom_density_ridges(aes(color = .data$Type, fill = .data$Type), alpha = 0.4) +
        coord_cartesian(xlim = c(-1, xlim)) +
        labs(y = y) +
        theme_bw()
    }

    height_factor <- 1.3

  } else if (format == 2) { # all batches in same row
    # For each marker, make the plot
    p <- list()
    for (c in 1:length(markers)) {

      p[[c]] <- df %>%
        ggplot(aes_string(x = markers[c], y = 'Type')) +
        ggridges::geom_density_ridges(aes(color = batch, fill = batch), alpha = 0.4) +
        coord_cartesian(xlim = c(-1, xlim)) +
        labs(y = y) +
        theme_bw()
    }

    height_factor <- 2
  }


  if (!is.null(filename)) {
    cowplot::save_plot(filename, cowplot::plot_grid(plotlist = p, ncol = ncol), base_width = 28, base_height = length(markers)/height_factor)
  } else {
    return(cowplot::plot_grid(plotlist = p, ncol = ncol))
  }
}



# @importFrom uwot umap
#' Dimensionality reduction plot
#'
#' @param df tibble of data to plot
#' @param name Name of the output plot
#' @param type Type of dimensionality reduction ("umap" or "pca")
#' @param plot The column to use for coloring
#' @param markers Markers to include in dimensionality reduction
#' @param seed For reproducibility
#' @param return_coord Return coordinates and not just the plot
#' @examples
#' uncor_umap <- plot_dimred(uncorrected, "Uncorrected", markers = markers)
#' @export
plot_dimred <- function(df,
                        name,
                        type = "umap",
                        plot = "batch",
                        markers = NULL,
                        seed = 473,
                        return_coord = FALSE) {

  # Check for missing packages
  if(type == "umap") missing_package("uwot", "CRAN")
  missing_package("ggridges", "CRAN")
  if(plot != "batch") missing_package("viridis", "CRAN")

  if (!(type %in% c('pca', 'umap'))) {
    stop("Error, please use either type = 'pca' or type = 'umap'.")
  }


  if(is.null(markers)){
    markers <- cyCombine::get_markers(df)
  }

  Batch <- df$batch %>%
    as.factor()
  df <- df %>%
    dplyr::select(dplyr::all_of(markers))

  if (type == "pca") {
    # Run PCA
    pca <- df %>%
      stats::prcomp(scale. = TRUE, center = TRUE)

    # Make dataframe with output
    if (plot == "batch") {
      df <- pca$x %>%
        tibble::as_tibble() %>%
        dplyr::mutate(Batch = as.factor(Batch))

    } else {
      df <- cbind.data.frame(pca$x, as.factor(Batch), df[, plot]); colnames(df)[(ncol(df)-1):ncol(df)] <- c("Batch", plot)
    }

  } else if (type == "umap") {
    # Run UMAP
    set.seed(seed)
    umap <- df %>%
      uwot::umap(n_neighbors = 15, min_dist = 0.2, metric = "euclidean")

    # Make dataframe with output
    if (plot == "batch") {
      colnames(umap) <- c("UMAP1", "UMAP2")
      df <- umap %>%
        tibble::as_tibble() %>%
        dplyr::mutate(Batch = Batch)

    } else {
      df <- cbind.data.frame(umap, Batch, df[,plot]); colnames(df) <- c("UMAP1", "UMAP2", "Batch", plot)
    }
  }

  # Make the plot
  if (plot == "batch") {
    plot <- df %>%
      ggplot(aes_string(x = colnames(df)[1],
                        y = colnames(df)[2])) +
      geom_point(aes_string(color = "Batch"),
                 alpha = 0.3,
                 size = 0.4,
                 shape = 1) +
      guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) +
      theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
      ggtitle(paste(toupper(type), "-", name))
  } else {
    plot <- df %>%
      ggplot(aes_string(x = colnames(df)[1],
                        y = colnames(df)[2])) +
      geom_point(aes_string(color = plot),
                 alpha = 0.3,
                 size = 0.4) +
      #guides(colour = guide_legend(override.aes = list(alpha = 1, size = 1))) +
      theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
      ggtitle(paste(toupper(type), "-", name)) +
      scale_color_viridis_c()
  }


  if (return_coord) {
    if (type == 'pca') {
      return(list("plot" = plot, "dimred" = pca$x))
    } else {
      colnames(umap) <- c('UMAP1', 'UMAP2')
      return(list("plot" = plot, "dimred" = umap))
    }
  } else {
    return(plot)
  }

}


# @importFrom uwot umap
#' Dimensionality reduction plots - colored with labels, batches and marker expression
#' @inheritParams plot_dimred
#' @param out_dir Directory to put output figures
#' @export
plot_dimred_full <- function(df,
                             name,
                             type = "umap",
                             markers = NULL,
                             seed = 473,
                             out_dir = NULL) {

  if(type == "umap") missing_package("uwot", "CRAN")
  missing_package("ggridges", "CRAN")

  # Check out dir
  if (is.null(out_dir)) {
    stop('Error! Please speicify output directory.')
  }

  if(is.null(markers)){
    markers <- cyCombine::get_markers(df)
  }
  Batch <- df$batch %>%
    as.factor()

  Label <- df$label %>%
    as.factor()

  df <- df %>%
    dplyr::select(dplyr::all_of(markers))

  if (type == "pca") {
    # Run PCA
    pca <- df %>%
      prcomp(scale. = TRUE, center = TRUE)

    # Make dataframe with output
    df <- cbind.data.frame(pca$x, as.factor(Batch), as.factor(Label), df); colnames(df)[3:ncol(df)] <- c("Batch", "Label", markers)

  } else if (type == "umap") {
    # Run UMAP
    set.seed(seed)
    umap <- df %>%
      uwot::umap(n_neighbors = 15, min_dist = 0.2, metric = "euclidean")

    # Make dataframe with output
    df <- cbind.data.frame(umap, as.factor(Batch), as.factor(Label), df); colnames(df) <- c("UMAP1", "UMAP2", "Batch", "Label", markers)
  }

  # Make the plots
  batch_plot <- df %>%
    ggplot(aes_string(x = colnames(df)[1],
                      y = colnames(df)[2])) +
    geom_point(aes_string(color = "Batch"),
               alpha = 0.3,
               size = 0.4,
               shape = 1) +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle(paste(toupper(type), "-", name))

  label_plot <- df %>%
    ggplot(aes_string(x = colnames(df)[1],
                      y = colnames(df)[2])) +
    geom_point(aes_string(color = "Label"),
               alpha = 0.3,
               size = 0.4,
               shape = 1) +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle(paste(toupper(type), "-", name))

  # Saving plots
  cyCombine::plot_save_two(batch_plot, label_plot, paste0(out_dir, '/UMAP_batches_labels.png'))

  # Marker plots
  cyCombine::check_make_dir(paste0(out_dir, '/UMAP_markers'))

  marker_plots <- list()
  for (m in markers) {
    p <- ggplot(df, aes_string(x = colnames(df)[1],
                               y = colnames(df)[2])) +
      geom_point(aes_string(color = m),
                 alpha = 0.3,
                 size = 0.4) +
      scale_color_gradientn(m, colors = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))(50)) +
      theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
      ggtitle(paste(toupper(type), "-", name))

    suppressMessages(ggsave(p, filename = paste0(out_dir, '/UMAP_markers/UMAP_batches_labels', m, '.png')))
  }

}



# @importFrom uwot umap
#' Dimensionality reduction plots for a two-batch dataset before and after correction - colored by batch and label
#' @inheritParams plot_dimred
#' @export
plot_umap_labels <- function(uncorrected,
                             corrected,
                             filename = NULL,
                             seed = 473) {


  missing_package("uwot", "CRAN")
  missing_package("ggridges", "CRAN")
  missing_package("cowplot", "CRAN")


  # Check if two-batch data and extract batch names
  if (length(unique(c(unique(uncorrected$batch), unique(corrected$batch)))) == 2) {
    batch1 <- unique(uncorrected$batch)[1]
    batch2 <- unique(uncorrected$batch)[2]

  } else {
    stop('Error, data must contain only two batches.')

  }


  # Make a list format with all the plots to be made
  datasets <- list(list(df = uncorrected[uncorrected$batch == batch1,], plot = 'Label', nplots = 1, name = paste0('Uncorrected ', batch1)),
                   list(df = uncorrected[uncorrected$batch == batch2,], plot = 'Label', nplots = 1, name = paste0('Uncorrected ', batch2)),
                   list(df = uncorrected, plot = 'Batch', nplots = 3, name = paste0('Uncorrected ', batch1, ' + ', batch2)),
                   list(df = corrected, plot = 'Batch', nplots = 3, name = paste0('Corrected ', batch1, ' + ', batch2)),
                   list(df = corrected[corrected$batch == batch1,], plot = 'Label', nplots = 1, name = paste0('Corrected ', batch1)),
                   list(df = corrected[corrected$batch == batch2,], plot = 'Label', nplots = 1, name = paste0('Corrected ', batch2)))

  seed = seed
  plots <- plots2 <- list()
  plot_count = plot_count2 = 1
  markers = get_markers(uncorrected)
  for (i in 1:6) {

    df <- datasets[[i]][['df']] #%>% sample_n(1000)
    plot <- datasets[[i]][['plot']]
    nplots <- datasets[[i]][['nplots']]
    name <- datasets[[i]][['name']]

    Batch <- df$batch %>%
      as.factor()
    Label <- df$label %>%
      factor(levels = c(unique(corrected$label)))
    df <- df %>%
      dplyr::select(dplyr::all_of(markers))


    # Run UMAP
    set.seed(seed)
    umap <- df %>%
      uwot::umap(n_neighbors = 15, min_dist = 0.2, metric = "euclidean")


    colnames(umap) <- c("UMAP1", "UMAP2")
    df <- umap %>%
      tibble::as_tibble() %>%
      dplyr::mutate(Batch = Batch,
                    Label = Label)

    # Make the plot
    if (nplots == 1) {
      plot <- df %>%
        ggplot(aes_string(x = colnames(df)[1],
                          y = colnames(df)[2])) +
        geom_point(aes_string(color = plot),
                   alpha = 0.3,
                   size = 0.4,
                   shape = 1) +
        labs(plot) +
        guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) +
        theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
        ggtitle(paste("UMAP -", name))

      plots[[plot_count]] <- plot
      plot_count <- plot_count + 1

    } else if (nplots == 3) {

      # Together plot
      plot <- df %>%
        ggplot(aes_string(x = colnames(df)[1],
                          y = colnames(df)[2])) +
        geom_point(aes_string(color = plot),
                   alpha = 0.3,
                   size = 0.4,
                   shape = 1) +
        guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) +
        labs(plot) +
        theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
        ggtitle(paste("UMAP -", name))

      plots[[plot_count]] <- plot
      plot_count <- plot_count + 1

      # Faceted plots
      plot <- df %>%
        ggplot(aes_string(x = colnames(df)[1],
                          y = colnames(df)[2])) +
        facet_wrap(~Batch) +
        geom_point(aes_string(color = Label),
                   alpha = 0.3,
                   size = 0.4,
                   shape = 1) +
        guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) +
        labs(colour = 'Label') +
        theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
        ggtitle(paste("UMAP -", name))

      plots2[[plot_count2]] <- plot
      plot_count2 = plot_count2 + 1

    }
  }

  # Save the plots
  plot1 <- cowplot::plot_grid(plotlist = plots, scale = 0.9)
  plot2 <- cowplot::plot_grid(plotlist = plots2, scale = 0.9)

  if (!is.null(filename)) {
    cowplot::save_plot(filename = filename, plot1, base_width = 24, base_height = 14)
    cowplot::save_plot(filename = paste0(tools::file_path_sans_ext(filename), '_facets.', tools::file_ext(filename)), plot2, base_width = 32, base_height = 8)

  } else {
    return(list(plot1, plot2))

  }

}




#' Save two plots aligned with cowplot
#' @param plot1 Left plot in output
#' @param plot2 Right plot in output
#' @param filename Filename of output
#' @export
plot_save_two <- function(plot1, plot2, filename) {
  missing_package("cowplot", "CRAN")

  plot <- cowplot::plot_grid(plot1, plot2, align = "v", scale = 0.9)
  cowplot::save_plot(filename = filename, plot, base_width = 12, base_height = 6)
}


#' Save four plots aligned with cowplot
#' @param plot1 Top left plot in output
#' @param plot2 Top right plot in output
#' @param plot3 Lower left plot in output
#' @param plot4 Lower right plot in output
#' @param filename Filename of output
#' @export
plot_save_four <- function(plot1, plot2, plot3, plot4, filename) {
  missing_package("cowplot", "CRAN")

  plot <- cowplot::plot_grid(plot1, plot2, plot3, plot4, align = "v", scale = 0.9, nrow = 2)
  cowplot::save_plot(filename = filename, plot, base_width = 12, base_height = 12)
}
