#### Plotting functions ----


#' Density ridges for two sets
#'
#' Compare densities between each batch both before and after correction.
#' Thanks to GitHub user asongggg for suggesting storing markers on separate pages and for suggesting an approach.
#' @param format Plotting format (1 = 1 row per batch, 2 = all batches in same row.)
#' @param uncorrected tibble with uncorrected data
#' @param corrected tibble with corrected data
#' @param markers Character vector with the markers to plot
#' @param filename Output figure filename. If NULL plots are returned
#' @param y The column to stack densities with. Default: "batch". If set to "Type", the dataset_names will be used to stack densities.
#' @param xlims The limits of the x axis as a vector
#' @param dataset_names Change the names of the datasets from Uncorrected and Corrected to something else. Format: c("Uncorrected", "Corrected")
#' @param ncol Number of density plots in a single row of plots. Default: 6
#' @param directory Folder to store figures in.
#' @param markers_per_page Store figures on individual pages, with the number of markers as specified. Default: NULL.
#' @param return_format Allows returning the raw plots for manual adjustments. This argument is mainly used when plots are not stored in files.
#' @family plot
#' @examples
#' \dontrun{
#' plot_density(uncorrected, corrected, y = 'batch', filename = 'my/dir/batchcor_plot.pdf')
#' plot_density(imputed1, imputed2, y = 'Type', dataset_names = paste('Panel', 1:2),
#'              filename = 'my/dir/merging_plot.pdf')
#' }
#' @export
plot_density <- function(uncorrected, corrected, markers = NULL, directory = NULL,
                         filename = NULL, y = "batch", dataset_names = NULL,
                         xlims = c(-1, 10),
                         ncol = 6, format = 1, markers_per_page = NULL,
                         return_format = c("cowplot", "raw", "pdf")) {
  return_format <- match.arg(return_format)
  # Check for packages
  cyCombine:::missing_package("ggridges", "CRAN")
  cyCombine:::missing_package("ggplot2", "CRAN")
  if(return_format != "raw") cyCombine:::missing_package("cowplot", "CRAN")

  if(return_format == "pdf" & is.null(filename)){
    stop("Please provide a filename when returning as pdf.")
  }
  # Ensure ´markers_per_page´ is NULL when returning raw.
  if(return_format == "raw" | is.null(filename)) markers_per_page <- NULL

  # Get markers
  if (is.null(markers)) {
    markers <- cyCombine::get_markers(uncorrected)
  }

  # Combine directory and filename in filename
  if (!is.null(filename)) {
    filename <- ifelse(is.null(directory), filename, file.path(directory, filename))
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

  # Combine uncorrected and corrected data
  df <- corrected %>%
    dplyr::select(dplyr::all_of(markers)) %>%
    dplyr::mutate(Type = dataset_names[2], batch = batch2) %>%
    dplyr::bind_rows(uncorrected %>%
                       dplyr::select(dplyr::all_of(markers)) %>%
                       dplyr::mutate(Type = dataset_names[1], batch = batch1))

  # Determine number of plots and pages
  num_plots <- length(markers)
  num_pages <- ifelse(is.null(markers_per_page), 1, ceiling(num_plots / markers_per_page))

  # Create a list to store the plots
  plot_list <- list()

  for (page in 1:num_pages) {
    start_marker <- ifelse(is.null(markers_per_page), 1, (page - 1) * markers_per_page + 1)
    end_marker <- ifelse(is.null(markers_per_page), num_plots, min(page * markers_per_page, num_plots))

    for (c in start_marker:end_marker) {
      if (format == 1) {
        p <- df %>%
          ggplot2::ggplot(ggplot2::aes_string(x = paste0('`', markers[c], '`'), y = 'batch')) +
          ggridges::geom_density_ridges(ggplot2::aes(color = .data$Type, fill = .data$Type), alpha = 0.4) +
          ggplot2::coord_cartesian(xlim = xlims) +
          ggplot2::labs(y = y) +
          ggplot2::theme_bw()
      } else if (format == 2) {
        p <- df %>%
          ggplot2::ggplot(ggplot2::aes_string(x = paste0('`', markers[c], '`'), y = 'Type')) +
          ggridges::geom_density_ridges(ggplot2::aes(color = batch, fill = batch), alpha = 0.4) +
          ggplot2::coord_cartesian(xlim = xlims) +
          ggplot2::labs(y = y) +
          ggplot2::theme_bw()
      }
      plot_list[[c]] <- p
    }
    # Store or return plots
    if (return_format == "raw"){
      return(plot_list)
    }

    # Store shared legend from first plot
    if (page == 1){
      legend <- cowplot::get_legend(
        # create some space to the left of the legend
        plot_list[[1]] + ggplot2::theme(legend.box.margin = ggplot2::margin(0, 0, 0, 12))
      )
    }

    # Remove shared legend from plots
    for (i in 1:length(plot_list)) {
      plot_list[[i]] <- plot_list[[i]] + ggplot2::theme(legend.position="none")
    }

    if (is.null(filename)){
      return(cowplot::plot_grid(cowplot::plot_grid(plotlist = plot_list, ncol = ncol), legend, rel_widths = c(2*ncol,1)))
    } else{
      current_markers <- markers[start_marker:end_marker]
      marker_names <- paste(current_markers, collapse = "_")

      device <- tools::file_ext(filename)
      if (device == "") device <- "pdf"
      filename <- gsub(paste0(".", device), "", filename)
      file <- ifelse(num_pages > 1,
                     yes = file.path(sprintf("%s_page%d_%s.%s", filename, page, marker_names, device)),
                     no = paste0(filename,".", device))
      height_factor <- ifelse(format == 1, 1.3, 2)
      cowplot::save_plot(file, cowplot::plot_grid(cowplot::plot_grid(plotlist = plot_list[start_marker:end_marker], ncol = ncol), legend, rel_widths = c(2*ncol,1)), base_width = 20, base_height = length(markers)/height_factor)
    }
  }
}






#' Dimensionality reduction plot
#'
#' @param df tibble of data to plot
#' @param name Name of the output plot
#' @param type Type of dimensionality reduction ("umap" or "pca")
#' @param plot The column to use for coloring
#' @param markers Markers to include in dimensionality reduction
#' @param seed For reproducibility
#' @param return_coord Return coordinates and not just the plot
#' @family plot
#' @examples
#' \dontrun{
#' uncor_umap <- plot_dimred(uncorrected, "Uncorrected", markers = markers)
#' }
#' @export
plot_dimred <- function(df,
                        name,
                        type = c("umap", "pca"),
                        plot = "batch",
                        markers = NULL,
                        seed = 473,
                        return_coord = FALSE,
                        metric = "euclidean") {

  # Check for missing packages
  if (type == "umap") check_package("uwot", "CRAN")
  if (plot != "batch") check_package("viridis", "CRAN")
  check_package("ggridges", "CRAN")
  check_package("ggplot2", "CRAN")

  type <- match.arg(type)


  if(is.null(markers)){
    markers <- cyCombine::get_markers(df)
  }

  Batch <- df$batch %>%
    as.factor()
  # df <- df %>%
    # dplyr::select(dplyr::all_of(markers))

  if (type == "pca") {
    # Run PCA
    pca <- df %>%
      dplyr::select(dplyr::all_of(markers)) %>%
      stats::prcomp(scale. = TRUE, center = TRUE)

    # Make dataframe with output
    if (plot == "batch") {
      df <- pca$x %>%
        tibble::as_tibble() %>%
        dplyr::mutate(Batch = as.factor(Batch))

    } else {
      df <- cbind.data.frame(pca$x, as.factor(Batch), df[, plot])
      colnames(df)[(ncol(df)-1):ncol(df)] <- c("Batch", plot)
    }

  } else if (type == "umap") {
    # Run UMAP
    set.seed(seed)
    umap <- df %>%
      dplyr::select(dplyr::all_of(markers)) %>%
      uwot::umap(n_neighbors = 15, min_dist = 0.2, metric = metric)

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
      ggplot2::ggplot(ggplot2::aes_string(x = colnames(df)[1],
                                          y = colnames(df)[2])) +
      ggplot2::geom_point(ggplot2::aes_string(color = "Batch"),
                          alpha = 0.3,
                          size = 0.4,
                          shape = 1) +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(alpha = 1, size = 1))) +
      ggplot2::theme_bw() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::ggtitle(paste(toupper(type), "-", name))
  } else {
    plot <- df %>%
      ggplot2::ggplot(ggplot2::aes_string(x = colnames(df)[1],
                                          y = colnames(df)[2])) +
      ggplot2::geom_point(ggplot2::aes_string(color = plot),
                          alpha = 0.3,
                          size = 0.4) +
      #guides(colour = guide_legend(override.aes = list(alpha = 1, size = 1))) +
      ggplot2::theme_bw() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::ggtitle(paste(toupper(type), "-", name)) +
      ggplot2::scale_color_viridis_c()
  }


  if (return_coord) {
    return(list("plot" = plot, "dimred" = df))
  } else {
    return(plot)
  }

}


#' Dimensionality reduction plots - colored with labels, batches and marker expression
#' @inheritParams plot_dimred
#' @param out_dir Directory to put output figures
#'
#' @family plot
#' @export
plot_dimred_full <- function(df,
                             name,
                             type = "umap",
                             markers = NULL,
                             seed = 473,
                             out_dir = NULL) {

  if(type == "umap") cyCombine:::missing_package("uwot", "CRAN")
  check_package("ggridges", "CRAN")
  check_package("grDevices", "CRAN")
  check_package("RColorBrewer", "CRAN")
  check_package("ggplot2", "CRAN")

  # Check out dir
  if (is.null(out_dir)) {
    stop('Error! Please speicify output directory.')
  } else{
    # Create output directory if missing
    cyCombine:::check_make_dir(out_dir)
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
      stats::prcomp(scale. = TRUE, center = TRUE)

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
    ggplot2::ggplot(ggplot2::aes_string(x = colnames(df)[1],
                                        y = colnames(df)[2])) +
    ggplot2::geom_point(ggplot2::aes_string(color = "Batch"),
                        alpha = 0.3,
                        size = 0.4,
                        shape = 1) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(alpha = 1, size = 1))) +
    ggplot2::theme_bw() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::ggtitle(paste(toupper(type), "-", name))

  label_plot <- df %>%
    ggplot2::ggplot(ggplot2::aes_string(x = colnames(df)[1],
                                        y = colnames(df)[2])) +
    ggplot2::geom_point(ggplot2::aes_string(color = "Label"),
                        alpha = 0.3,
                        size = 0.4,
                        shape = 1) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(alpha = 1, size = 1))) +
    ggplot2::theme_bw() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::ggtitle(paste(toupper(type), "-", name))

  # Saving plots
  cyCombine::plot_save_two(batch_plot, label_plot, paste0(out_dir, '/UMAP_batches_labels.png'))

  # Marker plots
  cyCombine:::check_make_dir(paste0(out_dir, '/UMAP_markers'))

  marker_plots <- list()
  for (m in markers) {
    p <- ggplot2::ggplot(df, ggplot2::aes_string(x = colnames(df)[1],
                                                 y = colnames(df)[2])) +
      ggplot2::geom_point(ggplot2::aes_string(color = m),
                          alpha = 0.3,
                          size = 0.4) +
      ggplot2::scale_color_gradientn(m, colors = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))(50)) +
      ggplot2::theme_bw() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::ggtitle(paste(toupper(type), "-", name))

    suppressMessages(ggplot2::ggsave(p, filename = paste0(out_dir, '/UMAP_markers/UMAP_batches_labels', m, '.png')))
  }

}



#' Dimensionality reduction plots for a two-batch dataset before and after correction - colored by batch and label
#' @inheritParams plot_density
#' @inheritParams plot_dimred
#' @family plot
#' @export
plot_umap_labels <- function(uncorrected,
                             corrected,
                             filename = NULL,
                             seed = 473) {


  cyCombine:::missing_package("uwot", "CRAN")
  cyCombine:::missing_package("ggridges", "CRAN")
  cyCombine:::missing_package("cowplot", "CRAN")
  cyCombine:::missing_package("ggplot2", "CRAN")


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

  plots <- plots2 <- list()
  plot_count <- plot_count2 <- 1
  markers <- cyCombine::get_markers(uncorrected)
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
        ggplot2::ggplot(ggplot2::aes_string(x = colnames(df)[1],
                                            y = colnames(df)[2])) +
        ggplot2::geom_point(ggplot2::aes_string(color = plot),
                            alpha = 0.3,
                            size = 0.4,
                            shape = 1) +
        ggplot2::labs(plot) +
        ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(alpha = 1, size = 1))) +
        ggplot2::theme_bw() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::ggtitle(paste("UMAP -", name))

      plots[[plot_count]] <- plot
      plot_count <- plot_count + 1

    } else if (nplots == 3) {

      # Together plot
      plot <- df %>%
        ggplot2::ggplot(ggplot2::aes_string(x = colnames(df)[1],
                                            y = colnames(df)[2])) +
        ggplot2::geom_point(ggplot2::aes_string(color = plot),
                            alpha = 0.3,
                            size = 0.4,
                            shape = 1) +
        ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(alpha = 1, size = 1))) +
        ggplot2::labs(plot) +
        ggplot2::theme_bw() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::ggtitle(paste("UMAP -", name))

      plots[[plot_count]] <- plot
      plot_count <- plot_count + 1

      # Faceted plots
      plot <- df %>%
        ggplot2::ggplot(ggplot2::aes_string(x = colnames(df)[1],
                                            y = colnames(df)[2])) +
        ggplot2::facet_wrap(~Batch) +
        ggplot2::geom_point(ggplot2::aes_string(color = Label),
                            alpha = 0.3,
                            size = 0.4,
                            shape = 1) +
        ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(alpha = 1, size = 1))) +
        ggplot2::labs(colour = 'Label') +
        ggplot2::theme_bw() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::ggtitle(paste("UMAP -", name))

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
#' @family plot
#' @export
plot_save_two <- function(plot1, plot2, filename) {
  cyCombine:::missing_package("cowplot", "CRAN")

  plot <- cowplot::plot_grid(plot1, plot2, align = "v", scale = 0.9)
  cowplot::save_plot(filename = filename, plot, base_width = 12, base_height = 6)
}


#' Save four plots aligned with cowplot
#' @param plot1 Top left plot in output
#' @param plot2 Top right plot in output
#' @param plot3 Lower left plot in output
#' @param plot4 Lower right plot in output
#' @param filename Filename of output
#' @family plot
#' @export
plot_save_four <- function(plot1, plot2, plot3, plot4, filename) {
  cyCombine:::missing_package("cowplot", "CRAN")

  plot <- cowplot::plot_grid(plot1, plot2, plot3, plot4, align = "v", scale = 0.9, nrow = 2)
  cowplot::save_plot(filename = filename, plot, base_width = 12, base_height = 12)
}
