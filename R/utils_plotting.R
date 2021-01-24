#### Plotting functions ----

# @import ggplot2
# @import ggridges
# @import cowplot

# @importFrom dplyr select_if bind_rows rename
#' Density ridges for two sets
#' @import ggplot2
#' @examples 
#' plot_density(uncorrected, corrected, y = 'batch', filename = 'my/dir/batchcor_plot.pdf')
#' plot_density(imputed1, imputed2, y = 'Type', dataset_names = paste('Panel', 1:2), filename = 'my/dir/merging_plot.pdf')
#' @export
plot_density <- function(uncorrected, corrected, markers = NULL, filename = NULL, y = "batch", xlim = 10, dataset_names = NULL, ncol = 6) {

  # Check for packages
  missing_package("ggridges", "CRAN")
  missing_package("ggplot2", "CRAN")
  missing_package("cowplot", "CRAN")
  

  if (is.null(markers)) {
    markers <- uncorrected %>%
      get_markers()
  }
  
  if (is.null(dataset_names)) {
    dataset_names <- c('Uncorrected', 'Corrected')
  } else if (length(dataset_names) != 2) {
    dataset_names <- c('Dataset 1', 'Dataset 2')
  }

  if (y == 'Type') {
    batch1 <- dataset_names[1]
    batch2 <- dataset_names[2]
  } else {
    batch1 <- as.factor(uncorrected[[y]])
    batch2 <- as.factor(corrected[[y]])
  }
  
  
  uncorrected <- uncorrected %>%
    dplyr::select(all_of(markers)) %>%
    dplyr::mutate(Type = dataset_names[1],
                  batch = batch1)

  df <- corrected %>%
    dplyr::select(all_of(markers)) %>%
    dplyr::mutate(Type = dataset_names[2],
                  batch = batch2) %>%
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
      ggplot(aes_string(x = markers[c], y = 'batch')) +
      ggridges::geom_density_ridges(aes(color = .data$Type, fill = .data$Type), alpha = 0.4) +
      coord_cartesian(xlim = c(-1, xlim)) +
      labs(y = y) +
      theme_bw()
  }

  # p <- df %>%
  #   pivot_longer(cols = all_of(markers), names_to = "Marker") %>%
  #   ggplot(aes(x = value, y = batch)) +
  #   facet_wrap(~Marker, ncol = 6) +
  #   ggridges::geom_density_ridges(aes(color = .data$Type, fill = .data$Type), alpha = 0.4) +
  #   coord_cartesian(xlim = c(-1, xlim)) +
  #   labs(y = y) +
  #   theme_bw()


  # Save the plots
  # ggsave(filename = filename, plot = p,
  # device = "png", width = 28, height = 40)
  
  if (!is.null(filename)) {
    cowplot::save_plot(filename, cowplot::plot_grid(plotlist = p, ncol = ncol), base_width = 28, base_height = length(markers)/1.3)
  } else {
    return(cowplot::plot_grid(plotlist = p, ncol = ncol))
  }
}



# @importFrom uwot umap
#' Dimensionality reduction plot
#' @export
plot_dimred <- function(df, name, type = "pca", plot = "batch", markers = NULL, seed = 473, return_coord = F) {

  missing_package("uwot", "CRAN")
  missing_package("ggplot2", "CRAN")
  missing_package("ggridges", "CRAN")

  
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
      prcomp(scale. = TRUE, center = TRUE)

    # Make dataframe with output
    if (plot == "batch") {
      df <- pca$x %>%
        tibble::as_tibble() %>%
        dplyr::mutate(Batch = as.factor(Batch))
      #cbind.data.frame(pca$x, as.factor(batch_ids)); colnames(df)[ncol(df)] <- "Batch"
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
        dplyr::mutate(Batch = Batch) #%>%
        # dplyr::rename(UMAP1 = V1,
                      # UMAP2 = V2)
      # cbind.data.frame(umap, as.factor(batch_ids)); colnames(df) <- c("UMAP1", "UMAP2", "Batch")
    } else {
      df <- cbind.data.frame(umap, Batch, df[,plot]); colnames(df) <- c("UMAP1", "UMAP2", "Batch", plot)
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
      geom_point(aes_string(color = plot), alpha = 0.3, size = 0.4) +
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
#' @export
plot_dimred_full <- function(df, name, type = "pca", markers = NULL, seed = 473, out_dir = NULL) {

  missing_package("uwot", "CRAN")
  missing_package("ggplot2", "CRAN")
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
    ggplot(aes_string(x = colnames(df)[1], y = colnames(df)[2])) +
    geom_point(aes_string(color = "Batch"), alpha = 0.3, size = 0.4, shape = 1) +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle(paste(toupper(type), "-", name))

  label_plot <- df %>%
    ggplot(aes_string(x = colnames(df)[1], y = colnames(df)[2])) +
    geom_point(aes_string(color = "Label"), alpha = 0.3, size = 0.4, shape = 1) +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle(paste(toupper(type), "-", name))

  # Saving plots
  cyCombine::plot_save_two(batch_plot, label_plot, paste0(out_dir, '/UMAP_batches_labels.png'))

  # Marker plots
  cyCombine::check_make_dir(paste0(out_dir, '/UMAP_markers'))

  marker_plots <- list()
  for (m in markers) {
    p <- ggplot(df, aes_string(x = colnames(df)[1], y = colnames(df)[2])) +
      geom_point(aes_string(color = m), alpha = 0.3, size = 0.4) +
      scale_color_gradientn(m, colors = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))(50)) +
      theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
      ggtitle(paste(toupper(type), "-", name))

    suppressMessages(ggsave(p, filename = paste0(out_dir, '/UMAP_markers/UMAP_batches_labels', m, '.png')))
  }

}



# @importFrom uwot umap
#' Dimensionality reduction plots for a two-batch dataset before and after correction - colored by batch and label
#' @export
plot_umap_labels <- function(uncorrected, corrected, filename = NULL) {
  
  missing_package("uwot", "CRAN")
  missing_package("ggplot2", "CRAN")
  missing_package("ggridges", "CRAN")
  
  
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
  
  seed = 473
  plots <- plots2 <- list(); plot_count = plot_count2 = 1; markers = get_markers(uncorrected)
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
    
    
    firstup <- function(x) {
      substr(x, 1, 1) <- toupper(substr(x, 1, 1))
      x
    }
    
    # Make the plot
    if (nplots == 1) {
      plot <- df %>%
        ggplot(aes_string(x = colnames(df)[1], y = colnames(df)[2])) +
        geom_point(aes_string(color = plot), alpha = 0.3, size = 0.4, shape = 1) +
        labs(plot) +
        guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) +
        theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + 
        ggtitle(paste("UMAP -", name))
      
      plots[[plot_count]] <- plot
      plot_count = plot_count + 1
      
    } else if (nplots == 3) {
      
      # Together plot
      plot <- df %>%
        ggplot(aes_string(x = colnames(df)[1], y = colnames(df)[2])) +
        geom_point(aes_string(color = plot), alpha = 0.3, size = 0.4, shape = 1) +
        guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) +
        labs(plot) +
        theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
        ggtitle(paste("UMAP -", name))
      
      plots[[plot_count]] <- plot
      plot_count = plot_count + 1
      
      # Faceted plots
      plot <- df %>%
        ggplot(aes_string(x = colnames(df)[1], y = colnames(df)[2])) +
        facet_wrap(~Batch) +
        geom_point(aes_string(color = Label), alpha = 0.3, size = 0.4, shape = 1) +
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
