
#' Run analysis of a batch correction
#'
#' This function only runs under the assumption that the batch correction has been run with the cyCombine workflow in mind.
#'   Some result preparation can be necessary, if cyCombine was not used during batch correction.
#'   This function assumes that the uncorrected data us stored under the name \code{"{data_dir}/{tool}_{data}{variant}{uncorrected_extension}.RDS"} and
#'   the corrected data is stored with the name \code{"{data_dir}/{tool}_{data}{variant}{corrected_extension}.RDS"}.
#'   The analysis currently encompass marker-wise density plots, UMAP of corrected vs uncorrected, and Earth Movers Distance.
#'
#' @param tool The name of the tool used to batch correct
#' @param data The name of the data used
#' @param data_dir The location of the uncorrected and corrected data
#' @param uncorrected_extension The extension used to name the uncorrected data. Default: "_uncorrected
#' @param corrected_extension The extension used to name the corrected data. Default: "_corrected"
#' @param variant Optional: A parameter to set a variant name of an experiment
#' @param use_cycombine_uncor If TRUE, the uncorrected data made by cyCombine will be used.
#' @param restart If TRUE, the SOM grid will be calculated even if it has been computed and stored previously
#' @param md Optional: Metadata filename. Currently not useful
#' @param panel Optional: If given, it will be used to define markers. Otherwise the function \code{\link{get_markers}} will be used
#' @param markers Optional: Manually define markers to use in plots and performance metrics
#' @param celltype_col Optional: If the cell types are known, specify which column they are defined in. If NULL, a clustering will be run.
#' @param segment Optional: Run only a specific segment of the analysis. Options include: "emd", "density", "umap"
#' @param gridsize The gridsize to use when clustering. Only used if no celltype_col is given
#' @param seed The seed to use when creating the UMAP
#' @inheritParams create_som
#'
#' @examples
#' run_analysis(tool = "cycombine",
#'  data = "dfci1",
#'  data_dir = "_data")
#'  run_analysis(tool = "cycombine",
#'  data = "FR-FCM-ZY34",
#'  data_dir = "_data",
#'  variant = "_p3",
#'  panel = "/attachments/MC_panel3.xlsx")
#' @export
run_analysis <- function(tool,
                         data,
                         data_dir,
                         uncorrected_extension = "_uncorrected",
                         corrected_extension = "_corrected",
                         variant = NULL,
                         uncorrected_variant = NULL,
                         use_cycombine_uncor = FALSE,
                         restart = FALSE,
                         md = NULL,
                         panel = NULL,
                         markers = NULL,
                         celltype_col = NULL,
                         segment = "",
                         binSize = 0.1,
                         rlen = 10,
                         gridsize = 8,
                         seed = 473,
                         umap_size = 20000){
  # Load metadata
  if(!is.null(md)){
    if(endsWith(md, "csv")){
      md <- stringr::str_c(data_dir, "/", md) %>%
        readr::read_csv()
    } else{
      md <- stringr::str_c(data_dir, "/", md) %>%
        readxl::read_xlsx()
    }
  }
  # Load panel
  if(!is.null(panel)){
    if(endsWith(panel, "csv")){
      panel <- stringr::str_c(data_dir, "/", panel) %>%
        readr::read_csv()
    } else{
      panel <- stringr::str_c(data_dir, "/", panel) %>%
        readxl::read_xlsx()
    }
  }

  # dir + project
  project <- paste0(tool, "_", data, variant)
  projdir <- paste0(data_dir, "/", project)

  message("Loading data..")

  # Load data
  if(use_cycombine_uncor){
    uncorrected <- readRDS(paste0(data_dir, "/cycombine_", data, uncorrected_variant, uncorrected_extension, ".RDS"))
  } else{
    uncorrected <- readRDS(paste0(projdir, uncorrected_extension, ".RDS"))
  }
  corrected <- readRDS(paste0(projdir, corrected_extension, ".RDS"))


  # Get markers
  if(is.null(markers)){
    if(is.null(panel)){
      markers <- cyCombine::get_markers(uncorrected)
    } else{
      markers <- panel %>%
        dplyr::filter(marker_class %!in% c("none")) %>%
        dplyr::pull(antigen) %>%
        stringr::str_remove_all("[ _-]")
    }
  }


  # Plotting ----
  if (!dir.exists("figs")){
    dir.create("figs")
  }
  if(any(segment %in% c("", "density"))){

    message("Creating density plots..")
    # Density plots
    suppressMessages(
      cyCombine::plot_density(uncorrected = uncorrected,
                              corrected = corrected,
                              markers = markers,
                              filename = paste0("figs/", project, "_densities.png"))
    )
  }

  # UMAP
  if(any(segment %in% c("", "umap"))){
    message("Creating UMAPs..")

    # Down-sample
    set.seed(seed)
    uncorrected_sliced <- uncorrected %>%
      dplyr::slice_sample(n = umap_size)

    corrected_sliced <- corrected %>%
      dplyr::semi_join(uncorrected_sliced, by = "id")

    # UMAP plot uncorrected
    umap1 <- uncorrected_sliced %>%
      cyCombine::plot_dimred(name = "uncorrected", type = "umap", markers = markers)

    # UMAP plot corrected
    umap2 <- corrected_sliced %>%
      cyCombine::plot_dimred(name = "corrected", type = "umap", markers = markers)

    cyCombine::plot_save_two(umap1, umap2, filename = paste0("figs/", project, "_umap.png"))
  }

  # Evaluate ----
  if(any(segment %in% c("", "emd"))){

    if(is.null(celltype_col)){
      # Cluster and add labels
      if (file.exists(paste0(projdir, "_som.RDS")) & !restart){
        message("Loading SOMgrid..")
        labels <- readRDS(paste0(projdir, "_som.RDS"))
        # labels <- som_$unit.classif
      } else{
        labels <- corrected %>%
          cyCombine::create_som(seed = seed,
                                rlen = rlen,
                                xdim = gridsize,
                                ydim = gridsize,
                                markers = markers)
          saveRDS(labels, file = paste0(projdir, "_som.RDS"))
      }
      # Add labels
      corrected <- corrected %>%
        dplyr::mutate(som = labels)
      celltype_col <- "som"
    } else if (celltype_col %!in% colnames(corrected)) stop(paste0("Column '", celltype_col, "' was not found in 'Corrected'. Please set 'celltype_col' to NULL or a column containing the celltypes or cluster labels."))
    message("Adding labels to data..")

    if(celltype_col %!in% colnames(uncorrected)){
      uncorrected <- corrected %>%
        dplyr::select(id, all_of(celltype_col)) %>%
        dplyr::left_join(uncorrected, by = "id")
    }



    message("Evaluating Earth Movers Distance..")
    emd_val <- uncorrected %>%
      cyCombine::evaluate_emd(corrected,
                              binSize = binSize,
                              markers = markers,
                              cell_col = celltype_col)

    message("Saving results..")
    ggsave(filename = paste0("figs/", project, "_violin.png"),
           plot = emd_val$violin, device = "png")
    ggsave(filename = paste0("figs/", project, "_scatterplot.png"),
           plot = emd_val$scatterplot, device = "png")
    saveRDS(emd_val, file = paste0(projdir, "_emd.RDS"))
  }

  message("Done!")
}
