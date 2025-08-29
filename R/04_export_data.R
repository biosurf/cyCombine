#' Exporting a dataframe to SingleCellObject
#'
#' Conversion of dataframe into separate data and metadata objects for subsequent transformation.
#'  Rename variable names to fit the requirements of SCE-based tools
#'
#' @param df Tibble with expression values and metadata
#' @param markers Markers to include in exprs and counts object of SCE. If NULL, markers will be found using the \code{\link{get_markers}} function.
#' @param non_markers Non-markers to include as colData in SCE. If NULL, non_markers will be based on cyCombine::non_markers.
#' @param scatter The scatter columns to add to the SCE for FCS export.
#' @param sample_col It is the column name in the df that contains the sample names. Defaults to 'sample'.
#' @param panel Optional: Panel as a data.frame. Should have colnames Channel, Marker, Type unless otherwise specified in the panel_ args. Should be included if you want to store FCS files
#' @param panel_channel Optional: Only used if panel is given. It is the column name in the panel data that contains the channel names
#' @param panel_antigen Optional: Only used if panel is given. It is the column name in the panel data that contains the antigen names
#' @param panel_type Optional: Only used if panel is given. It is the column name in the panel data that contains the antigen types (none, state, type).
#'  "none" will be excluded from SCE. Set to NULL to disregard.
#' @param transform_cofactor The cofactor to use when reverse-transforming to raw counts
#' @param clean_names Cleans marker names.
#' @importFrom methods is
#' @family export
#' @examples
#' \dontrun{
#' sce <- df %>%
#'   df2SCE(markers = markers, non_markers = NULL, panel = panel)
#'   }
#' @export
df2SCE <- function(
    df,
    markers = NULL, non_markers = NULL,
    scatter = NULL,
    clean_names = FALSE,
    sample_col = "sample",
    panel = NULL, panel_channel = "Channel",
    panel_antigen = "Marker", panel_type = NULL,
    transform_cofactor = 5) {

  # Check for packages
  check_package("SingleCellExperiment", "Bioc")

  message("Converting dataframe to SingleCellExperiment object...")

  # Get the non markers and prepare column data
  if (is.null(non_markers)) {
    # Get non markers
    if (sample_col == "sample") {
      non_markers <- cyCombine::non_markers
    } else {
      non_markers <- c(cyCombine::non_markers, dplyr::all_of(sample_col))
    }
  }

  # Get the column data (from non markers)
  if (sample_col %in% colnames(df)) {
    colData <- df %>%
      dplyr::select(dplyr::any_of(non_markers))

    if (sample_col != 'sample_id') {
      colData <- colData %>%
        dplyr::rename(sample_id = sample_col)
    }

  } else {
    stop("Error, none of the non_markers/sample_col are available in the dataframe. You cannot make an SCE without sample names.")
  }

  # Get markers and check
  if (is.null(markers)) {
    # Get markers
    markers <- df %>%
      cyCombine::get_markers()
  }

  sapply(markers, function(x) {
    cyCombine:::check_colname(colnames(df), x, "df")})

  # Extract expression data and transpose to fit SCE format
  exprs <- t(dplyr::select(df, dplyr::all_of(union(scatter, markers))))


  # Prepare the experiment info table - first identify true meta data columns
  colCounting <- colData %>%
    dplyr::group_split(sample_id) %>%
    lapply(function(x) {(apply(x, 2, function(y) {length(table(y)) == 1}))})

  col_stability <- lapply(colCounting, function(z) {names(which(z))}) %>% unlist() %>% table()
  stable_cols <- names(which(col_stability == length(colCounting)))

  # Swap ordering
  stable_cols <- c("sample_id", stable_cols[stable_cols != "sample_id"])

  experiment_info <- colData %>%
    dplyr::select(dplyr::all_of(stable_cols)) %>%
    dplyr::group_by(sample_id) %>%
    dplyr::mutate(n_cells = dplyr::n()) %>%
    dplyr::distinct(sample_id, .keep_all = TRUE) %>%
    as.data.frame()



  # Prepare row data if available
  if (!is.null(panel) && is.data.frame(panel)) {
    # Check presence of panel's data names
    sapply(c(panel_channel, panel_antigen, panel_type), function(x) {
      cyCombine:::check_colname(colnames(panel), x, "panel")})

    rowData <- panel[match(rownames(exprs), panel[[panel_antigen]]), ]
    rm(panel)
    # Exclude none's
    if (!is(panel_type, "NULL")) {
      if ("none" %in% rowData[[panel_type]]){
        rowData <- rowData %>%
          dplyr::filter(.data[[panel_type]] != "none")
      }
    }


    # Change colnames to fit SCE standard
    if (panel_channel != "channel_name") {
      rowData <- rowData %>%
        dplyr::rename(channel_name = panel_channel)
    }
    if (panel_antigen != "marker_name") {
      rowData <- rowData %>%
        dplyr::rename(marker_name = panel_antigen)
    }
    if (panel_type != "marker_class" && !is(panel_type, "NULL")) {
      rowData <- rowData %>%
        dplyr::rename(marker_class = panel_type)
    }

    # Change marker names to exclude spaces and dashes
    if (clean_names) {
      rowData[["marker_name"]] <- rowData[["marker_name"]] %>%
        stringr::str_remove_all("^\\d+[A-Za-z]+_") %>%
        stringr::str_remove_all("[ _-]")
    }


    # Subset rowData to rows in exprs
    rowData <- rowData %>%
      dplyr::filter(.data[["marker_name"]] %in% rownames(exprs))
    rownames(rowData) <- rownames(exprs)[match(rownames(exprs), rowData[["marker_name"]])]

  } else {
    rowData <- NULL
    warning("To store as FCS files later, you should include panel information at this step.")
  }

  exprs <- exprs[c(scatter, markers), ]
  # Creating the SCE
  sce <- SingleCellExperiment::SingleCellExperiment(
    list(exprs = exprs,
         counts = exprs %>%
           t() %>%
           tibble::as_tibble() %>%
           transform_asinh(reverse = TRUE,
                           cofactor = transform_cofactor,
                           derand = FALSE,
                           .keep = TRUE,
                           markers = markers) %>%
           t()),
    colData = colData,
    rowData = rowData,
    metadata = list("experiment_info" = experiment_info)
  )
  message("Your SingleCellExperiment object is now created. The 'counts' assay contains reverse transformed expression values and 'exprs' contains expression values.")

  return(sce)
}

#' Convert SingelCellExperiment into flowSet and store as FCS files
#'
#' Wrapper function that makes it easier to go from a SCE to flowSet and written FCS files.
#'   The function uses CATALYST::sce2fcs and flowCore::write.flowSet to store the FCS files.
#'
#' @inheritParams CATALYST::sce2fcs
#' @inheritParams flowCore::write.flowSet
#' @param sce SingleCellExperiment to write to FCS files
#' @param outdir If given, the flowSet will be stored in FCS files
#' @param randomize (Default: FALSE) Logical determining whether counts are randomized for plotting purposes or not. Only works when assay = "counts" (default).
#' @family export
#' @examples
#' \dontrun{
#'  sce2FCS(sce, outdir = "fcs_files")
#'   }
#' @export
sce2FCS <- function(sce,
                    outdir = NULL,
                    split_by = "sample_id",
                    assay = "counts",
                    keep_dr = TRUE,
                    keep_cd = TRUE,
                    randomize = FALSE) {


  # Check CATALYST is installed
  check_package("CATALYST", repo = "Bioc")

  # If no panel information is available
  stopifnot("Your SingleCellExperiment should contain channel information to be stored.\n
            Consider rerunning df2sce with a panel to continue." = !is.null(CATALYST::channels(sce)))

  # Randomize counts
  if (randomize && assay == "counts") {
    SummarizedExperiment::assay(sce, "counts") <- SummarizedExperiment::assay(sce, "counts") %>%
      cyCombine:::randomize_matrix()
  } else if (randomize && assay != "counts") {
    warning("Please only use randomization with count values.")
  }

  # Convert to flowset
  fcs <- CATALYST::sce2fcs(sce,
                           keep_dr = keep_dr,
                           keep_cd = keep_cd,
                           split_by = split_by,
                           assay = assay)


  # Write FCS files
  if (!is.null(outdir)) {
    # Create output directory
    cyCombine:::check_make_dir(outdir)

    message("Writing FCS files to ", outdir)
    switch(class(fcs)[1],
           "flowFrame" = flowCore::write.FCS(fcs, filename = paste0(outdir, sce$sample_id[1], ".fcs")),
           "flowSet" = flowCore::write.flowSet(fcs, outdir = outdir))
  }

  return(fcs)
}


#' Function to remove specific keywords from a flowFrame
#'
#' @note flowCore parameters exprs keyword
#' @note Biobase pData
#' @noRd
remove_keywords <- function(ff, remove_kwd="flowCore_") {
  params <- flowCore::parameters(ff)
  pdata <- Biobase::pData(params)
  exp_mat <- flowCore::exprs(ff)
  old_kwds <- flowCore::keyword(ff)

  keyval <- list()
  unchanged_keywords <- flowCore::keyword(ff)
  keyval <- unchanged_keywords[!grepl("\\$P", names(unchanged_keywords))]

  new_pdata <- pdata

  for (i in 1:nrow(pdata)) {
    marker <- pdata$name[i]
    Rmax <- ceiling(max(exp_mat[, marker]))
    Rmin <- ceiling(min(exp_mat[, marker]))

    new_pdata[i, "name"] <- pdata[i, "name"]
    new_pdata[i, "desc"] <- pdata[i, "desc"]
    new_pdata[i, "range"] <- Rmax - Rmin + 1
    new_pdata[i, "minRange"] <- Rmin
    new_pdata[i, "maxRange"] <- Rmax
    keyval[[paste0("$P", i, "R")]] <- as.character(Rmax - Rmin + 1)
  }

  new_kwds <- names(old_kwds)[!grepl(remove_kwd, names(old_kwds))]

  for (key in new_kwds) {
    new_key <- key
    keyval[[new_key]] <- old_kwds[[new_key]]
    if (new_key == "r") {
      keyval[[new_key]] <- as.character(Rmax - Rmin + 1)
    }
  }
  keyval[["transformation"]] <- NULL

  new_ff <- ff
  flowCore::keyword(new_ff) <- keyval
  Biobase::pData(flowCore::parameters(new_ff)) <- new_pdata

  return(new_ff)
}

#' Process FCS files for FlowJo integration
#'
#' Process corrected and exported FCS files such that FlowJo can read them.
#'    Reading FCS files in R with flowCore adds 'flowCore_' keywords that are
#'    incompatible with FlowJo.
#'    The `process_fcs_files` function reads a directory of corrected files
#'    and exports them without the flowCore keywords.
#'
#' @param input_folder path to input folder
#' @param output_folder path to output folder
#'
#' @examples
#' \dontrun{
#'  process_fcs_files(
#'  input_folder = "path/to/input_folder",
#'  output_folder = "path/to/output_folder")
#'   }
#'
#' @note flowCore write.FCS
#'
#' @export
process_fcs_files <- function(input_folder, output_folder) {

  cyCombine:::missing_package("flowCore", repo = "Bioc")
  cyCombine:::missing_package("Biobase", repo = "Bioc")


  files <- list.files(input_folder, pattern = "\\.(fcs|FCS)$", full.names = TRUE)
  system(paste("mkdir -p", output_folder))

  for (file in files) {
    ff <- flowCore::read.FCS(file)
    new_ff <- remove_keywords(ff, remove_kwd="flowCore_")
    new_filename <- file.path(output_folder, basename(file))
    flowCore::write.FCS(new_ff, new_filename)
  }
}

#' Exporting a dataframe to Seurat Object
#'
#' Conversion of dataframe into Seurat object for downstream analysis.
#' Rename variable names to fit the requirements of Seurat-based tools
#'
#' @param df Tibble with expression values and metadata
#' @param markers Markers to include in data layer of Seurat object. If NULL, markers will be found using the \code{\link{get_markers}} function.
#' @param non_markers Non-markers to include as metadata in Seurat object. If NULL, non_markers will be based on cyCombine::non_markers.
#' @param scatter The scatter columns to add to the Seurat object.
#' @param sample_col It is the column name in the df that contains the sample names. Defaults to 'sample'.
#' @param panel Optional: Panel as a data.frame. Should have colnames Channel, Marker, Type unless otherwise specified in the panel_ args.
#' @param panel_channel Optional: Only used if panel is given. It is the column name in the panel data that contains the channel names
#' @param panel_antigen Optional: Only used if panel is given. It is the column name in the panel data that contains the antigen names
#' @param panel_type Optional: Only used if panel is given. It is the column name in the panel data that contains the antigen types (none, state, type).
#'  "none" will be excluded from Seurat object. Set to NULL to disregard.
#' @param transform_cofactor The cofactor to use when reverse-transforming to raw counts
#' @param clean_names Cleans marker names.
#' @param assay_name Name for the assay in the Seurat object. Defaults to 'RNA'.
#' @importFrom methods is
#' @family export
#' @examples
#' \dontrun{
#' seurat_obj <- df %>%
#'   df2Seurat(markers = markers, non_markers = NULL, panel = panel)
#'   }
#' @export
df2Seurat <- function(
    df,
    markers = NULL,
    non_markers = NULL,
    scatter = NULL,
    clean_names = FALSE,
    sample_col = "sample",
    panel = NULL, panel_channel = "Channel",
    panel_antigen = "Marker", panel_type = NULL,
    transform_cofactor = 5,
    assay_name = "RNA") {

  # Check for packages
  check_package("Seurat", "CRAN")

  message("Converting dataframe to Seurat object...")

  # Get the non markers and prepare metadata
  if (is.null(non_markers)) {
    # Get non markers
    if (sample_col == "sample") {
      non_markers <- cyCombine::non_markers
    } else {
      non_markers <- c(cyCombine::non_markers, dplyr::all_of(sample_col))
    }
  }

  # Get the metadata (from non markers)
  if (sample_col %in% colnames(df)) {
    metadata <- df %>%
      dplyr::select(dplyr::any_of(non_markers))

    if (sample_col != 'sample_id') {
      metadata <- metadata %>%
        dplyr::rename(sample_id = sample_col)
    }
    metadata <- as.data.frame(metadata)

  } else {
    stop("Error, none of the non_markers/sample_col are available in the dataframe. You cannot make a Seurat object without sample names.")
  }

  # Get markers and check
  if (is.null(markers)) {
    # Get markers
    markers <- cyCombine::get_markers(df)
  }

  sapply(markers, function(x) {
    cyCombine:::check_colname(colnames(df), x, "df")})

  # Extract expression data and transpose to fit Seurat format (features x cells)
  exprs <- t(dplyr::select(df, dplyr::all_of(union(scatter, markers))))

  # Create counts matrix (reverse transformed)
  counts <- exprs %>%
    t() %>%
    tibble::as_tibble() %>%
    transform_asinh(reverse = TRUE,
                    cofactor = transform_cofactor,
                    derand = FALSE,
                    .keep = TRUE,
                    markers = markers) %>%
    t()

  # Prepare the experiment info table - first identify true meta data columns
  colCounting <- metadata %>%
    dplyr::group_split(sample_id) %>%
    lapply(function(x) {(apply(x, 2, function(y) {length(table(y)) == 1}))})

  col_stability <- lapply(colCounting, function(z) {names(which(z))}) %>% unlist() %>% table()
  stable_cols <- names(which(col_stability == length(colCounting)))

  # Swap ordering
  stable_cols <- c("sample_id", stable_cols[stable_cols != "sample_id"])

  experiment_info <- metadata %>%
    dplyr::select(dplyr::all_of(stable_cols)) %>%
    dplyr::group_by(sample_id) %>%
    dplyr::mutate(n_cells = dplyr::n()) %>%
    dplyr::distinct(sample_id, .keep_all = TRUE) %>%
    as.data.frame()

  # Prepare feature data if panel is available
  if (!is.null(panel) && is.data.frame(panel)) {
    # Check presence of panel's data names
    sapply(c(panel_channel, panel_antigen, panel_type), function(x) {
      cyCombine:::check_colname(colnames(panel), x, "panel")})

    feature_data <- panel
    rm(panel)

    # Exclude none's
    if (!is(panel_type, "NULL")) {
      if ("none" %in% feature_data[[panel_type]]){
        feature_data <- feature_data %>%
          dplyr::filter(.data[[panel_type]] != "none")
      }
    }

    # Change colnames to fit standard naming
    if (panel_channel != "channel_name") {
      feature_data <- feature_data %>%
        dplyr::rename(channel_name = panel_channel)
    }
    if (panel_antigen != "marker_name") {
      feature_data <- feature_data %>%
        dplyr::rename(marker_name = panel_antigen)
    }
    if (panel_type != "marker_class" && !is(panel_type, "NULL")) {
      feature_data <- feature_data %>%
        dplyr::rename(marker_class = panel_type)
    }

    # Change marker names to exclude spaces and dashes
    if (clean_names) {
      feature_data[["marker_name"]] <- feature_data[["marker_name"]] %>%
        stringr::str_remove_all("^\\d+[A-Za-z]+_") %>%
        stringr::str_remove_all("[ _-]")
    }

    # Subset feature_data to rows in exprs
    feature_data <- feature_data %>%
      dplyr::filter(.data[["marker_name"]] %in% rownames(exprs))
    rownames(feature_data) <- rownames(exprs)[match(rownames(exprs), feature_data[["marker_name"]])]

  } else {
    feature_data <- NULL
    warning("Panel information not provided. Feature metadata will be minimal.")
  }

  # Set cell names (barcodes) for Seurat
  cell_names <- paste0("Cell_", seq_len(ncol(exprs)))
  colnames(exprs) <- cell_names
  colnames(counts) <- cell_names
  rownames(metadata) <- cell_names

  # Creating the Seurat object
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = counts,
    meta.data = metadata,
    assay = assay_name
  )

  # Add expression data to data layer (Seurat v5 terminology)
  seurat_obj <- Seurat::SetAssayData(
    object = seurat_obj,
    layer = "data",
    new.data = exprs,
    assay = assay_name
  )

  # Add feature metadata if available
  if (!is.null(feature_data)) {
    # Ensure feature names match
    feature_data <- feature_data[rownames(seurat_obj), , drop = FALSE]
    # Use Seurat v5 method for adding feature metadata
    seurat_obj[[assay_name]] <- Seurat::AddMetaData(
      object = seurat_obj[[assay_name]],
      metadata = feature_data
    )
  }

  # Add experiment info to misc slot
  seurat_obj@misc[["experiment_info"]] <- experiment_info

  message("Your Seurat object is now created. The 'counts' layer contains reverse transformed expression values and 'data' layer contains expression values.")

  return(seurat_obj)
}

