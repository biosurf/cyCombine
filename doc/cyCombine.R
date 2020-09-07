## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
if(FALSE){
  library(cyCombine)
  library(magrittr)

  ### From raw fcs files ----
  data_dir <- paste0("~/Documents/thesis/raw/fcs")
  markers <- c("CD20", "CD3", "CD27", "CD45RA", "CD279", "CD5", "CD19", "CD14", "CD45RO", "GranzymeA", "GranzymeK", "FCRL6", "CD355", "CD152", "CD69", "CD33", "CD4", "CD337", "CD8", "CD197", "LAG3", "CD56", "CD137", "CD161", "FoxP3", "CD80", "CD270", "CD275", "CD134", "CD278", "CD127", "KLRG1", "CD25", "HLADR", "TBet", "XCL1")
  fcs_data <- data_dir %>%
    compile_fcs(down_sample = FALSE)
  fcs_preprocessed <- fcs_data %>%
    transform_asinh(markers)
  # fcs_corrected <- fcs_preprocessed %>%
    # batch_correct(markers)
  som2 <- fcs_preprocessed %>%
    scale_expr() %>%
    create_som()

  corrected_data <- fcs_preprocessed %>%
    correct_data(som_classes = som2$unit.classif)


  ### Preprocessing ----
  # load("data/raw/DFCI_panel2_data.Rdata")
  panel2 <- panel2_data$fcs_raw %>%
    prepare_flowset(batch_ids = panel2_data$batch_ids,
                    sample_ids = panel2_data$sample_ids)
  preprocessed <- panel2 %>%
    transform_asinh(panel2_data$markers)

  preprocessed_data <- fcs_data %>%
    preprocess(markers = markers)

  save(preprocessed_data, file = "data/01_panel2_preprocessed.Rdata")


  ### Batch correction ----

  # Load preprocessed data
  load("data/01_panel2_preprocessed.Rdata")
  # Run batch correction
  corrected_data <- preprocessed_data %>%
    batch_correct(markers = markers)
  # Save result
  save(corrected_data, file = "data/02_panel2_corrected.Rdata")

  som <- preprocessed_data %>%
    scale_expr() %>%
    create_som()

  corrected_data <- preprocessed_data %>%
    correct_data(som_classes = som$unit.classif)


  ### Plotting ----
  density_plots(uncorrected = preprocessed_data,
                corrected = corrected_data,
                markers = markers,
                filename = 'figs/02_panel2_densities_withcovar.png')

  # PCA plot uncorrected
  pca1 <- preprocessed_data %>%
    dimred_plot('uncorrected', type = 'pca')


  # UMAP plot uncorrected
  umap1 <- preprocessed_data %>%
    dimred_plot('uncorrected', type = 'umap')



  # PCA plot corrected
  pca2 <- corrected_data %>%
    dimred_plot('corrected', type = 'pca')
  save_two_plots(pca1, pca2, filename = 'figs/02_panel2_pca.png')


  # UMAP plot corrected
  umap2 <- dimred_plot(corrected_data, 'corrected', type = 'umap')
  save_two_plots(umap1, umap2, filename = 'figs/02_umap.png')

}
