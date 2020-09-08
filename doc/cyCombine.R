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
  data_dir <- "~/Documents/thesis/raw/fcs"
  markers <- c("CD20", "CD3", "CD27", "CD45RA", "CD279", "CD5", "CD19", "CD14", "CD45RO", "GranzymeA", "GranzymeK", "FCRL6", "CD355", "CD152", "CD69", "CD33", "CD4", "CD337", "CD8", "CD197", "LAG3", "CD56", "CD137", "CD161", "FoxP3", "CD80", "CD270", "CD275", "CD134", "CD278", "CD127", "KLRG1", "CD25", "HLADR", "TBet", "XCL1")
  fcs_data <- preprocess(data_dir = data_dir,
                         markers = markers,
                         sample_size = 100000,
                         seed = 473)


  raw_flowset <- data_dir %>%
    compile_fcs()
  fcs_data <- raw_flowset$fcs_raw %>%
    convert_flowset(batch_ids = raw_flowset$batch_ids,
                    sample_ids = raw_flowset$sample_ids,
                    down_sample = TRUE,
                    sample_size = 300000,
                    seed = 473)

  fcs_preprocessed <- fcs_data %>%
    transform_asinh(markers)
  # fcs_corrected <- fcs_preprocessed %>%
    # batch_correct()
  som2 <- fcs_data %>%
    scale_expr() %>%
    create_som()

  corrected_data <- fcs_data %>%
    correct_data(som_classes = som2$unit.classif)


  ### Preprocessing ----
  # load("data/raw/DFCI_panel1_data.Rdata")

  panel2 <- panel1_data$fcs_raw %>%
    convert_flowset(sample_ids = panel1_data$sample_ids,
                    batch_ids = panel1_data$batch_ids,
                    down_sample = TRUE,
                    sample_size = 100000,
                    seed = 473)

  # panel2 <- panel2_data$fcs_raw %>%
  #   prepare_flowset(batch_ids = panel2_data$batch_ids,
  #                   sample_ids = panel2_data$sample_ids)
  preprocessed <- panel2 %>%
    transform_asinh(panel1_data$all_markers)

  # preprocessed_data <- fcs_data %>%
  #   preprocess(markers = markers)

  save(preprocessed, file = "data/01_panel2_preprocessed.Rdata")


  ### Batch correction ----

  # Load preprocessed data
  load("data/01_panel2_preprocessed.Rdata")
  # Run batch correction
  corrected <- preprocessed %>%
    batch_correct()
  # Save result
  save(corrected_data, file = "data/02_panel2_corrected.Rdata")

  # som <- preprocessed_data %>%
  #   scale_expr() %>%
  #   create_som()
  #
  # corrected_data <- preprocessed_data %>%
  #   correct_data(som_classes = som$unit.classif)


  ### Plotting ----
  density_plots(uncorrected = preprocessed,
                corrected = corrected_data,
                markers = panel2_data$all_markers,
                filename = 'figs/02_panel2_2_densities_withcovar.png')

  # PCA plot uncorrected
  pca1 <- preprocessed %>%
    dimred_plot('uncorrected', type = 'pca')


  # UMAP plot uncorrected
  umap1 <- preprocessed_data %>%
    dimred_plot('uncorrected', type = 'umap')



  # PCA plot corrected
  pca2 <- corrected_data %>%
    dimred_plot('corrected', type = 'pca')
  save_two_plots(pca1, pca2, filename = 'figs/02_raw_pca.png')


  # UMAP plot corrected
  umap2 <- dimred_plot(corrected_data, 'corrected', type = 'umap')
  save_two_plots(umap1, umap2, filename = 'figs/02_raw_umap.png')

}
