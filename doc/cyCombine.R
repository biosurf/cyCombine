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
  fcs_preprocessed <- preprocess(data_dir = data_dir,
                         meta_filename = "CyTOF samples cohort.xlsx",
                         markers = markers,
                         down_sample = TRUE,
                         sample_size = 300000,
                         seed = 473,
                         cofactor = 5)
  # Batch correct
  fcs_corrected <- fcs_preprocessed %>%
    batch_correct()

  # Save result
  save(fcs_corrected, file = "data/02_fcs_corrected.Rdata")



  ### From .Rdata file ----
  load("data/raw/DFCI_panel1_data.Rdata")

  preprocessed <- panel1_data$fcs_raw %>%
    convert_flowset(sample_ids = panel1_data$sample_ids,
                    batch_ids = panel1_data$batch_ids,
                    down_sample = TRUE,
                    sample_size = 100000,
                    seed = 473) %>%
    transform_asinh(panel1_data$all_markers)

  # som <- preprocessed %>%
  #   scale_expr() %>%
  #   create_som()
  #
  # corrected <- preprocessed %>%
  #   correct_data(som_classes = som$unit.classif)
  #
  # corrected2 <- preprocessed %>%
  #   correct_data_prev(som_classes = som$unit.classif)

  # Run batch correction
  corrected <- preprocessed %>%
    batch_correct()
  # Save result
  save(corrected, file = "data/02_panel1_corrected.Rdata")


  ### Plotting ----
  plot_density(uncorrected = preprocessed,
                corrected = corrected,
                markers = panel1_data$all_markers,
                filename = 'figs/02_panel1_densities_withcovar.png')

  # PCA plot uncorrected
  pca1 <- preprocessed %>%
    plot_dimred('uncorrected', type = 'pca')


  # UMAP plot uncorrected
  umap1 <- preprocessed %>%
    plot_dimred('uncorrected', type = 'umap')



  # PCA plot corrected
  pca2 <- corrected %>%
    plot_dimred('corrected', type = 'pca')
  plot_save_two(pca1, pca2, filename = 'figs/02_raw_pca.png')


  # UMAP plot corrected
  umap2 <- plot_dimred(corrected, 'corrected', type = 'umap')
  plot_save_two(umap1, umap2, filename = 'figs/02_raw_umap.png')




  ### Evaluate performance ----
  load("data/02_corrected.Rdata")
  load("data/01_preprocessed.Rdata")

  # Run clustering
  corrected <- corrected %>%
    dplyr::mutate(label = run_flowsom(.))
  preprocessed <- preprocessed %>%
    dplyr::mutate(label = run_flowsom(.))

  # Compute LISI score
  lisi_cor <- corrected %>%
    dplyr::slice_sample(n = 50000) %>%
    evaluate_lisi()
  lisi_prep <- preprocessed %>%
    dplyr::slice_sample(n = 50000) %>%
    evaluate_lisi()

  # Compute EMD reduction
  emd_val <- preprocessed %>%
    evaluate_emd(corrected)

}
