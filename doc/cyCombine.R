## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
if(FALSE){
  pkgload::load_all()
  library(cyCombine)
  library(magrittr)

  ### From raw fcs files ----
  system.time({
  data_dir <- "~/Documents/thesis/raw/fcs"
  markers <- c("CD20", "CD3", "CD27", "CD45RA", "CD279", "CD5", "CD19", "CD14", "CD45RO", "GranzymeA", "GranzymeK", "FCRL6", "CD355", "CD152", "CD69", "CD33", "CD4", "CD337", "CD8", "CD197", "LAG3", "CD56", "CD137", "CD161", "FoxP3", "CD80", "CD270", "CD275", "CD134", "CD278", "CD127", "KLRG1", "CD25", "HLADR", "TBet", "XCL1")
  preprocessed <- preprocess(data_dir = data_dir,
                         meta_filename = "CyTOF samples cohort.xlsx",
                         markers = markers,
                         down_sample = TRUE,
                         sample_size = 700000,
                         seed = 473,
                         cofactor = 5)

  save(preprocessed, file = "data/01_preprocessed_700k.Rdata")
  # Batch correct

  corrected <- preprocessed %>%
    batch_correct()

  # Save result
  save(corrected, file = "data/02_corrected_700k.Rdata")
  })

  # 100: 1.8 min
  # 300: 4.2 min
  # 500: 6.7 min
  # 700: 9.2 min

  ### From .Rdata file ----
  load("data/raw/DFCI_panel1_data.Rdata")

  preprocessed <- panel1_data$fcs_raw %>%
    convert_flowset(sample_ids = panel1_data$sample_ids,
                    batch_ids = panel1_data$batch_ids,
                    down_sample = TRUE,
                    sample_size = 100000,
                    seed = 473) %>%
    transform_asinh(panel1_data$all_markers)

  som <- preprocessed %>%
    scale_expr() %>%
    create_som(seed = 473)
  #
  corrected <- preprocessed %>%
    correct_data(som_classes = som$unit.classif)
  #



  # Run batch correction
  corrected <- preprocessed %>%
    batch_correct()
  # Save result
  save(corrected, file = "data/02_panel1_corrected.Rdata")


  ### Plotting ----
  plot_density(uncorrected = preprocessed,
                corrected = corrected,
                markers = markers,
                filename = 'figs/02_densities_withcovar_500k.png')

  # PCA plot uncorrected
  pca1 <- preprocessed %>%
    plot_dimred('uncorrected', type = 'pca')


  # UMAP plot uncorrected
  umap1 <- preprocessed %>%
    plot_dimred('uncorrected', type = 'umap')



  # PCA plot corrected
  pca2 <- corrected %>%
    plot_dimred('corrected', type = 'pca')
  plot_save_two(pca1, pca2, filename = 'figs/02_pca_500k.png')


  # UMAP plot corrected
  umap2 <- plot_dimred(corrected, 'corrected', type = 'umap')
  plot_save_two(umap1, umap2, filename = 'figs/02_raw_umap.png')




  ### Evaluate performance ----
  load("data/02_corrected_500k.Rdata")
  load("data/02_preprocessed_500k.Rdata")

  # Run clustering
  corrected <- corrected %>%
    dplyr::mutate(label = run_flowsom(., k = 10))
  preprocessed <- preprocessed %>%
    dplyr::mutate(label = run_flowsom(., k = 10))

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
