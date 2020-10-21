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

  # test <- data_dir %>%
  #   compile_fcs(meta_filename = "CyTOF samples cohort.xlsx",sample_col = "FCS_name")
  preprocessed <- preprocess(data_dir = data_dir,
                             meta_filename = "CyTOF samples cohort.xlsx",
                             markers = markers,
                             down_sample = TRUE,
                             sample_size = 100000,
                             seed = 474,
                             cofactor = 5)

  save(preprocessed, file = "_data/01_preprocessed_100k.Rdata")
  # Batch correct

  corrected <- preprocessed %>%
    batch_correct()

  # Save result
  save(corrected, file = "_data/02_corrected_100k.Rdata")
  })

  # 100k: 1.8 min
  # 300k: 4.2 min
  # 500k: 6.7 min
  # 700k: 9.2 min

  ### Panel1 ----
  load("_data/raw/DFCI_panel2_data.Rdata")

  markers <- panel2_data$all_markers
  preprocessed <- panel2_data$fcs_raw %>%
    convert_flowset(sample_ids = panel2_data$sample_ids,
                    batch_ids = panel2_data$batch_ids,
                    down_sample = TRUE,
                    sample_size = 700000,
                    seed = 473) %>%
    transform_asinh(panel2_data$all_markers)

  save(preprocessed, file = "_data/01_dfci2_preprocessed_700k.Rdata")
  som2 <- preprocessed %>%
    scale_expr() %>%
    create_som(seed = 473)
  # #
  # corrected <- preprocessed %>%
  #   correct_data(som_classes = som3$unit.classif)
  #



  # Run batch correction
  corrected <- preprocessed %>%
    batch_correct()
  # Save result
  save(corrected, file = "_data/02_dfci2_corrected_700k.Rdata")

  .
  ### Plotting ----

  plot_density(uncorrected = preprocessed,
                corrected = corrected,
                markers = markers,
                filename = "figs/03_dfci2_densities_withcovar_700k.png")

  # Down-sample
  preprocessed_sliced <- preprocessed %>%
    slice_sample(n = 10000)

  corrected_sliced <- corrected %>%
    semi_join(preprocessed_sliced, by = "id")

  # PCA plot uncorrected
  pca1 <- preprocessed_sliced %>%
    plot_dimred("uncorrected", type = "pca")


  # PCA plot corrected
  pca2 <- corrected_sliced %>%
    plot_dimred("corrected", type = "pca")
  plot_save_two(pca1, pca2, filename = "figs/02_dfci2_pca_700k.png")




  # UMAP plot uncorrected
  umap1 <- preprocessed_sliced %>%
    plot_dimred(name = "uncorrected", type = "umap")

  # UMAP plot corrected
  umap2 <- corrected_sliced %>%
    plot_dimred(name = "corrected", type = "umap")

  plot_save_two(umap1, umap2, filename = "figs/02_dfci2_umap_700k.png")
  # library(patchwork)
  # umap1 + umap2



  ### Evaluate performance ----
  # load("_data/02_corrected_700k.Rdata")
  # load("_data/02_preprocessed_700k.Rdata")
  load("_data/01_dfci2_preprocessed_700k.Rdata")
  load("_data/02_dfci2_corrected_700k.Rdata")

  som3 <- corrected %>%
    # scale_expr() %>%
    create_som(seed = 473, xdim = 5)

  # Run clustering
  corrected <- corrected %>%
    dplyr::mutate(label = som$unit.classif)

  threshold <- corrected %>%
    count(label) %>%
    filter(n < 3000)
  corrected_t <- corrected %>%
    anti_join(threshold, by = "label")

  preprocessed <- preprocessed %>% select(-label)
  #run_flowsom(., k = 10))
  preprocessed <- corrected %>%
    dplyr::select(id, label) %>%
    dplyr::left_join(preprocessed, by = "id")# %>%
    # anti_join(threshold, by = "label")
  preprocessed_t <- preprocessed %>%
    anti_join(threshold, by = "label")

    # dplyr::mutate(label = run_flowsom(., k = 100))

  # Compute EMD reduction
  emd_val <- preprocessed %>%
    cyCombine::evaluate_emd(corrected)
  emd_val2 <- preprocessed %>%
    cyCombine::evaluate_emd2(corrected)
  save(emd_val, emd_val2, file = "_data/04_dfci2_emd_700k.Rdata")
  # Compute LISI score

  # Down-sample
  # preprocessed_sliced <- preprocessed %>%
  #   slice_sample(n = 10000)
  #
  # corrected_sliced <- corrected %>%
  #   semi_join(preprocessed_umap, by = "id")
  #
  # lisi_cor <- corrected_sliced %>%
  #   evaluate_lisi()
  # lisi_prep <- preprocessed_sliced %>%
  #   evaluate_lisi()

### CD4 ----

  library(dplyr)
  umap <- preprocessed_sliced %>%
    mutate(outlier = case_when(label == 34 ~ "Label 34",
                               label== 23 ~ "Label 23",
                               TRUE ~ "other"))
set.seed(473)
umap_d <- umap %>%
  select(get_markers(.)) %>%
  uwot::umap(n_neighbors = 15, min_dist = 0.2, metric = "euclidean")

colnames(umap_d) <- c("UMAP1", "UMAP2")
df <- umap_d %>%
  tibble::as_tibble() %>%
  dplyr::mutate(outlier = umap$outlier)

plot <- df %>%
  ggplot(aes_string(x = colnames(df)[1], y = colnames(df)[2])) +
  geom_point(aes_string(color = "outlier"), alpha = 0.3, size = 0.4, shape = 1) +
  viridis::scale_color_viridis(discrete = TRUE) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle(paste(toupper(type), "-", name))
plot
}
