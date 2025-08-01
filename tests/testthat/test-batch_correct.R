data_dir <- system.file("extdata", package = "cyCombine")
df <- prepare_data(data_dir, transform = FALSE)
markers <- get_markers(df)
df <- df[, markers]
df$batch <- c(rep(1, 5000), rep(2, 5000))
df$sample <- df$batch

sce <- suppressWarnings(df2SCE(df))
seu <- suppressWarnings(df2Seurat(df))

test_that("batch correct works", {
  corrected <- batch_correct(df, xdim = 1, ydim = 1, mc.cores = 2, pb = FALSE, seed = 333)
  corrected2 <- cyCombine(df, xdim = 1, ydim = 1, mc.cores = 2, pb = FALSE, seed = 333)
  corrected3 <- batch_correct(sce, xdim = 1, ydim = 1, mc.cores = 2, pb = FALSE, seed = 333)
  corrected4 <- suppressWarnings(batch_correct(seu, xdim = 1, ydim = 1, mc.cores = 2, pb = FALSE, seed = 333))

  expect_true(all(corrected$label == corrected2$label))
  expect_true(all(corrected$CD19 == corrected2$CD19))
  expect_true(all(corrected$label == corrected3$label))
  expect_true(all(corrected$CD19 == corrected3$CD19))
  expect_true(all(corrected$label == corrected4$label))
  expect_true(sum(corrected$CD19 - SeuratObject::LayerData(corrected4, "data")["CD19",]) < 1)

})
