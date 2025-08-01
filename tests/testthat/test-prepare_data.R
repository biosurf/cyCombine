data_dir <- system.file("extdata", package = "cyDefine")
test_that("prepare data functions works", {
  fs <- compile_fcs(data_dir)
  expect_is(fs, "flowSet")
  df <- convert_flowset(fs)
  expect_is(df, "data.frame")
  df <- transform_asinh(df)
  expect_is(df, "data.frame")
  expect_false(any(is.na(df)))
})

test_that("prepare data wrapper works", {
  df <- prepare_data(data_dir, transform = FALSE)
  expect_is(df, "data.frame")
  expect_false(any(is.na(df)))
})
