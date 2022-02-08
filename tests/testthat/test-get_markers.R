test_that("get_markers finds markers", {
  df <- tibble::tibble("A" = 1, "id" = 2, "B" = 3)
  expect_equal(get_markers(df), c("A", "B"))
})
