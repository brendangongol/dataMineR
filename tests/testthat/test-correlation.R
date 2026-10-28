test_that("flattenCorrMatrix uses matrix column names and validates dimensions", {
  correlation <- matrix(c(1, 0.5, 0.5, 1), 2, 2,
                        dimnames = list(c("row_a", "row_b"), c("col_a", "col_b")))
  p_values <- matrix(c(0, 0.1, 0.1, 0), 2, 2)

  result <- flattenCorrMatrix(correlation, p_values)

  expect_equal(result$row, "row_a")
  expect_equal(result$column, "col_b")
  expect_error(flattenCorrMatrix(correlation, matrix(1, 1, 1)), "same dimensions")
})
