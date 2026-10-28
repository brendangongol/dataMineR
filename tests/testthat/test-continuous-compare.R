test_that("continuous_compare handles empty and ordinary predictor sets", {
  first <- data.table::data.table(x = 1:5)
  second <- data.table::data.table(x = 4:8)

  empty <- continuous_compare(first, second, "a-b", character())
  expect_s3_class(empty, "data.table")
  expect_equal(nrow(empty), 0)

  result <- continuous_compare(first, second, "a-b", "x")
  expect_equal(result$parameter, "x")
  expect_false(result$skipped)
  expect_true(is.finite(result$CliffsDelta))
  expect_true(result$mann_whitney_u_pvalue >= 0 &&
                result$mann_whitney_u_pvalue <= 1)
})
