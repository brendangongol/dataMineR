test_that("CategoricalVariableCreate finds only binary non-blank columns", {
  dt <- data.table::data.table(
    binary = c(1, 0, 1, NA),
    multi = c("a", "b", "c", "a"),
    blank = c("", "", NA, "")
  )

  result <- CategoricalVariableCreate(dt, names(dt))

  expect_equal(result$CategoricalVariable, "binary")
  expect_equal(result$Table1Subset, "0")
  expect_equal(result$Table2Subset, "1")
})

test_that("CategoricalVariableCreate validates requested columns", {
  expect_error(
    CategoricalVariableCreate(data.frame(x = 1:3), c("x", "missing")),
    "missing"
  )
})
