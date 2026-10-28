test_that("InitialTablebg rejects invalid weights", {
  dt <- data.table::data.table(group = c(0, 1), value = c(1, 2), weight = c(1, -1))
  expect_error(
    InitialTablebg(dt, GroupBy = "group", weight_col = "weight"),
    "non-negative"
  )
})
