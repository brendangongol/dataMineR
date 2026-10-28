test_that("ComparisonSetup generates valid combinations without an empty model", {
  expect_equal(ComparisonSetup(3, 2), c("1_2", "1_3", "2_3"))
  expect_equal(length(ComparisonSetup(4, "ALL")), 2^4 - 1)
  expect_false(any(ComparisonSetup(3, "ALL") == "0_0_0"))
  expect_error(ComparisonSetup(2, 3), "between 1 and n")
})

test_that("RegressionCompile does not require a categorical reference for linear models", {
  set.seed(1)
  dt <- data.table::data.table(
    x = rnorm(40),
    y = rnorm(40),
    z = rnorm(40)
  )
  original_names <- names(dt)

  result <- RegressionCompile(
    DT = dt,
    dependents = c("y", "z"),
    independents = c("x", "z"),
    Fpath = tempdir(),
    ReturnSignificant = FALSE,
    removeAdditionalDependents = TRUE,
    type = "Linear",
    WriteFiles = FALSE
  )

  expect_gt(nrow(result), 0)
  expect_true(all(result$independent == "x"))
  expect_equal(names(dt), original_names)
})

test_that("RegressionStabilityR is reproducible when Seed is supplied", {
  set.seed(10)
  dt <- data.table::data.table(
    x = rnorm(60),
    y = rnorm(60)
  )
  dt[, y := 2 * x + y / 10]

  first <- RegressionStabilityR(
    dt, dependents = "y", independents = "x",
    Nboots = 4, Fraction = 0.75, Seed = 42
  )
  second <- RegressionStabilityR(
    dt, dependents = "y", independents = "x",
    Nboots = 4, Fraction = 0.75, Seed = 42
  )

  expect_equal(first, second)
  expect_error(
    RegressionStabilityR(dt, "y", "x", Nboots = 1, Seed = NA_real_),
    "Seed must be NULL or a finite integer"
  )
})
