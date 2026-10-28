test_that("automatic family inference is conservative", {
  expect_identical(dataMineR:::PS_family(c(0, 1, 0)), "binomial")
  expect_identical(dataMineR:::PS_family(c(0, 1, 2)), "gaussian")
  expect_identical(dataMineR:::PS_family(c(0, 1, 2), "poisson"), "poisson")
})

test_that("propensity models report incomplete inputs before fitting", {
  d <- data.frame(exposure = c(0, 1, 0, 1), x = c(1, NA, 3, 4), w = 1)
  expect_error(
    PSFit(d, "exposure", "x", method = "glm", weights = "w"),
    "Problem columns: x \\(1\\)"
  )
})

test_that("compact matching results retain diagnostics and setup identity", {
  result <- data.table::data.table(
    ModelID = "m1", ExposureVariable = "a", IndependentVariable = "a + x",
    MatchingVariable = "x", Method = "glm", Type = "Logistic",
    Ratio = 2, PSCaliperSD = 0.1, CovariateCalipers = "x = 2",
    status = "failed", error = "deliberate failure"
  )
  compact <- dataMineR:::Matching_SelectColumns(result, detail = FALSE)
  expect_true(all(c("status", "error", "ExposureVariable", "IndependentVariable",
                    "Ratio", "PSCaliperSD", "CovariateCalipers") %in% names(compact)))
  expect_identical(compact$error, "deliberate failure")
})

test_that("overview subset rules are validated rather than evaluated", {
  d <- data.frame(p = c(0.01, 0.10, NA))
  expect_equal(dataMineR:::Matching_ApplySubset(d, "p", "<0.05", "rule")$p, 0.01)
  expect_error(
    dataMineR:::Matching_ApplySubset(d, "p", "<0.05; stop('run')", "rule"),
    "one comparison operator"
  )
})

test_that("LassoCompile records dependent failures and continues", {
  d <- data.frame(
    y1 = rep(c(0, 1), 10), y2 = rep(c(1, 0), 10), x = seq_len(20),
    w = 1, cluster = rep(seq_len(5), each = 4), strata = rep(1:2, each = 10)
  )
  result <- suppressMessages(
    LassoCompile(
      d, dependents = c("y1", "y2"), predictors = "x", weights = "w",
      cluster = "cluster", strata = "strata", B = 0,
      excluded_reasons = FALSE, checkpoint = FALSE
    )
  )
  expect_equal(nrow(result$models), 2)
  expect_true(all(result$models$status == "failed"))
  expect_true(all(grepl("At least two nonconstant numeric predictors", result$models$error)))
})

test_that("exposure-specific LASSO roles are not mixed", {
  results <- list(table = data.frame(
    dependent = c("outcome", "exposure_a", "exposure_b"),
    independent = "x", selected = c(TRUE, TRUE, TRUE), stringsAsFactors = FALSE
  ))
  roles <- data.frame(
    exposure = c("exposure_a", "exposure_b"), predictor = "x",
    role = c("confounder (strong)", "near-duplicate of exposure"),
    selected_for_exposure = c("exposure_a", "exposure_b"),
    stringsAsFactors = FALSE
  )
  d <- data.frame(outcome = c(0, 1), exposure_a = c(0, 1), exposure_b = c(1, 0), x = c(1, 2))
  setup <- MatchingSetupFromLasso(
    results, roles, exposure = "exposure_a", DT = d, dependents = "outcome",
    type = "Logistic"
  )
  expect_identical(setup$MatchingVariable, "x")
})

test_that("wlasso defaults no longer fail on a vector condition", {
  d <- data.frame(y = c(0, 1, 0, 1), x1 = 1:4, x2 = 4:1,
                  cluster = 1:4, strata = 1, w = 1)
  expect_error(
    svyVarSel_wlasso_adjusted(
      d, col.y = "y", col.x = c("x1", "x2"), cluster = "cluster",
      strata = "strata", weights = "w", family = "binomial", R = 0
    ),
    "R=0 lower than 1",
    fixed = TRUE
  )
})
