test_that("ComparisonOverviewBubblePlot accepts compiled results directly", {
  results <- data.table::data.table(
    comparison = c("group_a", "group_a", "group_b", "group_b"),
    parameter = c("x", "y", "x", "y"),
    p_value = c(0.01, 0.20, 0.03, 0.40),
    effect = c(-0.5, 0.1, 0.7, -0.2)
  )

  overview <- ComparisonOverviewBubblePlot(
    ResultsDT = results,
    SizeVariable = "p_value",
    ColorVariable = "effect",
    ClusterMethod = "none",
    ClusterAxes = "none",
    ShowMarginalBars = FALSE
  )

  expect_s3_class(overview$Plot, "ggplot")
  expect_equal(nrow(overview$Table), 2L)
  expect_error(
    ComparisonOverviewBubblePlot(
      ResultsDT = results[, .(comparison, parameter, p_value)],
      SizeVariable = "p_value",
      ColorVariable = "effect"
    ),
    "ResultsDT is missing required columns"
  )
})
