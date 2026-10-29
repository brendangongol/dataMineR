#' Write a workbook and confirm it was written
#'
#' Turns write warnings into errors and checks that the file exists with every
#' sheet.
#'
#' @param sheets Named list of data frames.
#' @param file Output path.
#' @param colWidths Passed to \code{openxlsx::write.xlsx}.
#' @return Invisibly, \code{file}.
#' @keywords internal
#' @noRd
Lasso_SaveWorkbook <- function(sheets, file, colWidths = "auto") {
  if (dir.exists(file)) stop("Workbook destination is a directory: ", file)
  withCallingHandlers(openxlsx::write.xlsx(sheets, file=file, colWidths=colWidths),
    warning=function(w) stop(conditionMessage(w), call.=FALSE))
  if (!file.exists(file) || !all(names(sheets) %in% openxlsx::getSheetNames(file)))
    stop("Workbook export verification failed: ", file)
  invisible(file)
}

#' Record a dependent that could not be fitted
#'
#' @param dep,fam Dependent and family.
#' @param predictors Candidate predictors for that dependent.
#' @param error Error message.
#' @return A result shaped like \code{Lasso_OneDependent()} so compilation can
#'   continue with the remaining dependents.
#' @keywords internal
#' @noRd
Lasso_FailedDependent <- function(dep, fam, predictors, error) {
  shown <- if (length(predictors)) predictors else NA_character_
  table <- data.frame(
    dependent = dep, independent = shown, family = fam, selected = FALSE,
    lasso_coef = NA_real_, selected_percent = NA_real_, stability_lower = NA_real_,
    stability_upper = NA_real_, Estimate = NA_real_, `Std. Error` = NA_real_,
    `Pr(>|t|)` = NA_real_, p_fdr = NA_real_, effect_measure = NA_character_,
    effect = NA_real_, effect_lower = NA_real_, effect_upper = NA_real_,
    shrinkage_ratio = NA_real_, shift_in_SE = NA_real_, sign_change = NA,
    refit_reliable = NA, excluded_reason = NA_character_, status = "failed",
    error = error, check.names = FALSE, stringsAsFactors = FALSE)
  model <- data.frame(dependent = dep, family = fam, status = "failed",
                      error = error, notes = paste0("fit: ", error),
                      stringsAsFactors = FALSE)
  list(table = table, model = model, final_model = data.frame(),
       lambda_path = data.frame(), warnings = data.frame(step = "fit", message = error),
       bootstrap_strata = data.frame(), calibration = NULL, fit = NULL,
       selection = NULL, excluded = NULL, bootstrap = NULL, refit = NULL,
       evaluation = NULL, rows = integer(0))
}

#' Weighted LASSO for several dependents
#'
#' Runs the weighted LASSO chain for each dependent and compiles the results
#' like \code{RegressionCompile()}: (1) tune the lambda grid with
#' \code{svyVarSel_wlasso_tune()}; (2) choose lambda with
#' \code{svyVarSel_lambda_select()}; (3) explain each excluded predictor
#' (\code{svyVarSel_excluded_reasons()}); (4) coefficient stability bootstrap
#' (\code{svyVarSel_coef_bootstrap()}); (5) unpenalized refit of the selected
#' predictors (\code{svyVarSel_unpenalized_refit()}); (6) validation on held-out
#' clusters (\code{evaluate_prediction_model()}, binomial dependents only).
#' Include the exposure among the dependents to get both halves of double
#' selection. A dependent is never its own predictor, and each dependent uses
#' its own complete records. A step that fails is recorded in the Model sheet's
#' \code{notes}; the other steps and dependents still run.
#'
#' @param DT Data frame with the dependents, predictors and design columns.
#' @param dependents Outcome(s) and, for double selection, the exposure(s).
#' @param predictors Candidate predictors.
#' @param family Optional named vector of families by dependent, for example
#'   \code{c(LOS = "poisson")}; unspecified 0/1 outcomes use binomial and all
#'   other outcomes use gaussian. Poisson outcomes must be named explicitly.
#' @param weights Name of the weight column.
#' @param cluster Name of the cluster column.
#' @param strata Name of the strata column.
#' @param validation Optional held-out records for step 6.
#' @param method,k,R Replicate-weight cross-validation settings, passed to
#'   \code{svyVarSel_wlasso_adjusted()}.
#' @param lambda_rule \code{"min"} or \code{"1se"}.
#' @param max_rounds Maximum lambda-grid refinements.
#' @param seed Random seed.
#' @param excluded_reasons Logical; run step 3.
#' @param B Bootstrap refits for step 4 (0 skips it).
#' @param B_validation Bootstrap draws for the validation intervals.
#' @param Fpath Optional output folder: one
#'   \code{LASSO_Sub_DataPredictions_<dependent>.xlsx} per dependent and
#'   \code{LASSO_README.xlsx}.
#' @param checkpoint Logical; with \code{Fpath}, save each finished dependent
#'   and reuse it on a rerun when the data and settings are unchanged.
#' @param singleton \code{"fixed"} or \code{"fail"}: handling of
#'   single-cluster strata in the bootstrap.
#' @param exclude_predictors Optional list, named by dependent, of predictors
#'   to leave out of that dependent's LASSO (for example near-duplicates of an
#'   exposure found in a first run). \code{NULL} offers every predictor to every
#'   LASSO.
#' @return Invisibly, a list: \code{table} (one row per dependent x candidate
#'   predictor: selection, LASSO coefficient, bootstrap stability, unpenalized
#'   effect and confidence interval, shrinkage, refit reliability, exclusion
#'   reason), \code{models} (one row per dependent: data, lambda, stability,
#'   refit comparison, validation, notes) and \code{fits} (per dependent: the
#'   fit and the result of each step, and the rows of \code{DT} used).
#' @seealso \code{\link{LassoOverviewPlot}}, \code{\link{MatchingSetupFromLasso}}
#' @examples
#' \dontrun{
#' res <- LassoCompile(training, dependents = c("Death_01", "LOS", "sud_any_01"),
#'   predictors = preds, family = c(LOS = "poisson"), weights = "DISCWT",
#'   cluster = "cluster_id", strata = "strata_id", validation = validation,
#'   Fpath = tempdir())
#' res$models
#' }
#' @export
LassoCompile <- function(DT, dependents, predictors, family = NULL, weights, cluster, strata, validation = NULL,
                         method = "dCV", k = 5, R = 1, lambda_rule = c("min", "1se"), max_rounds = 3, seed = 1621,
                         excluded_reasons = TRUE, B = 200, B_validation = 500, Fpath = NULL, checkpoint = TRUE,
                         singleton = c("fixed", "fail"), exclude_predictors = NULL) {
  singleton <- match.arg(singleton)
  if (is.null(exclude_predictors)) exclude_predictors <- list()
  if (length(exclude_predictors) && (is.null(names(exclude_predictors)) || any(!names(exclude_predictors) %in% dependents)))
    stop("exclude_predictors must be a list named by dependent.")
  if (!is.list(exclude_predictors) || any(vapply(exclude_predictors, function(x) !is.character(x) || anyNA(x) || any(!x %in% predictors), logical(1))))
    stop("Each exclude_predictors entry must contain existing candidate predictor names.")
  lambda_rule <- match.arg(lambda_rule)
  DT <- as.data.frame(DT)
  design_cols <- c(weights, cluster, strata)
  missing_cols <- setdiff(c(dependents, predictors, design_cols), names(DT))
  if (length(missing_cols)) { stop("Columns not in DT: ", paste(missing_cols, collapse = ", ")) }
  if (!is.null(family) && (is.null(names(family)) || !all(names(family) %in% dependents))) {
    stop("family must be named by dependent, e.g. c(LOS = \"poisson\").")
  }
  settings <- list(report_version = 2L, singleton = singleton, exclude_predictors = exclude_predictors, predictors = predictors, weights = weights, cluster = cluster, strata = strata, method = method, k = k, R = R,
                   lambda_rule = lambda_rule, max_rounds = max_rounds, seed = seed, excluded_reasons = excluded_reasons, B = B,
                   B_validation = B_validation,
                   validation = if (is.null(validation)) NULL else Matching_Fingerprint(validation, NULL, list(), c(dependents, predictors, design_cols))$data)
  cp <- Matching_CheckpointDir(Fpath, "LASSO", checkpoint)
  results <- list()
  for (n_dep in seq_along(dependents)) {
    dep <- dependents[n_dep]
    fam <- PS_family(DT[[dep]], if (dep %in% names(family)) family[[dep]] else NULL)
    fp <- Matching_Fingerprint(DT, data.frame(dependent = dep, family = fam), settings, c(dep, predictors, design_cols))
    saved <- Matching_CheckpointLoad(cp, n_dep, dep, fp)
    if (!is.null(saved)) {
      message(sprintf("LASSO %d of %d: %s (from checkpoint)", n_dep, length(dependents), dep))
      results[[dep]] <- saved$rows
      if (!"status" %in% names(results[[dep]]$table)) results[[dep]]$table$status <- "ok"
      if (!"error" %in% names(results[[dep]]$table)) results[[dep]]$table$error <- ""
      if (!"status" %in% names(results[[dep]]$model)) results[[dep]]$model$status <- "ok"
      if (!"error" %in% names(results[[dep]]$model)) results[[dep]]$model$error <- ""
      next
    }
    message(sprintf("LASSO %d of %d: %s (%s)", n_dep, length(dependents), dep, fam))
    dep_predictors <- setdiff(predictors, c(dep, exclude_predictors[[dep]]))
    attempt <- tryCatch(
      Lasso_OneDependent(DT, dep, fam, dep_predictors, weights, cluster, strata, validation,
                         method, k, R, lambda_rule, max_rounds, seed, excluded_reasons, B, B_validation, singleton),
      error = function(e) e)
    if (inherits(attempt, "error")) {
      message("  ", dep, " failed: ", conditionMessage(attempt), "; continuing with the next dependent.")
      results[[dep]] <- Lasso_FailedDependent(dep, fam, dep_predictors, conditionMessage(attempt))
    } else {
      results[[dep]] <- attempt
    }
    results[[dep]]$model$excluded_by_user <- paste(exclude_predictors[[dep]], collapse="; ")
    if (!inherits(attempt, "error"))
      Matching_CheckpointSave(cp, n_dep, dep, fp, results[[dep]], NULL, NULL)
  }
  table <- as.data.frame(data.table::rbindlist(lapply(results, `[[`, "table"), fill = TRUE)); rownames(table) <- NULL
  models <- as.data.frame(data.table::rbindlist(lapply(results, `[[`, "model"), fill = TRUE)); rownames(models) <- NULL
  if (!is.null(Fpath)) {
    dir.create(Fpath, recursive = TRUE, showWarnings = FALSE)
    export_ok <- TRUE
    for (dep in names(results)) {
      sheets <- list(RegressionAnalysis = results[[dep]]$table, Model = results[[dep]]$model,
        FinalModel = results[[dep]]$final_model, LambdaPath = results[[dep]]$lambda_path,
        Warnings = results[[dep]]$warnings, BootstrapStrata = results[[dep]]$bootstrap_strata)
      if (!is.null(results[[dep]]$evaluation)) { sheets$Calibration <- results[[dep]]$calibration }
      tryCatch(Lasso_SaveWorkbook(sheets, file = file.path(Fpath, paste0("LASSO_Sub_DataPredictions_", dep, ".xlsx")), colWidths = "auto"),
               error = function(e) { export_ok <<- FALSE; message(sprintf("  Could not write the file for '%s': %s", dep, conditionMessage(e))) })
    }
    tryCatch(Lasso_WriteReadme(table, models, Fpath), error=function(e) {
      export_ok <<- FALSE; message("Could not write README: ", conditionMessage(e))
    })
    all_fit_ok <- all(vapply(results, function(x) !any(x$model$status == "failed"), logical(1)))
    if (export_ok && all_fit_ok) Matching_CheckpointClear(cp)
    else if (!export_ok) warning("Some exports failed. Checkpoints retained; rerun after resolving the file error.")
    else warning("One or more dependents failed. Successful checkpoints were retained so a rerun can retry only failed dependents.")
  }
  fits <- lapply(results, function(r) r[c("fit", "selection", "excluded", "bootstrap", "refit", "evaluation", "rows", "final_model", "lambda_path", "warnings", "bootstrap_strata")])
  invisible(list(table = table, models = models, fits = fits))
}

#' Run the LASSO chain for one dependent
#'
#' Steps 1 to 6 of \code{LassoCompile()} for one dependent.
#'
#' @param DT Data frame.
#' @param dep Dependent.
#' @param fam Its family.
#' @param predictors Candidate predictors (without the dependent).
#' @param weights,cluster,strata Design columns.
#' @param validation Held-out records, or \code{NULL}.
#' @param method,k,R,lambda_rule,max_rounds,seed,excluded_reasons,B,B_validation,singleton
#'   As in \code{LassoCompile()}.
#' @return A list: \code{table}, \code{model}, \code{final_model},
#'   \code{lambda_path}, \code{warnings}, \code{bootstrap_strata},
#'   \code{calibration}, the step results and \code{rows}.
#' @keywords internal
#' @noRd
Lasso_OneDependent <- function(DT, dep, fam, predictors, weights, cluster, strata, validation, method, k, R, lambda_rule,
                               max_rounds, seed, excluded_reasons, B, B_validation, singleton = "fixed") {
  notes <- character(0)
  warnings <- data.frame(step=character(), message=character())
  capture <- function(name, expr) withCallingHandlers(expr, warning=function(w) {
    warnings <<- unique(rbind(warnings, data.frame(step=name, message=conditionMessage(w))))
    notes <<- unique(c(notes, paste0(name, ": ", conditionMessage(w))))
    invokeRestart("muffleWarning")
  })
  step <- function(name, expr) tryCatch(capture(name, expr), error = function(e) { notes <<- c(notes, sprintf("%s: %s", name, conditionMessage(e))); NULL })
  rows <- which(complete.cases(DT[, c(dep, predictors, weights, cluster, strata)]))
  d <- DT[rows, , drop = FALSE]
  if (!nrow(d)) stop("No complete records remain for ", dep, ".")
  if (any(!vapply(d[, predictors, drop = FALSE], is.numeric, logical(1))))
    stop("LASSO predictors must be numeric or dummy-encoded before fitting.")
  if (any(!is.finite(d[[weights]]) | d[[weights]] <= 0))
    stop("LASSO weights must be finite and positive.")
  constant <- predictors[vapply(predictors, function(v) length(unique(d[[v]])) < 2, logical(1))]
  if (length(constant)) { notes <- c(notes, paste("constant, left out:", paste(constant, collapse = ", "))) }
  preds <- setdiff(predictors, constant)
  if (length(preds) < 2L)
    stop("At least two nonconstant numeric predictors are required by glmnet; ", length(preds), " remain for ", dep, ".")
  fit <- capture("fit", suppressMessages(svyVarSel_wlasso_tune(data = d, col.y = dep, col.x = preds, family = fam, strata = strata, weights = weights,
                                                cluster = cluster, method = method, k = k, R = R, max_rounds = max_rounds, seed = seed,
                                                PrintSummary = FALSE)))
  sel <- capture("lambda selection", svyVarSel_lambda_select(fit, lambda_rule = lambda_rule, PrintSummary = FALSE))
  j <- sel$j
  beta <- fit$model$grid$beta[, j]
  selected <- names(beta)[beta != 0]
  excluded <- if (excluded_reasons && length(selected) < length(preds))
    step("excluded reasons", svyVarSel_excluded_reasons(fit, d, dep, weights, j = j, family = fam, PrintSummary = FALSE)) else NULL
  bootstrap <- if (B > 0 && length(selected))
    step("bootstrap", suppressMessages(svyVarSel_coef_bootstrap(fit, d, dep, weights, cluster, strata, lambda = sel$lambda, B = B,
                                                                seed = seed, family = fam, progress_every = 0, PrintSummary = FALSE, singleton = singleton))) else NULL
  refit <- if (length(selected))
    step("refit", svyVarSel_unpenalized_refit(fit, d, dep, weights, cluster, strata, lambda = sel$lambda, family = fam, PrintSummary = FALSE)) else NULL
  evaluation <- NULL
  if (!is.null(validation) && fam == "binomial") {
    v <- as.data.frame(validation)
    v <- v[complete.cases(v[, c(dep, preds, weights, cluster, strata)]), , drop = FALSE]
    at_j <- fit
    at_j$model$min <- list(a0 = fit$model$grid$a0[j], beta = fit$model$grid$beta[, j, drop = FALSE], df = fit$model$grid$df[j])
    evaluation <- step("validation", evaluate_prediction_model(at_j, v, training_clusters = unique(d[[cluster]]), outcome = dep,
                                                                                weights = weights, cluster = cluster, strata = strata,
                                                                                B = B_validation, seed = seed))
  } else if (!is.null(validation)) {
    notes <- c(notes, sprintf("validation: not available yet for %s dependents", fam))
  }
  tab <- data.frame(dependent = dep, independent = preds, family = fam, selected = preds %in% selected,
                    lasso_coef = unname(beta[preds]), stringsAsFactors = FALSE)
  b <- if (is.null(bootstrap)) NULL else bootstrap$table
  tab$selected_percent <- if (is.null(b)) NA_real_ else b$selected_percent[match(preds, b$variable)]
  tab$stability_lower <- if (is.null(b)) NA_real_ else b$stability_lower[match(preds, b$variable)]
  tab$stability_upper <- if (is.null(b)) NA_real_ else b$stability_upper[match(preds, b$variable)]
  rt <- if (is.null(refit)) NULL else refit$table
  at <- if (is.null(rt)) rep(NA_integer_, length(preds)) else match(preds, rt$variable)
  pick <- function(col) if (is.null(rt)) rep(NA_real_, length(preds)) else rt[[col]][at]
  tab$Estimate <- pick("unpenalized_coef")
  tab$`Std. Error` <- pick("std_error")
  tab$`Pr(>|t|)` <- pick("p_value")
  tab$p_fdr <- p.adjust(tab$`Pr(>|t|)`, "BH")
  tab$effect_measure <- ifelse(tab$selected, switch(fam, binomial = "Odds ratio", poisson = "Ratio of expected counts", gaussian = "Difference"), NA_character_)
  tab$effect <- pick(switch(fam, binomial = "odds_ratio", poisson = "rate_ratio", gaussian = "difference"))
  tab$effect_lower <- pick("lower_95")
  tab$effect_upper <- pick("upper_95")
  tab$shrinkage_ratio <- pick("shrinkage_ratio")
  tab$shift_in_SE <- pick("shift_in_SE")
  tab$sign_change <- if (is.null(rt)) NA else rt$sign_change[at]
  separated <- !is.null(refit) && isTRUE(refit$separation$detected)
  tab$refit_reliable <- ifelse(tab$selected & !is.null(refit), !separated, NA)
  if (separated) {
    notes <- c(notes, sprintf("separation in the unpenalized refit (%d records with fitted probability 0 or 1; extreme coefficients: %s): its effects, CIs and p-values are unreliable; use the LASSO selection only",
                              refit$separation$records_at_bounds,
                              if (length(refit$separation$variables)) paste(refit$separation$variables, collapse = ", ") else "none"))
  }
  es <- if (is.null(excluded)) NULL else excluded$summary
  tab$excluded_reason <- if (is.null(es) || !nrow(es)) NA_character_ else es$reason[match(preds, es$variable)]
  tab$status <- "ok"
  tab$error <- ""
  tab <- tab[order(!tab$selected, -abs(tab$lasso_coef)), ]
  rownames(tab) <- NULL
  ev <- if (is.null(evaluation)) NULL else evaluation$summary
  metric <- function(m, col) if (is.null(ev)) NA_real_ else ev[[col]][ev$metric == m]
  y <- d[[dep]]
  model <- data.frame(dependent = dep, family = fam, records = nrow(d), hospitals = length(unique(d[[cluster]])),
                      events = if (fam == "binomial") sum(y) else NA_real_, outcome_mean = weighted.mean(y, d[[weights]]),
                      candidates = length(preds), n_selected = length(selected), lambda_rule = lambda_rule, lambda = sel$lambda,
                      lambda_position = j, n_lambdas = length(sel$grid), lambda_min_at_grid_end = sel$at_boundary,
                      tuning_rounds = nrow(fit$tuning), cv_error = sel$cv_mean[j],
                      bootstrap_successful = if (is.null(bootstrap)) NA_real_ else bootstrap$successful, bootstrap_B = if (B > 0) B else NA_real_,
                      lp_slope = if (is.null(refit$comparison)) NA_real_ else refit$comparison$lp_slope,
                      lp_correlation = if (is.null(refit$comparison)) NA_real_ else refit$comparison$lp_correlation,
                      median_shrinkage = if (is.null(refit$comparison)) NA_real_ else refit$comparison$median_shrinkage,
                      validation_records = if (is.null(evaluation)) NA_real_ else evaluation$validation_records,
                      AUC = metric("AUC", "estimate"), AUC_lower = metric("AUC", "lower"), AUC_upper = metric("AUC", "upper"),
                      Brier = metric("Brier", "estimate"), Log_loss = metric("Log_loss", "estimate"),
                      Observed_expected = metric("Observed_expected", "estimate"),
                      Observed_expected_lower = metric("Observed_expected", "lower"), Observed_expected_upper = metric("Observed_expected", "upper"),
                      notes = paste(notes, collapse = "; "), stringsAsFactors = FALSE)
  model$status <- "ok"
  model$error <- ""
  cluster_table <- unique(d[, c(cluster, strata), drop=FALSE])
  counts <- table(cluster_table[[strata]])
  bootstrap_strata <- data.frame(stratum=names(counts), hospitals=as.integer(counts))
  bootstrap_strata$records <- vapply(bootstrap_strata$stratum,function(z) sum(d[[strata]]==z),integer(1))
  bootstrap_strata$weight <- vapply(bootstrap_strata$stratum,function(z) sum(d[[weights]][d[[strata]]==z]),numeric(1))
  bootstrap_strata$held_fixed <- bootstrap_strata$hospitals==1
  fixed <- bootstrap_strata$held_fixed
  model$singleton_handling <- singleton
  model$singleton_strata <- sum(fixed)
  model$singleton_records <- sum(bootstrap_strata$records[fixed])
  model$singleton_weight_share <- sum(bootstrap_strata$weight[fixed])/sum(bootstrap_strata$weight)
  model$bootstrap_interpretation <- "Fixed-lambda descriptive stability ranges, not confidence intervals; fixed singleton strata contribute no between-hospital resampling variation."
  model$lambda_min <- sel$grid[sel$j_min]; model$lambda_1se <- sel$grid[sel$j_1se]
  model$intercept <- as.numeric(fit$model$grid$a0[j])
  model$link <- switch(fam,binomial="logit",poisson="log",gaussian="identity")
  model$prediction_formula <- paste(dep,"~",if(length(selected)) paste(selected,collapse=" + ") else "1")
  model$prediction_equation <- paste0(model$link, "(E[", dep, "]) = ", format(model$intercept,digits=16),
    if(length(selected)) paste0(" + ", paste(sprintf("(%.16g * `%s`)",beta[selected],selected),collapse=" + ")) else "")
  final_model <- data.frame(dependent=dep,term=c("(Intercept)",selected),
    penalized_coefficient=c(model$intercept,as.numeric(beta[selected])),lambda=sel$lambda,
    lambda_rule=lambda_rule,family=fam,link=model$link)
  lambda_path <- data.frame(position=seq_along(sel$grid),lambda=sel$grid,cv_error=sel$cv_mean,
    cv_se=sel$cv_se,n_selected=colSums(as.matrix(fit$model$grid$beta)!=0),
    selected=seq_along(sel$grid)==j,minimum_error=seq_along(sel$grid)==sel$j_min,
    one_se=seq_along(sel$grid)==sel$j_1se)
  if(nrow(warnings)) message(dep, ": ", nrow(warnings), " warning(s) recorded in Model notes and Warnings sheet.")
  calibration <- if (is.null(evaluation)) NULL else data.frame(risk_band = rownames(evaluation$calibration), evaluation$calibration, row.names = NULL)
  list(table = tab, model = model, final_model=final_model, lambda_path=lambda_path, warnings=warnings,
       bootstrap_strata=bootstrap_strata, calibration = calibration, fit = fit, selection = sel, excluded = excluded,
       bootstrap = bootstrap, refit = refit, evaluation = evaluation, rows = rows)
}

#' Write LASSO_README.xlsx
#'
#' One row per column of the RegressionAnalysis, Model and other sheets:
#' sheet, column, what it is and how to read it.
#'
#' @param table The compiled LASSO table.
#' @param models The compiled model rows.
#' @param Fpath Output folder.
#' @return Invisibly, the guide rows written.
#' @keywords internal
#' @noRd
Lasso_WriteReadme <- function(table, models, Fpath) {
  d <- function(sheet, column, what, how = "") data.frame(Sheet = sheet, Column = column, `What it is` = what, `How to read it` = how,
                                                          check.names = FALSE, stringsAsFactors = FALSE)
  RA <- "RegressionAnalysis"; M <- "Model"; C <- "Calibration"
  guide <- rbind(
    d(RA, "dependent", "The variable the LASSO predicts: an outcome, or an exposure (what predicts the exposure).", "One file per dependent."),
    d(RA, "independent", "A candidate predictor.", "Every candidate is listed, selected or not."),
    d(RA, "family", "Model type: binomial (0/1), poisson (counts, e.g. days) or gaussian (continuous).", ""),
    d(RA, "selected", "TRUE if the LASSO kept this predictor at the chosen lambda.", "Selected rows come first, largest coefficient first."),
    d(RA, "lasso_coef", "LASSO coefficient at the chosen lambda (log odds for binomial, log rate for poisson, difference for gaussian), shrunk toward zero by the penalty.", "0 = left out. Its size is biased toward zero; use Estimate for effect sizes."),
    d(RA, "selected_percent", "Share (%) of bootstrap refits (hospitals resampled within strata) in which the predictor stayed in the model.", "Near 100 = stable choice; below about 70 = the selection could easily have gone the other way."),
    d(RA, "stability_lower", "2.5th percentile of the LASSO coefficient across the bootstrap refits.", "A stability range, not a confidence interval: it includes 0 when the predictor dropped out of some refits."),
    d(RA, "stability_upper", "97.5th percentile of the LASSO coefficient across the bootstrap refits.", "See stability_lower."),
    d(RA, "Estimate", "Coefficient from the unpenalized survey-weighted refit of the selected predictors (log odds, log rate or difference).", "Blank for predictors not selected."),
    d(RA, "Std. Error", "Standard error of Estimate (survey design-based).", ""),
    d(RA, "Pr(>|t|)", "p-value of Estimate in the unpenalized refit.", "Ignores that the predictors were chosen from the same data, so it is somewhat too small."),
    d(RA, "p_fdr", "Pr(>|t|) adjusted for the number of selected predictors of this dependent (Benjamini-Hochberg).", "Below 0.05 = significant after the adjustment."),
    d(RA, "effect_measure", "Scale of effect: Odds ratio (binomial), Rate ratio (poisson: ratio of expected counts) or Difference (gaussian).", ""),
    d(RA, "effect", "The unpenalized effect on the effect_measure scale (exp(Estimate) for ratios).", "Ratios: above 1 = higher outcome with the predictor (per unit for continuous predictors). Difference: above 0 = higher."),
    d(RA, "effect_lower", "Lower 95% confidence limit of effect.", "An interval that includes 1 (ratios) or 0 (differences) is compatible with no effect."),
    d(RA, "effect_upper", "Upper 95% confidence limit of effect.", "See effect_lower."),
    d(RA, "shrinkage_ratio", "lasso_coef / Estimate: how much the penalty shrank the coefficient.", "1 = no shrinkage; 0.8 = pulled 20% toward zero. Unstable when Estimate is near 0."),
    d(RA, "shift_in_SE", "How far the penalty pulled the coefficient toward zero, in standard errors of Estimate.", "Below about 0.2 = negligible next to its own uncertainty; negative = pushed away from zero."),
    d(RA, "sign_change", "TRUE if lasso_coef and Estimate point in opposite directions.", "Should be FALSE; TRUE means the predictor's direction depends on the penalty."),
    d(RA, "refit_reliable", "FALSE when the unpenalized refit of this dependent shows separation: some predictors (almost) perfectly predict a 0/1 dependent, or a coefficient exceeds 10 on the log scale. Blank for predictors not selected.",
      "If FALSE, Estimate, Std. Error, Pr(>|t|), p_fdr, effect and its limits are unreliable for every row of this dependent (the Model sheet notes name the predictors involved). The LASSO selection itself is still valid."),
    d(RA, "excluded_reason", "For predictors not selected: the first reason that applies (overlaps with named partners, small group, reproduced by the other predictors, no association, or explained by the other predictors).", "A description, not proof of why the LASSO left it out."),
    d(M, "dependent", "The variable the LASSO predicts.", ""),
    d(M, "family", "Model type: binomial, poisson or gaussian.", ""),
    d(M, "records", "Records used: complete for the dependent, predictors and survey design.", ""),
    d(M, "hospitals", "Hospitals (clusters) in those records.", ""),
    d(M, "events", "Binomial only: records with the outcome (unweighted).", "Fewer than about 10 per selected predictor makes the model unstable."),
    d(M, "outcome_mean", "Weighted mean of the dependent (the event rate for binomial).", ""),
    d(M, "candidates", "Candidate predictors (constant ones are left out and listed in notes).", ""),
    d(M, "n_selected", "Predictors kept at the chosen lambda.", ""),
    d(M, "lambda_rule", "min = lambda with the lowest cross-validated error; 1se = largest lambda within one standard error of it (fewer predictors).", ""),
    d(M, "lambda", "The chosen penalty.", ""),
    d(M, "lambda_position", "Position of lambda in the grid (1 = largest penalty).", ""),
    d(M, "n_lambdas", "Number of lambdas in the final grid.", ""),
    d(M, "lambda_min_at_grid_end", "TRUE if the lowest CV error is still at the smallest lambda tried.", "Should be FALSE; TRUE means the minimum was not reached: raise max_rounds."),
    d(M, "tuning_rounds", "Fits run while refining the lambda grid.", ""),
    d(M, "cv_error", "Cross-validated prediction error at lambda (deviance-type loss: log loss for binomial, Poisson deviance for poisson, squared error for gaussian).", "Compare between settings for the same dependent only."),
    d(M, "bootstrap_successful", "Bootstrap refits that succeeded.", "Should equal bootstrap_B; the step stops below 90%."),
    d(M, "bootstrap_B", "Bootstrap refits attempted.", ""),
    d(M, "lp_slope", "Slope of the unpenalized refit's linear predictor on the LASSO's.", "Above 1 = LASSO predictions are flatter (shrunk)."),
    d(M, "lp_correlation", "Correlation of the two linear predictors.", "Near 1 = both rank patients the same way."),
    d(M, "median_shrinkage", "Median shrinkage_ratio over the selected predictors.", ""),
    d(M, "validation_records", "Held-out records used for validation.", ""),
    d(M, "AUC", "Validation: area under the ROC curve (binomial only).", "0.5 = chance; 1 = perfect ranking."),
    d(M, "AUC_lower", "Lower limit of the AUC bootstrap interval (hospitals resampled).", ""),
    d(M, "AUC_upper", "Upper limit of the AUC bootstrap interval.", ""),
    d(M, "Brier", "Validation: weighted mean squared error of the predicted probabilities.", "Smaller is better."),
    d(M, "Log_loss", "Validation: weighted log loss.", "Smaller is better."),
    d(M, "Observed_expected", "Validation: weighted observed events / predicted events.", "1 = right overall level; above 1 = underprediction."),
    d(M, "Observed_expected_lower", "Lower limit of the observed/expected interval.", ""),
    d(M, "Observed_expected_upper", "Upper limit of the observed/expected interval.", ""),
    d(M, "notes", "Steps that failed or were skipped, predictors left out as constant, and separation in the unpenalized refit.", "Blank means no recorded notes; it does not establish validity of the model or causal interpretation."),
    d(C, "risk_band", "Validation records grouped by predicted risk (tenths).", ""),
    d(C, "records", "Records in the band.", ""),
    d(C, "predicted", "Weighted mean predicted risk in the band.", ""),
    d(C, "observed", "Weighted observed event rate in the band.", "Close to predicted = well calibrated."))
  guide <- rbind(guide,
    d("FinalModel", "penalized_coefficient", "Selected LASSO coefficients in original predictor units, including the intercept; these reproduce the penalized prediction model."),
    d("LambdaPath", "lambda", "Every fitted penalty with CV error, CV SE, model size and flags for selected/minimum/1-SE lambda."),
    d("Warnings", "message", "Captured warning text and analysis step; review even if fitting completed."),
    d("BootstrapStrata", "held_fixed", "TRUE identifies single-hospital strata; fixed bootstrap sampling cannot represent their between-hospital variation."),
    d(M, "prediction_equation", "Exact fitted penalized linear predictor at the selected lambda; inverse-logit for binomial, exp for Poisson, identity for Gaussian."),
    d(M, "singleton_weight_share", "Fraction of training sampling weight in single-hospital strata."),
    d(M, "lambda_min", "Lambda with the lowest CV error."),
    d(M, "lambda_1se", "Largest lambda within one CV standard error of the minimum."),
    d(M, "excluded_by_user", "Predictors explicitly excluded for this dependent before fitting."))
  present <- (guide$Sheet == RA & guide$Column %in% names(table)) | (guide$Sheet == M & guide$Column %in% names(models)) | guide$Sheet %in% c(C, "FinalModel", "LambdaPath", "Warnings", "BootstrapStrata")
  Lasso_SaveWorkbook(list(Columns = guide[present, ]), file = file.path(Fpath, "LASSO_README.xlsx"), colWidths = "auto")
  invisible(guide[present, ])
}

#' Overview of LassoCompile results in two or three panes
#'
#' One row per covariate, so each role can be traced to its source. An
#' exposure that is also a predictor of the other dependents is a row too,
#' pinned with the other exposures in an \code{Exposures} group at the top, in
#' the order of \code{exposures}, so the exposure-by-exposure correlations form
#' one block in pane 1. Pane 1
#' (from the data): the covariate's weighted correlation with each exposure,
#' from \code{MatchingAssumptionCheck()}. Pane 2 (from \code{LassoCompile()}): a
#' bubble where a dependent's LASSO kept the covariate; colour is the
#' standardized coefficient and size the bootstrap selection percentage. With
#' \code{ShowRefitMarks = TRUE}, an outline marks a significant unpenalized
#' refit and an X marks one made unreliable by separation.
#' Optional pane 3 (\code{ShowRoles = TRUE}): the role and the reason. Without
#' it the role groups label the rows on the right of pane 2. Roles, in order: near-duplicate of the
#' exposure (flagged alone); near-duplicate with other codes (flagged by the
#' joint rule); confounder, strong or weak (kept by an exposure LASSO and at
#' least one outcome LASSO); outcome only; exposure only (a possible
#' instrument). Roles are statistical: a mediator can look like a confounder,
#' so check timing.
#'
#' @param results A \code{LassoCompile()} result.
#' @param data The data given to \code{LassoCompile()}.
#' @param exposures Dependents that are exposures; \code{NULL} leaves out pane
#'   1 and the confounder roles.
#' @param weights Name of the weight column, or \code{NULL}.
#' @param SigColumn Results column used for the optional unpenalized-refit
#'   outline.
#' @param SigThreshold Threshold for \code{SigColumn}.
#' @param ShowRefitMarks Logical; add the unpenalized-refit outline and X marks.
#'   The default is \code{FALSE}.
#' @param WeakThreshold Absolute standardized exposure-LASSO coefficient below
#'   which a confounder is weak.
#' @param phi_threshold,joint_min_cor,share_threshold Passed to
#'   \code{MatchingAssumptionCheck()}.
#' @param ReferenceLine Position of pane 1's dashed lines (display only;
#'   \code{NULL} for none).
#' @param SizeBreaks Bootstrap percentages shown in the size legend.
#' @param ShowNotSelected Logical; also show covariates no LASSO kept.
#' @param ColorLimitQuantile Quantile of the absolute coefficient that caps the
#'   colour scale.
#' @param ClusterMethod Row order within each role: \code{"correlation"},
#'   \code{"signed"} or \code{"none"}.
#' @param ShowRoles Logical; add pane 3, the role and reason for every row.
#' @param Fpath Optional folder to save the plot to (PDF).
#' @param FileName File name used with \code{Fpath}.
#' @param PrintPlot Logical; print the plot.
#' @return Invisibly, a list: \code{plot}, \code{roles} (one overall row per
#'   covariate), \code{roles_by_exposure} (one row per exposure and covariate,
#'   used for exposure-specific matching recommendations), \code{table}
#'   (the plotted cells; \code{role} is the plotted row group, with the
#'   exposures in \code{"Exposures"}), \code{correlations} (every covariate,
#'   exposures included, against every exposure, strongest first: weighted
#'   \code{correlation}, \code{exposed_share_with} (the share of patients with
#'   the covariate who are exposed), the near-duplicate flags and
#'   \code{joint_r}; pane 1 prints no numbers) and \code{exposure_correlations}
#'   (the exposure-by-exposure rows of that table: \code{correlation},
#'   \code{other_share_exposed}, the near-duplicate flags and \code{joint_r};
#'   \code{NULL} without \code{exposures}). The function writes nothing except
#'   the optional plot file.
#' @examples
#' \dontrun{
#' ov <- LassoOverviewPlot(res, training, exposures = "sud_any_01",
#'   weights = "DISCWT")
#' ov$roles
#' }
#' @export
LassoOverviewPlot <- function(results, data, exposures = NULL, weights = NULL, SigColumn = "p_fdr", SigThreshold = 0.05,
                              ShowRefitMarks = FALSE,
                              WeakThreshold = 0.1, phi_threshold = 0.8, joint_min_cor = 0.3, share_threshold = 0.01,
                              ReferenceLine = phi_threshold, SizeBreaks = c(25, 50, 75, 100),
                              ShowNotSelected = FALSE, ColorLimitQuantile = 0.95,
                              ClusterMethod = c("correlation", "signed", "none"), ShowRoles = FALSE, Fpath = NULL,
                              FileName = "LASSO_SelectionOverview.pdf", PrintPlot = TRUE) {
  ClusterMethod <- match.arg(ClusterMethod)
  tab <- as.data.frame(results$table)
  data <- as.data.frame(data)
  if (!is.logical(ShowRefitMarks) || length(ShowRefitMarks) != 1L || is.na(ShowRefitMarks))
    stop("ShowRefitMarks must be TRUE or FALSE.")
  dependents <- unique(tab$dependent)
  if (length(setdiff(exposures, dependents))) { stop("exposures must be dependents of the LassoCompile results: ", paste(setdiff(exposures, dependents), collapse = ", ")) }
  if (ShowRefitMarks && !SigColumn %in% names(tab)) { stop("SigColumn '", SigColumn, "' is not a column of results$table.") }
  outcomes <- setdiff(dependents, exposures)
  tab$std_coef <- NA_real_
  for (dep in dependents) {
    at <- which(tab$dependent == dep)
    rows <- results$fits[[dep]]$rows
    w <- if (is.null(weights)) rep(1, length(rows)) else data[[weights]][rows]
    sds <- vapply(tab$independent[at], function(v) { x <- data[[v]][rows]; m <- weighted.mean(x, w); sqrt(sum(w * (x - m)^2) / sum(w)) }, numeric(1))
    tab$std_coef[at] <- tab$lasso_coef[at] * sds
  }
  tab$size <- ifelse(is.na(tab$selected_percent), 100, tab$selected_percent)
  # An exposure that is also a predictor stays a row, so its link to the other exposures shows (an exposure is never its own covariate).
  covariates <- unique(tab$independent)
  link <- do.call(rbind, lapply(exposures, function(ex) {
    rows <- results$fits[[ex]]$rows
    cv <- setdiff(covariates, ex)   # every covariate, including those left out of the exposure LASSO (exclude_predictors)
    chk <- MatchingAssumptionCheck(data[rows, , drop = FALSE], exposure = ex, covariates = cv, weights = weights, phi_threshold = phi_threshold,
                                   joint_min_cor = joint_min_cor, share_threshold = share_threshold, PrintSummary = FALSE)$table
    data.frame(exposure = ex, covariate = chk$variable, correlation = chk$correlation, exposed_share_with = chk$exposed_share_with,
               duplicate = chk$near_duplicate | chk$joint_duplicate, near_duplicate = chk$near_duplicate, joint_duplicate = chk$joint_duplicate,
               joint_r = suppressWarnings(as.numeric(sub(".*r of any = ([-0-9.]+).*", "\\1", chk$action))), stringsAsFactors = FALSE)
  }))
  chosen <- tab[tab$selected, ]
  sel_by <- function(p, deps) intersect(deps, chosen$dependent[chosen$independent == p])
  excluded_from <- function(p) { others <- setdiff(exposures, p)
    others[!vapply(others, function(ex) p %in% tab$independent[tab$dependent == ex], logical(1))] }
  exp_coef <- vapply(covariates, function(p) { s <- chosen[chosen$independent == p & chosen$dependent %in% exposures, ]
                                               if (nrow(s)) s$std_coef[which.max(abs(s$std_coef))] else NA_real_ }, numeric(1))
  roles <- do.call(rbind, lapply(covariates, function(p) {
    for_exp <- sel_by(p, exposures); for_out <- sel_by(p, outcomes)
    lk <- if (is.null(link)) NULL else link[link$covariate == p, , drop = FALSE]
    dup_alone <- !is.null(lk) && any(lk$near_duplicate)
    dup_joint <- !is.null(lk) && any(lk$joint_duplicate)
    ec <- exp_coef[[p]]
    role <- if (dup_alone) "near-duplicate of exposure" else if (dup_joint) "near-duplicate with other codes" else if (length(for_exp) && length(for_out)) {
      if (abs(ec) >= WeakThreshold) "confounder (strong)" else "confounder (weak)"
    } else if (length(for_out)) { if (is.null(exposures)) "selected" else "outcome only" } else if (length(for_exp)) "exposure only" else "not selected"
    outs <- sprintf("outcomes %d of %d", length(for_out), length(outcomes))
    reason <- switch(role,
      "near-duplicate of exposure" = , "near-duplicate with other codes" = {
        d1 <- lk[lk$duplicate, , drop = FALSE][1, ]
        partners <- setdiff(link$covariate[link$exposure == d1$exposure & link$joint_duplicate], p)
        paste0(sprintf("r = %.2f", d1$correlation), if (d1$near_duplicate) sprintf(" (>= %g)", phi_threshold) else sprintf(" alone (< %g)", phi_threshold),
               if (!is.na(d1$exposed_share_with)) sprintf("; %.0f%% with it are exposed", 100 * d1$exposed_share_with) else "",
               if (d1$joint_duplicate && length(partners)) sprintf("; with %s: r = %.2f", paste(partners, collapse = ", "), d1$joint_r) else "")
      },
      "confounder (strong)" = , "confounder (weak)" = sprintf("exposure LASSO %+.2f (%s |%.2f|); %s", ec, if (abs(ec) >= WeakThreshold) ">=" else "<", WeakThreshold, outs),
      "outcome only" = sprintf("%s; %s", if (length(excluded_from(p))) "left out of the exposure LASSO (exclude_predictors)" else "not selected for the exposure", outs),
      "exposure only" = sprintf("exposure LASSO %+.2f; no outcome: possible instrument", ec),
      "selected" = outs,
      "not selected by any LASSO")
    data.frame(predictor = p, role = role, reason = reason,
               correlation_with_exposure = if (is.null(lk) || !nrow(lk)) NA_real_ else lk$correlation[which.max(abs(lk$correlation))],
               exposed_share_with = if (is.null(lk) || all(is.na(lk$exposed_share_with))) NA_real_ else max(lk$exposed_share_with, na.rm = TRUE),
               exposure_lasso_std_coef = ec, selected_for_exposure = paste(for_exp, collapse = ", "),
               selected_for_outcomes = paste(for_out, collapse = ", "), n_outcomes = length(for_out), stringsAsFactors = FALSE)
  }))
  role_levels <- c("near-duplicate of exposure", "near-duplicate with other codes", "confounder (strong)", "confounder (weak)", "outcome only", "exposure only", "selected", "not selected")
  roles <- roles[order(match(roles$role, role_levels), -roles$n_outcomes), ]
  rownames(roles) <- NULL
  roles_by_exposure <- if (length(exposures)) do.call(rbind, lapply(exposures, function(ex) {
    do.call(rbind, lapply(setdiff(covariates, ex), function(p) {
      for_exp <- p %in% chosen$independent[chosen$dependent == ex]
      for_out <- sel_by(p, outcomes)
      lk <- link[link$exposure == ex & link$covariate == p, , drop = FALSE]
      dup_alone <- nrow(lk) && any(lk$near_duplicate)
      dup_joint <- nrow(lk) && any(lk$joint_duplicate)
      ec_row <- chosen[chosen$dependent == ex & chosen$independent == p, , drop = FALSE]
      ec <- if (nrow(ec_row)) ec_row$std_coef[1] else NA_real_
      role <- if (dup_alone) "near-duplicate of exposure" else if (dup_joint) "near-duplicate with other codes" else
        if (for_exp && length(for_out)) {
          if (is.finite(ec) && abs(ec) >= WeakThreshold) "confounder (strong)" else "confounder (weak)"
        } else if (length(for_out)) "outcome only" else if (for_exp) "exposure only" else "not selected"
      data.frame(exposure = ex, predictor = p, role = role,
                 correlation_with_exposure = if (nrow(lk)) lk$correlation[1] else NA_real_,
                 exposed_share_with = if (nrow(lk)) lk$exposed_share_with[1] else NA_real_,
                 exposure_lasso_std_coef = ec,
                 selected_for_exposure = if (for_exp) ex else "",
                 selected_for_outcomes = paste(for_out, collapse = ", "),
                 n_outcomes = length(for_out), stringsAsFactors = FALSE)
    }))
  })) else data.frame()
  # Exposures are always shown: their correlations with each other are the point of the Exposures group.
  keep <- if (ShowNotSelected) roles$predictor else roles$predictor[roles$role != "not selected" | roles$predictor %in% exposures]
  if (!length(keep)) { stop("No covariate was selected for any dependent.") }
  mat <- matrix(vapply(dependents, function(dep) { s <- tab[tab$dependent == dep, ]; s$std_coef[match(keep, s$independent)] }, numeric(length(keep))),
                nrow = length(keep), dimnames = list(keep, dependents))
  mat[is.na(mat)] <- 0
  order_rows <- function(v) {
    if (ClusterMethod == "none" || length(v) < 3) return(v)
    m <- mat[v, , drop = FALSE]
    d <- if (ClusterMethod == "signed") dist(sign(m), method = "manhattan") else {
      cc <- suppressWarnings(cor(t(m))); cc[!is.finite(cc)] <- 0; as.dist(1 - cc) }
    v[stats::hclust(d, method = "average")$order]
  }
  # Exposures form their own group at the top, in the order of `exposures`, so pane 1's columns line up with them.
  exposure_rows <- intersect(exposures, keep)
  y_lv <- c(exposure_rows, unlist(lapply(role_levels, function(r) order_rows(setdiff(intersect(roles$predictor[roles$role == r], keep), exposure_rows)))))
  role_group <- ifelse(roles$predictor %in% exposure_rows, "Exposures", roles$role)
  group_levels <- c("Exposures", role_levels)
  used_roles <- group_levels[group_levels %in% role_group[roles$predictor %in% keep]]
  as_rows <- function(df, col) {
    df$role <- factor(role_group[match(df[[col]], roles$predictor)], levels = used_roles)
    df[[col]] <- factor(df[[col]], levels = rev(y_lv))
    df
  }
  pane_theme <- ggplot2::theme_minimal() +
    ggplot2::theme(strip.text.y = ggplot2::element_blank(), panel.border = ggplot2::element_rect(fill = NA, colour = "grey80"),
                   plot.title = ggplot2::element_text(size = 10, face = "bold"), plot.subtitle = ggplot2::element_text(size = 8))
  cells <- as_rows(tab[tab$independent %in% keep, ], "independent")
  cells$.column <- factor(ifelse(cells$dependent %in% exposures, "Exposure", "Outcomes"), levels = c("Exposure", "Outcomes"))
  cells$dependent <- factor(cells$dependent, levels = c(exposures, outcomes))
  cells$reliable <- is.na(cells$refit_reliable) | cells$refit_reliable
  cells$sig <- if (ShowRefitMarks) cells$selected & cells$reliable & !is.na(cells[[SigColumn]]) & cells[[SigColumn]] < SigThreshold else FALSE
  shown <- cells[cells$selected, ]
  mark_sig <- sprintf("%s < %g", SigColumn, SigThreshold)
  mark_unreliable <- "Unreliable (separation)"
  left_out <- do.call(rbind, lapply(exposures, function(ex) {
    v <- setdiff(keep, c(ex, tab$independent[tab$dependent == ex]))
    if (length(v)) data.frame(dependent = ex, independent = v, stringsAsFactors = FALSE) else NULL
  }))
  if (is.null(left_out)) left_out <- data.frame(dependent = character(0), independent = character(0))
  left_out <- as_rows(left_out, "independent")
  left_out$.column <- factor(rep("Exposure", nrow(left_out)), levels = c("Exposure", "Outcomes"))
  left_out$dependent <- factor(left_out$dependent, levels = c(exposures, outcomes))
  lim <- stats::quantile(abs(shown$std_coef), ColorLimitQuantile, na.rm = TRUE)
  if (!is.finite(lim) || lim == 0) lim <- max(abs(shown$std_coef), 1e-6, na.rm = TRUE)
  p_sel <- ggplot2::ggplot(cells, ggplot2::aes(x = dependent, y = independent)) +
    ggplot2::geom_point(colour = "grey90", size = 1, show.legend = FALSE) +
    ggplot2::geom_point(data = shown, ggplot2::aes(size = size, fill = std_coef), shape = 21, colour = "grey55", stroke = 0.3) +
    ggplot2::geom_text(data = left_out, ggplot2::aes(label = "excl."), size = 2.3, colour = "grey35", show.legend = FALSE) +
    ggplot2::scale_fill_distiller(palette = "RdBu", direction = -1, limits = c(-lim, lim), oob = scales::squish,
                                  name = "Colour: standardized\nLASSO coefficient\n(red raises, blue lowers)") +
    ggplot2::scale_radius(range = c(0.8, 7), limits = c(0, 100), breaks = SizeBreaks, labels = paste0(SizeBreaks, "%"),
                          name = if (all(is.na(tab$selected_percent))) "Size: bootstrap skipped\n(all shown at 100%)" else "Size: % of bootstrap\nrefits keeping it",
                          guide = ggplot2::guide_legend(override.aes = list(shape = 21, fill = "grey75", colour = "grey40"))) +
    ggplot2::facet_grid(ggplot2::vars(role), ggplot2::vars(.column), scales = "free", space = "free",
                        labeller = ggplot2::labeller(role = ggplot2::label_wrap_gen(20))) +
    ggplot2::labs(x = NULL, y = NULL, title = "2. LASSO selection",
                  subtitle = "Bubble = kept by the LASSO for that column.\nexcl. = left out of that LASSO by the user") +
    pane_theme + ggplot2::theme(axis.text.y = if (is.null(exposures)) ggplot2::element_text() else ggplot2::element_blank(),
                                axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  if (ShowRefitMarks) {
    p_sel <- p_sel +
      ggplot2::geom_point(data = shown[shown$sig, ], ggplot2::aes(size = size, shape = mark_sig),
                          fill = NA, colour = "black", stroke = 1.1,
                          show.legend = c(size = FALSE, shape = TRUE)) +
      ggplot2::geom_point(data = shown[!shown$reliable, ], ggplot2::aes(shape = mark_unreliable),
                          colour = "black", size = 2,
                          show.legend = c(shape = TRUE, size = FALSE, fill = FALSE)) +
      ggplot2::scale_shape_manual(name = "Unpenalized refit",
                                  values = stats::setNames(c(21, 4), c(mark_sig, mark_unreliable)),
                                  guide = ggplot2::guide_legend(override.aes = list(size = 4, fill = NA, stroke = 1.1)))
  }
  # Without pane 3 the row groups are named on the right edge of pane 2.
  if (!ShowRoles) p_sel <- p_sel + ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0, face = "bold", hjust = 0))
  why <- roles[roles$predictor %in% keep, ]
  why$colour_role <- why$role
  # An exposure row sits in the Exposures group, so its role moves into the text.
  why$reason <- ifelse(why$predictor %in% exposure_rows, paste0(why$role, ": ", why$reason), why$reason)
  why <- as_rows(why, "predictor")
  why$.column <- "Role and reason"
  role_colours <- c("near-duplicate of exposure" = "#7F2704", "near-duplicate with other codes" = "#D94801", "confounder (strong)" = "#B2182B", "confounder (weak)" = "#F4A582",
                    "outcome only" = "#4393C3", "exposure only" = "#762A83", "selected" = "#4393C3", "not selected" = "grey60")
  p_why <- ggplot2::ggplot(why, ggplot2::aes(x = 0, y = predictor)) +
    ggplot2::geom_point(ggplot2::aes(colour = colour_role), shape = 15, size = 3, show.legend = FALSE) +
    ggplot2::geom_text(ggplot2::aes(x = 0.04, label = reason), hjust = 0, size = 2.6) +
    ggplot2::scale_colour_manual(values = role_colours) +
    ggplot2::scale_x_continuous(limits = c(-0.02, 1), expand = c(0, 0)) +
    ggplot2::facet_grid(ggplot2::vars(role), ggplot2::vars(.column), scales = "free_y", space = "free_y",
                        labeller = ggplot2::labeller(role = ggplot2::label_wrap_gen(20))) +
    ggplot2::labs(x = NULL, y = NULL, title = "3. Role", subtitle = "Assigned from panes 1 and 2") +
    pane_theme + ggplot2::theme(axis.text = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(),
                                strip.text.y = ggplot2::element_text(angle = 0, face = "bold", hjust = 0))
  if (is.null(exposures)) {
    p <- (if (ShowRoles) patchwork::wrap_plots(p_sel, p_why, widths = c(length(dependents), 4)) else p_sel) +
      patchwork::plot_annotation(title = "LASSO selection overview")
  } else {
    lk <- as_rows(link[link$covariate %in% keep, ], "covariate")
    lk$.column <- paste("Correlation with", lk$exposure)
    p_link <- ggplot2::ggplot(lk, ggplot2::aes(x = correlation, y = covariate)) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey60") +
      { if (length(ReferenceLine)) ggplot2::geom_vline(xintercept = c(-ReferenceLine, ReferenceLine), linetype = 2, colour = "grey40") } +
      ggplot2::geom_col(ggplot2::aes(fill = ifelse(near_duplicate, "alone", ifelse(joint_duplicate, "joint", "none"))), width = 0.6, show.legend = FALSE) +
      ggplot2::scale_fill_manual(values = c(none = "grey65", alone = "#7F2704", joint = "#D94801")) +
      ggplot2::scale_x_continuous(limits = c(-1, 1), breaks = sort(unique(c(-1, 0, 1, -ReferenceLine, ReferenceLine)))) +
      ggplot2::facet_grid(ggplot2::vars(role), ggplot2::vars(.column), scales = "free_y", space = "free_y") +
      ggplot2::labs(x = "weighted correlation (r)", y = NULL, title = "1. Link to the exposure (data)",
                    subtitle = paste0("Bar: r of the covariate with the exposure.\n",
                                      if (length(ReferenceLine)) sprintf("Dashed: |r| = %g (flag at |r| >= %g).\n", ReferenceLine, phi_threshold) else "",
                                      "Values: $correlations")) +
      pane_theme
    p <- (if (ShowRoles) patchwork::wrap_plots(p_link, p_sel, p_why, widths = c(max(3, 2 * length(exposures)), length(dependents), 6.5))
          else patchwork::wrap_plots(p_link, p_sel, widths = c(max(3, 2 * length(exposures)), length(dependents) + 1))) +
      patchwork::plot_annotation(title = "LASSO selection overview")
  }
  p <- p + patchwork::plot_layout(guides = "collect")
  if (!is.null(Fpath)) {
    dir.create(Fpath, recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(file.path(Fpath, FileName), plot = p, width = max(14, 0.9 * length(dependents) + 13), height = max(6, 0.22 * length(keep) + 3.5), limitsize = FALSE,
                    device = if (capabilities("cairo")) grDevices::cairo_pdf else "pdf")
  }
  if (PrintPlot) print(p)
  # Exposure-by-exposure rows of the data link (pane 1), directional: `other_share_exposed` is the share of patients
  # with the other exposure who are exposed. Returned for the caller to save; nothing is written here.
  exposure_correlations <- if (is.null(link)) NULL else {
    pairs <- link[link$covariate %in% exposures, c("exposure", "covariate", "correlation", "exposed_share_with", "near_duplicate", "joint_duplicate", "joint_r")]
    names(pairs) <- c("exposure", "other_exposure", "correlation", "other_share_exposed", "near_duplicate", "joint_duplicate", "joint_r")
    pairs <- pairs[order(match(pairs$exposure, exposures), match(pairs$other_exposure, exposures)), ]
    rownames(pairs) <- NULL
    pairs
  }
  # Every covariate-by-exposure correlation behind pane 1 (the plot prints no numbers).
  correlations <- if (is.null(link)) NULL else {
    out <- link[, c("exposure", "covariate", "correlation", "exposed_share_with", "near_duplicate", "joint_duplicate", "joint_r")]
    out <- out[order(match(out$exposure, exposures), -abs(out$correlation)), ]
    rownames(out) <- NULL
    out
  }
  invisible(list(plot = p, roles = roles, roles_by_exposure = roles_by_exposure, table = cells,
                 correlations = correlations, exposure_correlations = exposure_correlations))
}

#' Recommended matching models from LASSO double selection
#'
#' One setup row per outcome for the exposure, in the format of
#' \code{MatchingRegressionSetup()} (pass the result as its \code{recommended}
#' argument). \code{MatchingVariable} holds the confounders from
#' \code{LassoOverviewPlot()} and is the same for every outcome;
#' \code{IndependentVariable} holds the exposure plus, with
#' \code{outcome_adjust = TRUE}, the other covariates that outcome's LASSO kept.
#' Near-duplicates of the exposure, the exposure itself and \code{exclude} are
#' always left out.
#'
#' @param results A \code{LassoCompile()} result.
#' @param roles \code{LassoOverviewPlot()$roles}; when several exposures were
#'   reviewed, use \code{LassoOverviewPlot()$roles_by_exposure} so duplicate and
#'   strength classifications are specific to the requested exposure.
#' @param exposure The exposure; it must be a dependent of \code{results}.
#' @param DT The data used for matching (to resolve each outcome's Type).
#' @param dependents Outcomes; defaults to every LASSO dependent except the
#'   exposure.
#' @param method \code{"glm"}, \code{"svyglm"} or \code{"lm"}.
#' @param type Outcome model type, as in \code{MatchingRegressionSetup()}.
#' @param include_weak Logical; also match on weak confounders.
#' @param outcome_adjust Logical; add the outcome-only covariates to the
#'   outcome model.
#' @param exclude Covariates never used, for example ones measured after the
#'   exposure.
#' @param notes Logical; add a \code{Notes} column explaining each row.
#' @return A data frame of setup rows with ModelIDs \code{LASSO_1},
#'   \code{LASSO_2}, ...
#' @examples
#' \dontrun{
#' rec <- MatchingSetupFromLasso(res, ov$roles, exposure = "sud_any_01",
#'   DT = model_dt, dependents = c("Death_01", "LOS"))
#' }
#' @export
MatchingSetupFromLasso <- function(results, roles, exposure, DT, dependents = NULL, method = "glm", type = NULL,
                                   include_weak = TRUE, outcome_adjust = TRUE, exclude = NULL, notes = FALSE) {
  tab <- as.data.frame(results$table)
  DT <- as.data.frame(DT)
  if (length(exposure) != 1L || !exposure %in% tab$dependent) stop("exposure must be one dependent of the LassoCompile results.")
  if ("exposure" %in% names(roles)) {
    roles <- roles[roles$exposure == exposure, , drop = FALSE]
    if (!nrow(roles)) stop("roles has no rows for exposure ", exposure, ".")
  }
  if (is.null(dependents)) dependents <- setdiff(unique(tab$dependent), exposure)
  missing_deps <- setdiff(dependents, unique(tab$dependent))
  if (length(missing_deps)) stop("No LASSO results for: ", paste(missing_deps, collapse = ", "), ". Add them to LassoCompile's dependents.")
  if (!method %in% c("glm", "svyglm", "lm")) stop("method must be glm, svyglm or lm.")
  in_list <- function(x, value) vapply(strsplit(x, ", ", fixed = TRUE), function(s) value %in% s, logical(1))
  wanted <- if (include_weak) c("confounder (strong)", "confounder (weak)") else "confounder (strong)"
  conf <- roles[roles$role %in% wanted & in_list(roles$selected_for_exposure, exposure) & !roles$predictor %in% exclude, , drop = FALSE]
  if (!nrow(conf)) stop("No confounder was found for ", exposure, ": nothing to match on. Check LassoOverviewPlot()$roles.")
  matching <- conf$predictor
  dropped_dup <- roles$predictor[startsWith(roles$role, "near-duplicate")]   # alone or together with other codes
  conf_note <- paste(sprintf("%s (%s)", conf$predictor, sub("confounder \\((.*)\\)", "\\1", conf$role)), collapse = ", ")
  rows <- lapply(seq_along(dependents), function(i) {
    dep <- dependents[i]
    kept <- tab$independent[tab$dependent == dep & tab$selected]
    adjust <- if (outcome_adjust) setdiff(kept, c(matching, dropped_dup, exclude, exposure)) else character(0)
    left_out <- intersect(kept, c(dropped_dup, exclude))
    data.frame(ModelID = paste0("LASSO_", i), DependentVariable = dep,
               IndependentVariable = paste(c(exposure, adjust), collapse = " + "), MatchingVariable = paste(matching, collapse = " + "),
               Method = method, ExposureVariable = exposure, Ratio = NA_real_, PSCaliperSD = NA_real_, CovariateCalipers = NA_character_,
               Notes = paste0("LASSO double selection. Propensity score: confounders ", conf_note, ". ",
                              if (outcome_adjust) sprintf("Outcome model adds %d other covariate%s the %s LASSO kept. ", length(adjust), if (length(adjust) == 1) "" else "s", dep) else "",
                              if (length(left_out)) sprintf("Left out although the %s LASSO kept them: %s (near-duplicate of the exposure or excluded).", dep, paste(left_out, collapse = ", ")) else ""),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$Type <- PSResolveType(DT, out$DependentVariable, type)
  if (method == "lm" && any(out$Type != "Linear")) stop("method lm needs type Linear for every dependent.")
  out[, c("ModelID", "DependentVariable", "IndependentVariable", "MatchingVariable", "Method", "Type", "ExposureVariable",
          "Ratio", "PSCaliperSD", "CovariateCalipers", if (isTRUE(notes)) "Notes")]
}
