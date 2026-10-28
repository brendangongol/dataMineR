#' Split a covariate specification into its terms
#'
#' @param x Covariate expression.
#' @return A character vector of terms.
#' @keywords internal
#' @noRd
PS_SplitVars <- function(x) {
  attr(terms(PS_Formula(x)), "term.labels")
}

#' Outcome model type to family
#'
#' A named vector: Logistic to binomial, Linear to gaussian, Poisson to poisson.
#' @keywords internal
#' @noRd
PS_TypeFamily <- c(Logistic = "binomial", Linear = "gaussian", Poisson = "poisson")

#' Headline effect measure for a family
#'
#' @param family \code{"binomial"}, \code{"gaussian"} or \code{"poisson"}.
#' @return \code{"Odds ratio"}, \code{"Mean difference"} or
#'   \code{"Ratio of expected counts"}.
#' @keywords internal
#' @noRd
PS_MeasureFor <- function(family) c(binomial = "Odds ratio", gaussian = "Mean difference", poisson = "Ratio of expected counts")[[family]]

#' Check and complete a setup table
#'
#' Checks columns and values, and adds ModelID, Method, Type, the outcome-model
#' terms, the Formula column and flag columns.
#'
#' @param DT Data frame.
#' @param setup Setup table.
#' @param weights,cluster,strata Design columns.
#' @param method,type Defaults for rows without Method or Type.
#' @param balance_covariates Covariates whose balance is assessed.
#' @param flags Optional covariate flags.
#' @param verb \code{"matched"}, \code{"weighted"} or \code{"balanced"}, used in
#'   the Formula column.
#' @return A list with the prepared \code{setup} and \code{balance_covariates}.
#' @keywords internal
#' @noRd
PS_PrepareSetup <- function(DT, setup, weights, cluster, strata, method, type, balance_covariates, flags, verb) {
  setup <- as.data.frame(setup, stringsAsFactors = FALSE)
  if (!nrow(setup)) stop("setup must contain at least one model row.")
  need <- c("DependentVariable", "IndependentVariable", "MatchingVariable")
  if (!all(need %in% names(setup))) stop("setup needs the columns: ", paste(need, collapse = ", "))
  if (!"ModelID" %in% names(setup)) setup$ModelID <- seq_len(nrow(setup))
  if (anyNA(setup$ModelID) || any(!nzchar(trimws(as.character(setup$ModelID)))) || anyDuplicated(setup$ModelID))
    stop("ModelID must be unique and nonblank; correct the setup table before fitting.")
  for (nm in need) {
    setup[[nm]] <- trimws(as.character(setup[[nm]]))
    if (anyNA(setup[[nm]]) || any(!nzchar(setup[[nm]]))) stop(nm, " contains blank values.")
  }
  if (!"ExposureVariable" %in% names(setup)) setup$ExposureVariable <- NA_character_
  setup$.IndependentExpression <- setup$IndependentVariable
  setup$.OutcomeAdjust <- ""
  for (i in seq_len(nrow(setup))) {
    vars <- PS_BaseVars(setup$IndependentVariable[i])
    if (PS_Optional(setup$ExposureVariable[i])) {
      if (length(vars) != 1L) stop("Specify ExposureVariable for a multi-predictor IndependentVariable (ModelID ",setup$ModelID[i],").")
      setup$ExposureVariable[i] <- vars
    }
    exposure <- setup$ExposureVariable[i]
    terms <- PS_SplitVars(setup$IndependentVariable[i])
    if (!exposure %in% terms) stop("IndependentVariable must include ExposureVariable as a main effect.")
    if (!exposure %in% names(DT) || !all(DT[[exposure]] %in% c(0,1))) stop("ExposureVariable must be an existing 0/1 column.")
    setup$.OutcomeAdjust[i] <- paste(setdiff(terms, exposure),collapse=" + ")
  }
  setup$Formula <- sprintf("%s ~ %s | %s on: %s", setup$DependentVariable, setup$IndependentVariable, verb,
                           vapply(setup$MatchingVariable, function(x) paste(PS_SplitVars(x), collapse = " + "), character(1)))
  if (is.null(balance_covariates)) balance_covariates <- unique(unlist(lapply(setup$MatchingVariable, PS_BaseVars)))
  setup <- PSRegressionFlags(setup, flags)
  all_vars <- unique(c(setup$DependentVariable, unlist(lapply(setup$IndependentVariable, PS_BaseVars)),
    unlist(lapply(setup$MatchingVariable, PS_BaseVars)), weights, cluster, strata, balance_covariates))
  setup$IndependentVariable <- setup$ExposureVariable
  missing_cols <- setdiff(all_vars, names(DT))
  if (length(missing_cols)) stop("Columns not found in DT: ", paste(missing_cols, collapse = ", "))
  if (length(method) != 1L || is.na(method) || !method %in% c("glm", "svyglm", "lm"))
    stop("method must be glm, svyglm or lm.")
  defaults <- PSResolveType(DT, setup$DependentVariable, type)
  if (!"Method" %in% names(setup)) setup$Method <- method
  if (!"Type" %in% names(setup)) setup$Type <- defaults
  for (nm in c("Method", "Type")) {
    setup[[nm]] <- trimws(as.character(setup[[nm]]))
    if (anyNA(setup[[nm]]) || any(!nzchar(setup[[nm]])))
      stop(nm, " contains blank values; fill the setup column or remove it to use defaults.")
  }
  bad <- !setup$Method %in% c("glm", "svyglm", "lm")
  if (any(bad)) stop("Method must be glm, svyglm or lm (ModelID ", paste(setup$ModelID[bad], collapse = ", "), ").")
  bad <- !setup$Type %in% c("Logistic", "Linear", "Poisson")
  if (any(bad)) stop("Type must be Logistic, Linear or Poisson (ModelID ", paste(setup$ModelID[bad], collapse = ", "), ").")
  bad <- setup$Method == "lm" & setup$Type != "Linear"
  if (any(bad)) stop("Method lm needs Type Linear (ModelID ", paste(setup$ModelID[bad], collapse = ", "), ").")
  if (any(setup$Method == "svyglm") && (is.null(weights) || (is.null(cluster) && is.null(strata))))
    stop("Method svyglm needs weights and cluster and/or strata; use glm or lm without a survey design.")
  list(setup = setup, balance_covariates = balance_covariates)
}

#' Unadjusted (crude) effects for every setup row
#'
#' @param DT Data frame.
#' @param setup Prepared setup table.
#' @param weights,cluster,strata Design columns.
#' @param estimand Estimand used for the crude comparison.
#' @return A list: \code{table} (results by key) and \code{key} (one per row).
#' @keywords internal
#' @noRd
PS_CrudeEffects <- function(DT, setup, weights, cluster, strata, estimand) {
  crude <- list()
  crude_key <- paste(setup$IndependentVariable, setup$DependentVariable, setup$Method == "svyglm", setup$Type)
  for (ck in unique(crude_key)) {
    i <- which(crude_key == ck)[1]; ind <- setup$IndependentVariable[i]; svy <- setup$Method[i] == "svyglm"
    unmatched <- structure(list(approach = "none", estimand = estimand, target_population = if (estimand == "ATT") "Exposed population" else estimand,
                                data = transform(DT, .ps_weight = if (is.null(weights)) 1 else DT[[weights]]),
                                diagnostics = data.frame(review_status = "Unmatched comparison"), exposure = ind,
                                weights = weights, cluster = if (svy) cluster else NULL, strata = if (svy) strata else NULL), class = "PSMatchWeight")
    crude[[ck]] <- tryCatch(PSOutcomeCompile(list(crude = unmatched), setup$DependentVariable[i], family = PS_TypeFamily[[setup$Type[i]]],
                                             PrintSummary = FALSE)$table, error = function(e) NULL)
  }
  list(table = crude, key = crude_key)
}

#' Result rows for setup rows that could not be fitted
#'
#' @param setup Prepared setup table.
#' @param idx Row positions.
#' @param error Error message.
#' @return A data frame of failed rows.
#' @keywords internal
#' @noRd
PS_FailedRows <- function(setup, idx, error) {
  data.frame(ModelID = setup$ModelID[idx], dependent = setup$DependentVariable[idx], independent = setup$IndependentVariable[idx],
             MatchingVariable = setup$MatchingVariable[idx], Method = setup$Method[idx], Type = setup$Type[idx],
             formula = setup$Formula[idx], status = "failed", error = error, stringsAsFactors = FALSE)
}

#' Propensity model, matching and weighting summary columns
#'
#' @param DT Data frame.
#' @param fit A \code{PSFit()} result, or \code{NULL} (entropy balancing).
#' @param m A \code{PSMatchWeight()} or \code{EBWeight()} result.
#' @param ind Exposure column.
#' @param ps_method Propensity model method.
#' @param matching Optional matching-only values.
#' @param weighting Optional weighting-only values.
#' @param quality_label Prefix of the quality note column.
#' @return A one-row data frame.
#' @keywords internal
#' @noRd
PS_PropensitySummary <- function(DT, fit, m, ind, ps_method, matching = NULL, weighting = NULL, quality_label = "matching") {
  smd <- m$balance
  matched_smd <- smd[smd$in_propensity_model, , drop = FALSE]
  finite_max <- function(x) if (any(is.finite(x))) max(x[is.finite(x)]) else NA_real_
  finite_min <- function(x) if (any(is.finite(x))) min(x[is.finite(x)]) else NA_real_
  ps_coef <- NULL
  if (!is.null(fit)) {
    old <- options(survey.lonely.psu = "adjust")
    cf <- tryCatch(summary(fit$model)$coefficients, error = function(e) NULL)
    options(old)
    ps_coef <- if (is.null(cf)) NA_character_ else paste(sprintf("%s: b = %.3g, OR = %.3g, p = %.2g", rownames(cf)[-1],
      cf[-1, 1], exp(cf[-1, 1]), cf[-1, ncol(cf)]), collapse = "; ")
  }
  dg <- m$diagnostics; ex <- dg$group == "exposed"; un <- dg$group == "unexposed"
  info <- list(
    ps_model = if (!is.null(fit)) sprintf("%s: %s", ps_method, deparse1(fit$formula)), ps_model_coefficients = ps_coef,
    ps_AUC_in_sample_unweighted = if (!is.null(fit)) tryCatch(as.numeric(pROC::roc(DT[[ind]], fit$ps, levels = c(0, 1), direction = "<", quiet = TRUE)$auc), error = function(e) NA_real_),
    weights_model = weighting$weights_model,
    match_ratio = matching$match_ratio, caliper_sd = matching$caliper_sd, caliper_logit = matching$caliper_logit,
    n_exposed_before = dg$n_before[ex], n_controls_before = dg$n_before[un], n_matched_sets = matching$n_matched_sets,
    n_controls_unused = matching$n_controls_unused, controls_per_exposed = matching$controls_per_exposed,
    exposed_retained_weighted = matching$exposed_retained_weighted,
    ess_exposed_before = dg$ess_before[ex], ess_exposed_after = dg$ess_after[ex],
    ess_controls_before = dg$ess_before[un], ess_controls_after = dg$ess_after[un],
    estimand = weighting$estimand, trim = weighting$trim, moments = weighting$moments,
    ess_kept_exposed = weighting$ess_kept_exposed, ess_kept_controls = weighting$ess_kept_controls,
    ps_min_exposed = dg$ps_min[ex], ps_max_exposed = dg$ps_max[ex], ps_min_controls = dg$ps_min[un], ps_max_controls = dg$ps_max[un],
    outside_support_exposed = dg$outside_common_support[ex], outside_support_controls = dg$outside_common_support[un],
    ps_extreme_exposed = dg$fraction_ps_extreme[ex], ps_extreme_controls = dg$fraction_ps_extreme[un],
    max_var_ratio_after = finite_max(smd$var_ratio_after), min_var_ratio_after = finite_min(smd$var_ratio_after),
    worst_var_ratio_fold_after = if (all(is.na(smd$var_ratio_after))) NA_real_ else max(pmax(smd$var_ratio_after, 1/smd$var_ratio_after), na.rm = TRUE),
    max_smd_matched_before = finite_max(abs(matched_smd$smd_before)),
    max_smd_matched_after = finite_max(abs(matched_smd$smd_after)),
    balanced_matched_variables = nrow(matched_smd) > 0 && all(matched_smd$balanced),
    balanced_all_assessed = nrow(smd) > 0 && all(smd$balanced),
    n_assessed_covariates = nrow(smd), n_imbalanced_matched = sum(!matched_smd$balanced),
    weight_min_exposed = dg$weight_min[ex], weight_max_exposed = dg$weight_max[ex],
    weight_min_controls = dg$weight_min[un], weight_max_controls = dg$weight_max[un],
    controls_retained = matching$controls_retained, controls_retained_weighted = matching$controls_retained_weighted,
    AUC_note = "Both AUCs are apparent, unweighted discrimination; neither measures matching quality or validated prediction",
    quality_note = "Heuristic thresholds only; inspect all assessed covariates, distributions and overlap",
    max_ks_after = if (all(is.na(smd$ks_after))) NA_real_ else max(smd$ks_after, na.rm = TRUE),
    ps_fit_warnings = dg$fit_warnings[1])
  info <- as.data.frame(info[!vapply(info, is.null, logical(1))], stringsAsFactors = FALSE)   # drop the other approach's values
  names(info)[names(info) == "quality_note"] <- paste0(quality_label, "_quality_note")
  for (j in seq_len(nrow(smd))) {
    for (metric in c("smd_before", "smd_after", "mean_control_before", "mean_exposed_before",
                     "mean_control_after", "mean_exposed_after", "var_ratio_before", "var_ratio_after",
                     "ks_before", "ks_after", "balanced", "in_propensity_model"))
      info[[paste0(metric, "_", smd$variable[j])]] <- smd[[metric]][j]
  }
  info
}

#' One result row for a setup row
#'
#' @param DT Data frame.
#' @param setup Prepared setup table.
#' @param i Row position.
#' @param eff \code{PSOutcomeCompile()} results by Type.
#' @param m Matching or weighting result.
#' @param fit Propensity fit, or \code{NULL}.
#' @param crude \code{PS_CrudeEffects()} result.
#' @param covs Matching covariates.
#' @param balance_covariates Covariates whose balance is assessed.
#' @param info \code{PS_PropensitySummary()} row.
#' @param retention Matching-only counts, or \code{NULL}.
#' @param quality One-column quality verdict.
#' @return A list: \code{row} and \code{outcome_balance}.
#' @keywords internal
#' @noRd
PS_EffectRow <- function(DT, setup, i, eff, m, fit, crude, covs, balance_covariates, info, retention, quality) {
  dep <- setup$DependentVariable[i]; tp <- setup$Type[i]; ind <- setup$IndependentVariable[i]; meth <- setup$Method[i]
  eff <- eff[[tp]]
  r <- eff$table[eff$table$dependent == dep, ]
  smd <- m$balance
  base <- data.frame(ModelID = setup$ModelID[i], dependent = dep, independent = ind, MatchingVariable = setup$MatchingVariable[i],
                     Method = meth, Type = tp, formula = setup$Formula[i], stringsAsFactors = FALSE)
  ob <- PSOutcomeBalance(DT, m, dep, setdiff(unique(c(balance_covariates, PS_BaseVars(covs))), ind),
                         outcome_predictors = PS_BaseVars(setup$.OutcomeAdjust[i]))
  outcome_balance <- data.frame(ModelID = setup$ModelID[i], dependent = dep, independent = ind, MatchingVariable = setup$MatchingVariable[i], Method = meth, Type = tp, ob$balance)
  full <- smd[match(outcome_balance$variable, smd$variable), , drop = FALSE]
  outcome_balance$full_smd_after <- full$smd_after
  outcome_balance$full_mean_control_after <- full$mean_control_after
  outcome_balance$full_mean_exposed_after <- full$mean_exposed_after
  outcome_balance$full_balanced <- full$balanced
  outcome_balance$in_propensity_model <- full$in_propensity_model
  if (!nrow(r) || is.na(r$family[1])) {
    row <- do.call(cbind, Filter(Negate(is.null), list(base,
      data.frame(status = "failed", error = if (nrow(r)) r$error[1] else "Outcome result missing", stringsAsFactors = FALSE), quality,
      data.frame(target_population = m$target_population, stringsAsFactors = FALSE), retention,
      data.frame(n_missing_outcome = sum(!complete.cases(m$data[, unique(c(dep, PS_BaseVars(setup$.OutcomeAdjust[i]))), drop=FALSE]))), info, ob$summary)))
    return(list(row = row, outcome_balance = outcome_balance))
  }
  fam <- r$family[1]; want <- PS_MeasureFor(fam)
  all_r <- r
  r <- r[r$measure == want, ][1, ]
  ct <- crude$table[[crude$key[i]]]
  cr <- if (is.null(ct)) data.frame(estimate = NA, lower = NA, upper = NA, p_value = NA) else ct[ct$measure == want, ][1, ]
  other <- function(meas, col) { v <- all_r[all_r$measure == meas, col]; if (length(v)) v[1] else NA_real_ }
  md <- m$data[complete.cases(m$data[, unique(c(dep, PS_BaseVars(setup$.OutcomeAdjust[i]))), drop=FALSE]), ]
  binary_dep <- all(md[[dep]] %in% c(0, 1))
  model <- eff$models[[1]][[dep]]
  outcome_info <- data.frame(
    outcome_model = sprintf("%s (%s): %s", if (meth == "lm") "glm (Gaussian identity; lm option)" else class(model)[1], fam, deparse1(formula(model))),
    outcome_exposed = r$mean_exposed, outcome_controls = r$mean_unexposed,
    events_exposed = if (binary_dep) sum(md[[dep]][md[[ind]] == 1]) else NA_real_,
    events_controls = if (binary_dep) sum(md[[dep]][md[[ind]] == 0]) else NA_real_,
    RiskDifference = other("Risk difference", "estimate"), `RiskDifference 2.5 %` = other("Risk difference", "lower"),
    `RiskDifference 97.5 %` = other("Risk difference", "upper"),
    RiskRatio = other("Risk ratio", "estimate"), `RiskRatio 2.5 %` = other("Risk ratio", "lower"), `RiskRatio 97.5 %` = other("Risk ratio", "upper"),
    MeanDifference = if (fam == "binomial") NA_real_ else other("Mean difference", "estimate"),
    `MeanDifference 2.5 %` = if (fam == "binomial") NA_real_ else other("Mean difference", "lower"),
    `MeanDifference 97.5 %` = if (fam == "binomial") NA_real_ else other("Mean difference", "upper"),
    `Survey-only 2.5 %` = r$survey_lower, `Survey-only 97.5 %` = r$survey_upper,
    inference = r$inference, inference_df = r$inference_df, outcome_fit_warnings = r$fit_warnings,
    target_population = r$target_population, inference_note = r$inference_note,
    n_missing_outcome = r$n_missing_outcome, outcome_review_status = r$review_status,
    outcome_measure = want,
    EffectEstimate = r$estimate, EffectLower95 = r$lower, EffectUpper95 = r$upper,
    effect_scale = if (fam == "gaussian") "Difference" else "Ratio",
    coefficient_interval_scale = if (fam == "gaussian") "Difference" else "Log ratio",
    CrudeLower = cr$lower, CrudeUpper = cr$upper, CrudeP = cr$p_value, check.names = FALSE, stringsAsFactors = FALSE)
  tr <- if (fam != "gaussian") log else identity
  zq <- qt(0.975, df = r$inference_df)
  est <- tr(r$estimate); se <- (tr(r$upper) - tr(r$lower)) / (2 * zq)
  auc <- if (fam == "binomial") {
    tryCatch(as.numeric(suppressMessages(pROC::roc(model$y, fitted(model), levels = c(0, 1), direction = "<", quiet = TRUE))$auc), error = function(e) NA_real_)
  } else NA_real_
  row <- do.call(cbind, Filter(Negate(is.null), list(base,   # retention is NULL for weighting
               data.frame(Estimate = est, `Std. Error` = se, `t value` = est / se, `Pr(>|t|)` = r$p_value, outcome_AUC_in_sample_unweighted = auc,
                          OddsRatio = if (fam == "binomial") r$estimate else NA_real_,
                          logOddsRatio = if (fam == "binomial") est else NA_real_,
                          OddsRatioP = if (fam == "binomial") r$p_value else NA_real_,
                          ExpectedCountRatio = if (fam == "poisson") r$estimate else NA_real_,
                          names = ind, NcountsNoNA = r$records, `2.5 %` = tr(r$lower), `97.5 %` = tr(r$upper),
                          CrudeRatio = cr$estimate, n_exposed = r$n_exposed, n_controls = r$n_control, check.names = FALSE, stringsAsFactors = FALSE),
               retention,
               as.data.frame(Filter(Negate(is.null), list(ps_at_bounds = fit$n_clipped,   # no ps_at_bounds without a propensity score
                          max_smd_before = max(abs(smd$smd_before), na.rm = TRUE), max_smd_after = max(abs(smd$smd_after), na.rm = TRUE),
                          balanced = all(smd$balanced), n_imbalanced = sum(!smd$balanced)))),
               quality,
               data.frame(status = r$status, error = r$error, stringsAsFactors = FALSE),
               outcome_info, info, ob$summary)))
  list(row = row, outcome_balance = outcome_balance)
}

#' Compile result rows and add multiple-testing corrections
#'
#' Bonferroni and FDR within each dependent and across all rows, in the column
#' order of \code{RegressionCompile()}.
#'
#' @param rows List of result data frames.
#' @param setup Prepared setup table.
#' @param extra_numeric Numeric columns to add when every fit failed.
#' @return A data.table of results.
#' @keywords internal
#' @noRd
Matching_CompileResults <- function(rows, setup, extra_numeric = NULL) {
  columns <- unique(unlist(lapply(rows, names)))
  res <- data.table::rbindlist(lapply(rows, function(r) { for (nm in setdiff(columns, names(r))) r[[nm]] <- NA; r[, columns, drop = FALSE] }), fill = TRUE)
  # Always expose numeric result columns, including when every fit failed.
  for (nm in c("Estimate", "Std. Error", "t value", "Pr(>|t|)", "2.5 %", "97.5 %",
               "OddsRatio", "ExpectedCountRatio", "NcountsNoNA",
               "n_exposed", "n_controls", extra_numeric,
               "max_smd_before", "max_smd_after", "max_smd_matched_before", "max_smd_matched_after",
               "max_var_ratio_after", "min_var_ratio_after", "max_ks_after", "n_imbalanced",
               "ess_exposed_before", "ess_exposed_after", "ess_controls_before", "ess_controls_after",
               "n_missing_outcome", "inference_df", "EffectEstimate", "EffectLower95", "EffectUpper95"))
    if (!nm %in% names(res)) res[[nm]] <- rep(NA_real_, nrow(res))
  res <- res[order(match(res$ModelID, setup$ModelID))]
  for (nm in c("flagged_covariates", "covariate_flags", "has_covariate_flags"))
    res[[nm]] <- setup[[nm]][match(res$ModelID, setup$ModelID)]
  p <- res[["Pr(>|t|)"]]
  res$Significance <- ifelse(is.na(p), NA, ifelse(p < 0.05, "Significant", "Not-Significant"))
  res$p.adjust <- ave(p, res$dependent, FUN = function(x) p.adjust(x, "bonferroni"))
  res$BonfSignificance <- ifelse(res$p.adjust < 0.05, "Significant", "Not-Significant")
  res$FDR <- ave(p, res$dependent, FUN = function(x) p.adjust(x, "BH"))
  res$BHSignificance <- ifelse(res$FDR < 0.05, "Significant", "Not-Significant")
  res$p_bonferroni <- res$p.adjust; res$p_fdr <- res$FDR
  res$sig_unadjusted <- p < 0.05; res$sig_bonferroni <- res$p_bonferroni < 0.05; res$sig_fdr <- res$p_fdr < 0.05
  res$p_fdr_overall <- p.adjust(p, "BH"); res$p_bonferroni_overall <- p.adjust(p, "bonferroni"); res$sig_fdr_overall <- res$p_fdr_overall < 0.05
  regression_cols <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)", "outcome_AUC_in_sample_unweighted", "OddsRatio", "logOddsRatio", "OddsRatioP", "ExpectedCountRatio",
                       "names", "formula", "NcountsNoNA", "2.5 %", "97.5 %", "Significance", "p.adjust", "BonfSignificance", "FDR",
                       "BHSignificance", "p_bonferroni", "p_fdr", "sig_unadjusted", "sig_bonferroni", "sig_fdr", "dependent",
                       "independent", "p_fdr_overall", "p_bonferroni_overall", "sig_fdr_overall")
  data.table::setcolorder(res, c(intersect(regression_cols, names(res)), setdiff(names(res), regression_cols)))
  res
}

#' Keep the compact result columns
#'
#' @param res Compiled results.
#' @param detail Logical; \code{TRUE} keeps every column.
#' @return The results with the core columns, keeping the balance attributes.
#' @keywords internal
#' @noRd
Matching_SelectColumns <- function(res, detail = FALSE) {
  if (isTRUE(detail)) return(res)
  core <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)", "outcome_AUC_in_sample_unweighted", "OddsRatio", "logOddsRatio", "OddsRatioP",
            "ExpectedCountRatio", "names", "formula", "NcountsNoNA", "2.5 %", "97.5 %", "Significance",
            "p.adjust", "BonfSignificance", "FDR", "BHSignificance", "p_bonferroni", "p_fdr", "sig_unadjusted", "sig_bonferroni", "sig_fdr",
            "dependent", "independent", "p_fdr_overall", "p_bonferroni_overall", "sig_fdr_overall",
            # what was fitted
            "ModelID", "ExposureVariable", "IndependentVariable", "MatchingVariable", "Method", "Type",
            "Ratio", "PSCaliperSD", "CovariateCalipers", "estimand", "trim", "moments",
            "status", "error",
            # effect measures on their natural scale
            "outcome_exposed", "outcome_controls", "RiskDifference", "RiskDifference 2.5 %", "RiskDifference 97.5 %",
            "RiskRatio", "RiskRatio 2.5 %", "RiskRatio 97.5 %", "MeanDifference", "MeanDifference 2.5 %", "MeanDifference 97.5 %",
            "CrudeRatio", "CrudeLower", "CrudeUpper",
            # quality checks
            "n_exposed", "n_controls", "n_exposed_unmatched", "exposed_retained", "ess_exposed_after", "ess_controls_after",
            "ps_at_bounds", "outside_support_exposed", "outside_support_controls", "weight_max_exposed", "weight_max_controls",
            "max_smd_after", "n_imbalanced", "matching_quality", "weighting_quality", "balancing_quality")
  out <- data.table::as.data.table(as.data.frame(res)[, intersect(core, names(res)), drop = FALSE])
  for (nm in c("balance", "outcome_balance", "flags")) attr(out, nm) <- attr(res, nm)
  out
}

#' Compact Balance sheet for one outcome
#'
#' One row per covariate per ModelID: balance on the whole matched or weighted
#' sample and on the records with this outcome.
#'
#' @param sub Result rows for the outcome.
#' @param balance Balance on the whole sample.
#' @param outcome_balance Balance on each outcome's records.
#' @param dep Outcome.
#' @return A data frame.
#' @keywords internal
#' @noRd
Matching_BalanceSheet <- function(sub, balance, outcome_balance, dep) {
  empty <- data.frame()
  if (is.null(outcome_balance) || !nrow(outcome_balance) || !"dependent" %in% names(outcome_balance)) return(empty)
  sub <- as.data.frame(sub)
  ob <- as.data.frame(outcome_balance)
  ob <- ob[ob$dependent == dep, , drop = FALSE]
  if (!nrow(ob)) return(empty)
  if (!"ModelID" %in% names(ob)) ob$ModelID <- NA
  model_key <- function(d) paste(d$independent, d$MatchingVariable, d$Method, d$Type, sep = "\r")
  unknown <- is.na(ob$ModelID)
  if (any(unknown)) ob$ModelID[unknown] <- sub$ModelID[match(model_key(ob)[unknown], model_key(sub))]
  ob <- ob[ob$ModelID %in% sub$ModelID & !duplicated(paste(ob$ModelID, ob$variable, sep = "\r")), , drop = FALSE]
  if (!nrow(ob)) return(empty)
  from_balance <- c(full_smd_after = "smd_after", full_mean_control_after = "mean_control_after",
                    full_mean_exposed_after = "mean_exposed_after", full_balanced = "balanced",
                    in_propensity_model = "in_propensity_model")
  b <- if (is.null(balance) || !nrow(balance)) NULL else as.data.frame(balance)
  if (!is.null(b)) at <- match(paste(ob$independent, ob$MatchingVariable, ob$Method, ob$variable, sep = "\r"),
                               paste(b$independent, b$MatchingVariable, b$Method, b$variable, sep = "\r"))
  for (nm in names(from_balance)) {
    if (!nm %in% names(ob)) ob[[nm]] <- NA
    gap <- is.na(ob[[nm]])
    if (!is.null(b) && any(gap)) ob[[nm]][gap] <- b[[from_balance[[nm]]]][at[gap]]
  }
  note <- if ("diagnostic_error" %in% names(ob)) as.character(ob$diagnostic_error) else rep("", nrow(ob))
  note[is.na(note)] <- ""
  out <- data.frame(ModelID = ob$ModelID, variable = ob$variable, in_propensity_model = ob$in_propensity_model,
                    mean_exposed_after = ob$full_mean_exposed_after, mean_control_after = ob$full_mean_control_after,
                    smd_before = ob$smd_before, smd_after = ob$full_smd_after, balanced = ob$full_balanced,
                    smd_after_outcome_sample = ob$smd_after,
                    balanced_outcome_sample = ifelse(nzchar(note), NA, ob$balanced), stringsAsFactors = FALSE)
  if (any(nzchar(note))) out$outcome_sample_note <- note
  out <- out[order(match(out$ModelID, sub$ModelID)), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Definitions of every result column, used by the README
#'
#' @return A data frame: Sheet, Group, Column, Meaning, HowToRead.
#' @keywords internal
#' @noRd
Matching_ColumnGuide <- function() {
  rows <- list()
  add <- function(sheet, group, column, what, how = "")
    rows[[length(rows) + 1L]] <<- data.frame(Sheet = sheet, Group = group, Column = column, Meaning = what, HowToRead = how,
                                              stringsAsFactors = FALSE)
  RA <- "RegressionAnalysis"
  g_id <- "1. Which model this row is"; g_eff <- "2. Effect of the exposure on the outcome"
  g_prec <- "3. How precise the effect is"; g_sig <- "4. Statistical significance"
  g_crude <- "5. Comparison with no adjustment"; g_n <- "6. How many patients were used"
  g_bal <- "7. Did the adjustment work? (balance and verdict)"; g_over <- "8. Overlap, weights and model fit"
  g_stat <- "9. Warnings and status"; g_set <- "10. Settings (detail = TRUE)"
  add(RA, g_id, "ModelID", "Row number from your setup table. One ModelID is one combination of outcome, exposure, matching or weighting variables, method and type.",
      "Use it to find the same model in the Balance sheet and in Matching_RegressionSetup.xlsx.")
  add(RA, g_id, "dependent", "The outcome being analyzed.", "The same in every row of a workbook; each outcome has its own workbook.")
  add(RA, g_id, "independent", "The exposure, coded 0/1 (1 = exposed, 0 = comparison patients).",
      "Every effect compares exposed with comparison patients. A ratio above 1, or a positive difference, means the outcome is higher in the exposed group.")
  add(RA, g_id, "names", "The exposure term of the outcome model; repeats independent.", "Can be ignored.")
  add(RA, g_id, "MatchingVariable", "The covariate(s) used to build the propensity score or the weights; several are joined with +.",
      "Balance is only targeted for these. Covariates that are not listed can still be imbalanced, so check every row of the Balance sheet.")
  add(RA, g_id, "Method", "How the models were fitted: glm (ordinary model, survey weights used as case weights), svyglm (survey design model that uses the cluster and strata) or lm (linear model).",
      "svyglm accounts for survey clustering and strata; glm does not.")
  add(RA, g_id, "Type", "The kind of outcome model: Logistic (yes/no outcome), Linear (continuous outcome) or Poisson (count outcome such as length of stay).",
      "Decides which effect column is filled: OddsRatio for Logistic, MeanDifference for Linear, ExpectedCountRatio for Poisson.")
  add(RA, g_id, "formula", "The model in words: outcome ~ exposure | matched, weighted or balanced on: the variables used.",
      "Read it as: the outcome is compared between exposed and comparison patients after adjusting for the listed variables.")
  add(RA, g_id, "estimand", "Which patients the effect describes: ATT (patients who were exposed), ATE (all patients) or ATO (patients with the most overlap between the groups).",
      "ATT answers: what was the effect among those actually exposed? ATE answers: what would change if everyone were exposed rather than no one? Only compare results that share an estimand.")
  add(RA, g_id, "trim", "Weighting only: weights above this quantile were set equal to it. Blank means no cap.",
      "Capping limits the influence of a few very large weights but can leave some imbalance.")
  add(RA, g_id, "moments", "Entropy balancing only: 1 means the group averages were made equal; 2 means averages and variances were made equal.",
      "2 balances more, but is harder to achieve and can leave fewer effective patients.")
  add(RA, g_eff, "OddsRatio", "Logistic outcomes only: odds of the outcome if the target patients were exposed, divided by the odds if they were not (equals exp(Estimate)).",
      "1 = no difference. 1.5 = 50% higher odds in the exposed group; 0.8 = 20% lower odds. Odds are not risks: when the outcome is common the odds ratio looks further from 1 than the risk ratio does.")
  add(RA, g_eff, "logOddsRatio", "Logistic outcomes only: natural logarithm of OddsRatio (same as Estimate).", "0 = no difference.")
  add(RA, g_eff, "OddsRatioP", "Logistic outcomes only: p-value for the odds ratio (same as Pr(>|t|)).", "See Pr(>|t|).")
  add(RA, g_eff, "ExpectedCountRatio", "Poisson outcomes (counts, such as days or events): expected count if exposed divided by expected count if not exposed.",
      "1 = no difference. 1.2 = 20% higher expected count in the exposed group; 0.9 = 10% lower.")
  add(RA, g_eff, "Estimate", "The effect on the scale the model uses: the natural logarithm of the ratio for Logistic and Poisson outcomes, or the plain difference in the outcome's units for Linear outcomes.",
      "0 = no difference between the groups. Positive = higher outcome in the exposed group; negative = lower. For ratio outcomes use exp(Estimate) to get the ratio (shown in OddsRatio or ExpectedCountRatio).")
  add(RA, g_eff, "outcome_exposed", "Average outcome expected for the target patients if all of them were exposed, estimated from the outcome model. A risk (proportion) for yes/no outcomes.",
      "One of the two numbers the effect measures compare. It should look plausible next to the outcome rate you observe in exposed patients.")
  add(RA, g_eff, "outcome_controls", "Average outcome expected for the target patients if none of them were exposed, estimated from the outcome model. A risk (proportion) for yes/no outcomes.",
      "The other number the effect measures compare.")
  add(RA, g_eff, "RiskDifference", "Logistic outcomes only: outcome_exposed minus outcome_controls, in proportion units.",
      "0 = no difference. 0.02 = 2 percentage points more people with the outcome in the exposed group.")
  add(RA, g_eff, "RiskDifference 2.5 %", "Lower 95% confidence limit of RiskDifference.", "If the interval includes 0 the data are compatible with no difference.")
  add(RA, g_eff, "RiskDifference 97.5 %", "Upper 95% confidence limit of RiskDifference.", "")
  add(RA, g_eff, "RiskRatio", "Logistic outcomes only: outcome_exposed divided by outcome_controls.",
      "1 = no difference. 1.3 = a 30% higher risk in the exposed group. Easier to interpret than OddsRatio when the outcome is common.")
  add(RA, g_eff, "RiskRatio 2.5 %", "Lower 95% confidence limit of RiskRatio.", "If the interval includes 1 the data are compatible with no difference.")
  add(RA, g_eff, "RiskRatio 97.5 %", "Upper 95% confidence limit of RiskRatio.", "")
  add(RA, g_eff, "MeanDifference", "Linear and Poisson outcomes: outcome_exposed minus outcome_controls, in the outcome's own units.",
      "0 = no difference. 2.5 = the exposed group averages 2.5 units more than the comparison group.")
  add(RA, g_eff, "MeanDifference 2.5 %", "Lower 95% confidence limit of MeanDifference.", "If the interval includes 0 the data are compatible with no difference.")
  add(RA, g_eff, "MeanDifference 97.5 %", "Upper 95% confidence limit of MeanDifference.", "")
  add(RA, g_prec, "Std. Error", "Standard error of Estimate. It treats the matching or weights as fixed.", "Smaller = more precise. Estimate +/- 1.96 x Std. Error gives the 95% confidence limits.")
  add(RA, g_prec, "t value", "Estimate divided by Std. Error.", "Beyond about +/-2 corresponds to p below 0.05 in large samples. Use Pr(>|t|) for the p-value itself.")
  add(RA, g_prec, "2.5 %", "Lower 95% confidence limit of Estimate, on the same scale as Estimate.",
      "The range of effects compatible with the data. For ratio outcomes this is a logarithm: use exp() of the two limits (in Excel, =EXP(cell)) to get limits for the ratio. If the interval includes 0 the data are compatible with no difference.")
  add(RA, g_prec, "97.5 %", "Upper 95% confidence limit of Estimate, on the same scale as Estimate.", "See 2.5 %.")
  add(RA, g_sig, "Pr(>|t|)", "Unadjusted p-value for the hypothesis that the exposure has no effect.",
      "Small values (conventionally below 0.05) mean the data would be unusual if there were no effect. It says nothing about how large or important the effect is. If you ran many models use the adjusted columns below.")
  add(RA, g_sig, "Significance", "Text label: Significant when Pr(>|t|) is below 0.05, otherwise Not-Significant.", "A convenience label only; see the cautions under Pr(>|t|).")
  add(RA, g_sig, "p.adjust", "Bonferroni-adjusted p-value: Pr(>|t|) multiplied by the number of models for this outcome (at most 1).",
      "Very strict. Use it when even one false positive would be costly. Compare with 0.05.")
  add(RA, g_sig, "BonfSignificance", "Text label: Significant when p.adjust is below 0.05.", "")
  add(RA, g_sig, "FDR", "Benjamini-Hochberg false discovery rate adjusted p-value, across the models for this outcome.",
      "Below 0.05 means roughly 5% or fewer of the results called significant are expected to be false findings. Less strict than Bonferroni; the usual choice when screening many models.")
  add(RA, g_sig, "BHSignificance", "Text label: Significant when FDR is below 0.05.", "")
  add(RA, g_sig, "p_bonferroni", "Same as p.adjust (name used by code that reads the results).", "")
  add(RA, g_sig, "p_fdr", "Same as FDR (name used by code that reads the results).", "")
  add(RA, g_sig, "sig_unadjusted", "TRUE when Pr(>|t|) is below 0.05.", "")
  add(RA, g_sig, "sig_bonferroni", "TRUE when p_bonferroni is below 0.05.", "")
  add(RA, g_sig, "sig_fdr", "TRUE when p_fdr is below 0.05.", "")
  add(RA, g_sig, "p_fdr_overall", "Benjamini-Hochberg adjusted p-value across every model in the whole analysis (all outcomes), not only this outcome.",
      "Use this when you want one correction for the entire set of results.")
  add(RA, g_sig, "p_bonferroni_overall", "Bonferroni-adjusted p-value across every model in the whole analysis.", "The strictest correction.")
  add(RA, g_sig, "sig_fdr_overall", "TRUE when p_fdr_overall is below 0.05.", "")
  add(RA, g_crude, "CrudeRatio", "The same effect measure with no matching or weighting (survey weights only): an odds ratio for Logistic, a ratio of expected counts for Poisson, a difference for Linear.",
      "Compare with the adjusted estimate. A large change shows the groups differed in ways that matter; a similar value suggests the adjustment changed little.")
  add(RA, g_crude, "CrudeLower", "Lower 95% confidence limit of CrudeRatio.", "")
  add(RA, g_crude, "CrudeUpper", "Upper 95% confidence limit of CrudeRatio.", "")
  add(RA, g_n, "NcountsNoNA", "Patients used to estimate the effect: analysis sample patients with the outcome recorded.", "Compare with the size of your data to see how many patients were left out.")
  add(RA, g_n, "n_exposed", "Exposed patients used to estimate the effect.", "Small numbers give imprecise effects.")
  add(RA, g_n, "n_controls", "Comparison patients used to estimate the effect.", "For matching this is the matched controls only, so it can be smaller than the number of controls you started with.")
  add(RA, g_n, "n_exposed_unmatched", "Matching only: exposed patients for whom no comparison patient was found within the caliper. They are left out.",
      "Should be small. If large, the result describes only the exposed patients who could be matched.")
  add(RA, g_n, "exposed_retained", "Matching only: share of exposed patients that were matched (0 to 1).", "Close to 1 is best. Below 0.8 the matching is rated poor.")
  add(RA, g_n, "ess_exposed_after", "Effective sample size of the exposed group after matching or weighting: the number of equally weighted patients that would carry the same information.",
      "Compare with n_exposed. If it is much smaller (for example under a quarter) a few patients carry large weights and the estimate is unstable.")
  add(RA, g_n, "ess_controls_after", "Effective sample size of the comparison group after matching or weighting.", "See ess_exposed_after.")
  add(RA, g_bal, "max_smd_after", "Largest absolute standardized mean difference (SMD) across all assessed covariates after matching or weighting.",
      "Below 0.1 for every covariate is the usual target. The Balance sheet shows which covariate is responsible.")
  add(RA, g_bal, "n_imbalanced", "Number of assessed covariates whose absolute SMD is still 0.1 or more after matching or weighting.", "0 is the goal.")
  add(RA, g_bal, "outcome_max_smd_after", "max_smd_after recomputed on only the patients who have this outcome recorded.",
      "Patients with a missing outcome are dropped from the effect estimate, which can upset the balance. If this is larger than max_smd_after, balance is worse for this outcome than for the whole analysis sample.")
  add(RA, g_bal, "matching_quality", "Automatic verdict on the matching, combining the share of exposed patients matched, propensity scores at the bounds and covariate balance.",
      "Starts with 'poor' when fewer than 80% of exposed patients were matched or any score is exactly 0 or 1: do not rely on the effect until you know why. 'fair' means matching worked but at least one covariate still has SMD 0.1 or more: check the Balance sheet. 'SMD/retention thresholds passed' means all automatic checks passed. These are rules of thumb, not proof that the comparison is valid.")
  add(RA, g_bal, "weighting_quality", "Automatic verdict on the weighting, combining the effective sample size kept, propensity scores at the bounds and covariate balance.",
      "Starts with 'poor' when either group keeps under 25% of its effective sample size or any score is exactly 0 or 1: a few patients dominate. 'fair' means at least one covariate still has SMD 0.1 or more: check the Balance sheet. 'SMD/ESS thresholds passed' means all automatic checks passed. These are rules of thumb, not proof that the comparison is valid.")
  add(RA, g_bal, "balancing_quality", "Automatic verdict on the entropy balancing, combining the effective sample size kept and the balance reached.",
      "Starts with 'poor' when either group keeps under 25% of its effective sample size, or a balanced covariate still has SMD 0.01 or more (usually failed convergence; ps_fit_warnings, shown with detail = TRUE, gives the reason). 'fair' means at least one other covariate has SMD 0.1 or more. 'SMD/ESS thresholds passed' means all automatic checks passed.")
  add(RA, g_over, "ps_at_bounds", "Number of patients whose propensity score is exactly 0 or 1.",
      "Should be 0. Anything above 0 means a matching variable almost perfectly predicts the exposure, so some patients have no comparable counterpart.")
  add(RA, g_over, "outside_support_exposed", "Exposed patients whose propensity score lies outside the range of scores shared with the comparison group.", "Should be small. Large counts mean the groups barely overlap and the effect for those patients is extrapolated.")
  add(RA, g_over, "outside_support_controls", "Comparison patients whose propensity score lies outside the range of scores shared with the exposed group.", "See outside_support_exposed.")
  add(RA, g_over, "weight_max_exposed", "Largest analysis weight among exposed patients (balancing weight x survey weight).", "Compare with a typical weight. A very large maximum means a few patients dominate the result.")
  add(RA, g_over, "weight_max_controls", "Largest analysis weight among comparison patients.", "See weight_max_exposed.")
  add(RA, g_over, "outcome_AUC_in_sample_unweighted", "Logistic outcomes only: area under the ROC curve of the outcome model, in the analysis sample, unweighted.",
      "0.5 = no better than chance; 1 = perfect separation. It describes how well the model separates patients with and without the outcome. It does not show that the effect estimate is valid.")
  add(RA, g_stat, "ps_fit_warnings", "Warnings raised while fitting the propensity score, matching or balancing.", "Blank = none. Convergence warnings mean the scores or weights may be unreliable.")
  add(RA, g_stat, "flagged_covariates", "Matching variables that the pre-fit screen (Flags sheet) marked for review for this exposure.",
      "Blank = none flagged. A flagged variable may nearly define the exposure or be measured after it; see the Flags sheet.")
  add(RA, g_stat, "status", "ok when the model was fitted, failed otherwise.", "If failed, no effect was estimated; read error.")
  add(RA, g_stat, "error", "Reason a model failed.", "Blank when status is ok.")
  add(RA, g_eff, "EffectEstimate", "The headline effect on its natural scale: a ratio (odds ratio, ratio of expected counts) or, for Linear outcomes, a difference.", "Same number as OddsRatio, ExpectedCountRatio or MeanDifference.")
  add(RA, g_eff, "EffectLower95", "Lower 95% confidence limit of EffectEstimate.", "If the interval includes 1 (ratio) or 0 (difference) the data are compatible with no difference.")
  add(RA, g_eff, "EffectUpper95", "Upper 95% confidence limit of EffectEstimate.", "")
  add(RA, g_eff, "outcome_measure", "Name of the headline effect measure.", "")
  add(RA, g_eff, "effect_scale", "Ratio or Difference.", "Tells you whether the headline effect is read against 1 (Ratio) or 0 (Difference).")
  add(RA, g_eff, "coefficient_interval_scale", "Scale of Estimate and its limits: Log ratio or Difference.", "")
  add(RA, g_eff, "events_exposed", "Yes/no outcomes: exposed patients with the outcome in the analysis sample (unweighted count).", "")
  add(RA, g_eff, "events_controls", "Yes/no outcomes: comparison patients with the outcome in the analysis sample (unweighted count).", "Very few events make estimates unstable.")
  add(RA, g_prec, "Survey-only 2.5 %", "Lower 95% limit of the headline effect using the model variance only, without the matched set clustering.", "Compare with the main interval; a much narrower interval shows the matched set clustering matters.")
  add(RA, g_prec, "Survey-only 97.5 %", "Upper 95% limit of the headline effect using the model variance only.", "")
  add(RA, g_prec, "inference", "How the standard errors were computed.", "")
  add(RA, g_prec, "inference_df", "Degrees of freedom used for the confidence intervals and p-values.", "Small values give wider intervals.")
  add(RA, g_prec, "inference_note", "What the confidence intervals do not include.", "The intervals assume the matching or weights are fixed, so they ignore uncertainty from estimating the propensity score and choosing variables.")
  add(RA, g_crude, "CrudeP", "p-value of CrudeRatio.", "")
  add(RA, g_n, "n_missing_outcome", "Analysis sample patients without the outcome recorded (left out of the effect model).", "If large, check outcome_max_smd_after.")
  add(RA, g_n, "n_exposed_before", "Exposed patients before matching or weighting.", "")
  add(RA, g_n, "n_controls_before", "Comparison patients before matching or weighting.", "")
  add(RA, g_n, "ess_exposed_before", "Effective sample size of the exposed group before matching or weighting.", "")
  add(RA, g_n, "ess_controls_before", "Effective sample size of the comparison group before matching or weighting.", "")
  add(RA, g_n, "ess_kept_exposed", "ess_exposed_after divided by ess_exposed_before.", "Close to 1 is best. Under 0.25 the weighting is rated poor.")
  add(RA, g_n, "ess_kept_controls", "ess_controls_after divided by ess_controls_before.", "Under 0.25 the weighting is rated poor.")
  add(RA, g_n, "match_ratio", "Matching only: number of comparison patients sought per exposed patient.", "")
  add(RA, g_n, "n_matched_sets", "Matching only: number of matched sets (one exposed patient and their controls).", "")
  add(RA, g_n, "n_controls_unused", "Matching only: comparison patients who were not matched.", "")
  add(RA, g_n, "controls_per_exposed", "Matching only: matched comparison patients per matched exposed patient.", "Below match_ratio means some exposed patients had fewer controls available.")
  add(RA, g_n, "exposed_retained_weighted", "Matching only: survey-weighted share of exposed patients that were matched.", "")
  add(RA, g_n, "controls_retained", "Matching only: share of comparison patients that were matched.", "")
  add(RA, g_n, "controls_retained_weighted", "Matching only: survey-weighted share of comparison patients that were matched.", "")
  add(RA, g_bal, "max_smd_before", "Largest absolute SMD across all assessed covariates before matching or weighting.", "Shows how different the groups were to begin with.")
  add(RA, g_bal, "balanced", "TRUE when every assessed covariate has absolute SMD below 0.1 after matching or weighting.", "")
  add(RA, g_bal, "balanced_all_assessed", "TRUE when every assessed covariate has absolute SMD below 0.1 after matching or weighting.", "")
  add(RA, g_bal, "balanced_matched_variables", "TRUE when every covariate used to build the score or weights has absolute SMD below 0.1 after.", "")
  add(RA, g_bal, "n_assessed_covariates", "Number of covariates assessed for balance.", "")
  add(RA, g_bal, "n_imbalanced_matched", "Covariates used to build the score or weights that still have absolute SMD of 0.1 or more.", "")
  add(RA, g_bal, "max_smd_matched_before", "Largest absolute SMD before, among the covariates used to build the score or weights.", "")
  add(RA, g_bal, "max_smd_matched_after", "Largest absolute SMD after, among the covariates used to build the score or weights.", "")
  add(RA, g_bal, "max_var_ratio_after", "Largest variance ratio (exposed / comparison) after matching or weighting.", "1 = same spread. Ratios far from 1 (a common rule of thumb is outside about 0.5 to 2) mean the groups differ in spread.")
  add(RA, g_bal, "min_var_ratio_after", "Smallest variance ratio after matching or weighting.", "See max_var_ratio_after.")
  add(RA, g_bal, "worst_var_ratio_fold_after", "The variance ratio furthest from 1 after, expressed as a fold (always 1 or more).", "1 = identical spread; 2 = one group has twice the variance.")
  add(RA, g_bal, "max_ks_after", "Largest Kolmogorov-Smirnov statistic after: the biggest gap between the two groups' distributions of a covariate.", "0 = identical distributions; larger = more different.")
  add(RA, g_bal, "outcome_n_imbalanced", "n_imbalanced recomputed on the patients with this outcome recorded.", "")
  add(RA, g_bal, "outcome_ess_exposed", "Effective sample size of the exposed group among patients with this outcome recorded.", "")
  add(RA, g_bal, "outcome_ess_controls", "Effective sample size of the comparison group among patients with this outcome recorded.", "")
  add(RA, g_bal, "outcome_balanced_all_assessed", "balanced_all_assessed recomputed on the patients with this outcome recorded.", "")
  add(RA, g_bal, "outcome_balance_note", "Balance verdict on the patients with this outcome recorded.", "")
  add(RA, g_bal, "matching_quality_note", "Reminder that the verdict uses rule of thumb thresholds.", "")
  add(RA, g_bal, "weighting_quality_note", "Reminder that the verdict uses rule of thumb thresholds.", "")
  add(RA, g_bal, "balancing_quality_note", "Reminder that the verdict uses rule of thumb thresholds.", "")
  add(RA, g_over, "ps_model", "Propensity score model (method and formula).", "")
  add(RA, g_over, "ps_model_coefficients", "Propensity score model coefficients: estimate (b), odds ratio (OR) and p-value for each matching variable.", "Shows which variables most strongly predict the exposure.")
  add(RA, g_over, "ps_AUC_in_sample_unweighted", "Area under the ROC curve of the propensity score for the exposure (unweighted, in sample).",
      "How well the matching variables predict who is exposed. Close to 0.5 means the groups were already similar; close to 1 means they differ a lot and overlap little. It does not measure the quality of the matching.")
  add(RA, g_over, "AUC_note", "Reminder on interpreting the AUC columns.", "")
  add(RA, g_over, "weights_model", "Entropy balancing specification (moments and formula).", "")
  add(RA, g_over, "caliper_sd", "Matching only: caliper in standard deviations of the logit propensity score.", "The largest score difference allowed for a match; smaller is stricter.")
  add(RA, g_over, "caliper_logit", "Matching only: caliper on the logit propensity score scale.", "")
  add(RA, g_over, "ps_min_exposed", "Smallest propensity score among the exposed.", "")
  add(RA, g_over, "ps_max_exposed", "Largest propensity score among the exposed.", "")
  add(RA, g_over, "ps_min_controls", "Smallest propensity score among comparison patients.", "")
  add(RA, g_over, "ps_max_controls", "Largest propensity score among comparison patients.", "The two groups' ranges should largely overlap.")
  add(RA, g_over, "ps_extreme_exposed", "Share of exposed patients with a propensity score below 0.01 or above 0.99.", "Should be near 0.")
  add(RA, g_over, "ps_extreme_controls", "Share of comparison patients with a propensity score below 0.01 or above 0.99.", "Should be near 0.")
  add(RA, g_over, "weight_min_exposed", "Smallest analysis weight among exposed patients.", "")
  add(RA, g_over, "weight_min_controls", "Smallest analysis weight among comparison patients.", "")
  add(RA, g_over, "outcome_model", "Outcome model that was fitted (class, family and formula).", "")
  add(RA, g_stat, "outcome_fit_warnings", "Warnings raised while fitting the outcome model.", "Blank = none.")
  add(RA, g_stat, "outcome_review_status", "Review notes for the outcome model.", "")
  add(RA, g_stat, "covariate_flags", "Each flagged matching variable with its reason (see the Flags sheet).", "")
  add(RA, g_stat, "has_covariate_flags", "TRUE when any matching variable is flagged.", "")
  add(RA, g_set, "target_population", "The population the effect describes.", "With matching this can be 'Retained matched exposed population', meaning only exposed patients who were matched.")
  add(RA, g_set, "ExposureVariable", "The 0/1 exposure column.", "")
  add(RA, g_set, "IndependentVariable", "IndependentVariable as written in the setup table (it may include outcome model terms).", "")
  add(RA, g_set, "Ratio", "Matching only: the Ratio from the setup table.", "")
  add(RA, g_set, "PSCaliperSD", "Matching only: the PSCaliperSD from the setup table.", "")
  add(RA, g_set, "CovariateCalipers", "Matching only: the CovariateCalipers from the setup table.", "")
  add(RA, g_set, "matching_settings_note", "Reminder that Ratio, PSCaliperSD and CovariateCalipers apply to matching only.", "")
  B <- "Balance"; gb <- "Covariate balance"
  add(B, gb, "ModelID", "The model this row belongs to (the same ModelID as in RegressionAnalysis).", "Rows with the same ModelID form one balance table; each covariate appears once per model.")
  add(B, gb, "variable", "The covariate being checked. A categorical covariate appears once per level.", "Every assessed covariate is listed, not only the matching variables.")
  add(B, gb, "in_propensity_model", "TRUE when the covariate was used to build the propensity score or the weights; FALSE when it was only assessed.",
      "FALSE covariates were not targeted, so they are the most likely to stay imbalanced. Imbalance on them is a sign the matching variables miss something.")
  add(B, gb, "mean_exposed_after", "Average of the covariate in the exposed group after matching or weighting. For a 0/1 covariate this is the proportion.", "Should be close to mean_control_after. For a 0/1 covariate, 0.45 means 45% of the group have it.")
  add(B, gb, "mean_control_after", "Average of the covariate in the comparison group after matching or weighting. For a 0/1 covariate this is the proportion.", "Should be close to mean_exposed_after.")
  add(B, gb, "smd_before", "Standardized mean difference between exposed and comparison patients before any matching or weighting: exposed minus comparison, in standard deviation units.",
      "Shows how different the groups were to begin with. 0 = identical averages. Positive = higher in the exposed group. Compare with smd_after.")
  add(B, gb, "smd_after", "Standardized mean difference after matching or weighting, in the whole analysis sample.",
      "Closer to 0 is better. An absolute value below 0.1 is the usual target. It should be much smaller than smd_before.")
  add(B, gb, "balanced", "TRUE when the absolute smd_after is below 0.1.", "FALSE means the groups still differ on this covariate, so the effect may partly reflect that difference. Also FALSE when the SMD could not be computed.")
  add(B, gb, "smd_after_outcome_sample", "Same as smd_after, recomputed using only the patients who have this outcome recorded.",
      "Patients with a missing outcome are left out of the effect estimate, which can upset the balance. Use this column rather than smd_after to judge balance for this particular outcome. Blank = could not be computed (see outcome_sample_note).")
  add(B, gb, "balanced_outcome_sample", "TRUE when the absolute smd_after_outcome_sample is below 0.1.", "FALSE means the covariate is imbalanced among the patients that actually contribute to this outcome's effect. Blank = could not be computed.")
  add(B, gb, "outcome_sample_note", "Why balance among the patients with this outcome recorded could not be computed.", "Appears only when at least one row could not be computed.")
  gd <- "Covariate balance (detail = TRUE)"
  add(B, gd, "independent", "The exposure (0/1).", "")
  add(B, gd, "MatchingVariable", "Covariate(s) the score or weights were built from.", "")
  add(B, gd, "Method", "glm, svyglm or lm.", "")
  add(B, gd, "mean_control_before", "Average of the covariate in the comparison group before matching or weighting (proportion for 0/1 covariates).", "")
  add(B, gd, "mean_exposed_before", "Average of the covariate in the exposed group before matching or weighting (proportion for 0/1 covariates).", "")
  add(B, gd, "var_ratio_before", "Variance ratio, exposed / comparison, before (continuous covariates).", "1 = same spread.")
  add(B, gd, "var_ratio_after", "Variance ratio after.", "Close to 1 = similar spread. A common rule of thumb is to worry outside about 0.5 to 2.")
  add(B, gd, "ks_before", "Kolmogorov-Smirnov statistic before: the largest gap between the two groups' distributions.", "0 = identical distributions; larger = more different.")
  add(B, gd, "ks_after", "Kolmogorov-Smirnov statistic after.", "Should be smaller than ks_before.")
  add(B, gd, "exact_match", "TRUE when exact matching was required on the covariate.", "")
  gg <- "Group diagnostics (detail = TRUE)"
  group_metric <- list(
    n_before = c("Patients before matching or weighting.", ""),
    n_after = c("Patients in the analysis sample.", ""),
    n_unmatched = c("Patients left out of the analysis sample.", ""),
    retained_fraction = c("Share of patients kept in the analysis sample.", "Close to 1 is best."),
    weighted_retained_fraction = c("Survey-weighted share of patients kept.", ""),
    ess_before = c("Effective sample size before.", ""),
    ess_after = c("Effective sample size after.", "Much smaller than n_after means a few patients carry large weights."),
    weight_min = c("Smallest analysis weight.", ""),
    weight_max = c("Largest analysis weight.", "A very large maximum means a few patients dominate."),
    ps_min = c("Smallest propensity score.", ""), ps_max = c("Largest propensity score.", ""),
    fraction_ps_extreme = c("Share of propensity scores below 0.01 or above 0.99.", "Should be near 0."),
    outside_common_support = c("Patients outside the propensity score range shared by both groups.", "Should be small."),
    max_abs_smd_after = c("Largest absolute SMD after, across assessed covariates (the same in both group columns).", "Below 0.1 is the usual target."),
    n_balance_flags = c("Covariates with absolute SMD of 0.1 or more after (the same in both group columns).", "0 is the goal."),
    n_clipped_scores = c("Propensity scores at 0 or 1 (the same in both group columns).", "Should be 0."),
    fit_warnings = c("Fitting warnings (the same in both group columns).", "Blank = none."),
    review_status = c("Overall balance verdict (the same in both group columns).", ""))
  for (g in c("exposed", "unexposed")) for (m in names(group_metric))
    add("Group", gg, paste0(g, "_", m), paste(if (g == "exposed") "Exposed group:" else "Comparison group:", group_metric[[m]][1]), group_metric[[m]][2])
  O <- "OutcomeBalance"; go <- "Balance on this outcome's patients (detail = TRUE)"
  for (cn in c("dependent", "independent", "MatchingVariable", "Method", "Type"))
    add(O, go, cn, paste("Same as the", cn, "column in RegressionAnalysis."), "")
  add(O, go, "ModelID", "The model this row belongs to (the same ModelID as in RegressionAnalysis).", "")
  add(O, go, "variable", "The covariate being checked. A categorical covariate appears once per level.", "")
  for (m in c("smd_before", "smd_after", "var_ratio_before", "var_ratio_after", "ks_before", "ks_after", "mean_control_before",
              "mean_exposed_before", "mean_control_after", "mean_exposed_after", "balanced")) {
    src <- Filter(function(r) r$Column == m, rows[vapply(rows, function(r) r$Sheet %in% c(B), logical(1))])[[1]]
    add(O, go, m, paste("On the patients who have this outcome recorded:", src$Meaning), src$HowToRead)
  }
  add(O, go, "diagnostic_error", "Why balance could not be computed on this outcome's patients.", "Blank = computed.")
  Fl <- "Flags"; gf <- "Pre-fit screen of candidate covariates"
  add(Fl, gf, "variable", "The candidate covariate that was screened.", "")
  add(Fl, gf, "missing_n", "Records with a missing or non-finite value.", "Many missing values reduce the usable sample.")
  add(Fl, gf, "constant", "TRUE when every non-missing record has the same value.", "A constant covariate carries no information and cannot be balanced.")
  add(Fl, gf, "binary", "TRUE for 0/1 covariates.", "")
  add(Fl, gf, "correlation", "Correlation between the covariate and the exposure (survey-weighted). For 0/1 covariates this is the phi coefficient.",
      "Values near +1 or -1 mean the covariate almost is the exposure. The screen marks an absolute value of 0.8 or more as a near duplicate.")
  add(Fl, gf, "n_with", "0/1 covariates: number of records with the covariate = 1.", "")
  add(Fl, gf, "exposed_share_with", "0/1 covariates: share of patients who are exposed among those with the covariate.", "A share near 0 or 1 means almost everyone with the covariate is in one exposure group, so there are no comparable patients in the other group.")
  add(Fl, gf, "exposed_share_without", "0/1 covariates: share of patients who are exposed among those without the covariate.", "See exposed_share_with.")
  add(Fl, gf, "near_duplicate", "TRUE when the absolute correlation with the exposure is 0.8 or more (default threshold).", "The covariate nearly is the exposure; propensity scores will sit near 0 or 1.")
  add(Fl, gf, "extreme_level", "TRUE when a level with at least 20 records is almost entirely one exposure group (within 1% of all or none; with ATT only the 'all exposed' end counts). Default thresholds.", "There are no comparable patients in the other group for that level.")
  add(Fl, gf, "post_exposure", "TRUE when you listed the covariate as measured after the exposure.",
      "Adjusting for something the exposure causes can hide part of the real effect, so such a covariate should usually not be matched or weighted on. The data cannot detect this; only your list sets it.")
  add(Fl, gf, "action", "Suggested review. 'keep' means no concern; otherwise the reason for concern.", "Descriptive only. Nothing is removed automatically.")
  add(Fl, gf, "joint_duplicate", "TRUE when this 0/1 covariate, combined with others, nearly reproduces the exposure.", "For example several drug or alcohol codes that together define a substance use flag.")
  add(Fl, gf, "unsupported_levels", "Levels found in one exposure group only.", "Patients in such a level have no comparable patients in the other group.")
  add(Fl, gf, "outside_control_range", "Numeric covariates: exposed patients whose value is beyond the lowest or highest value among comparison patients.", "These patients cannot be matched on this covariate.")
  add(Fl, gf, "eb_mean_outside_support", "TRUE when the average entropy balancing must reach is outside the range of the comparison patients.", "Entropy balancing cannot reach that target.")
  add(Fl, gf, "flagged", "TRUE when action is not 'keep'.", "")
  add(Fl, gf, "exposure", "The exposure the check refers to.", "")
  add(Fl, gf, "estimand", "The estimand the check refers to (ATT or ATE).", "")
  do.call(rbind, rows)
}

#' Glossary of statistical terms, used by the README
#'
#' @return A data frame: Term, Meaning.
#' @keywords internal
#' @noRd
Matching_Glossary <- function() {
  g <- function(term, meaning) data.frame(Term = term, Meaning = meaning, stringsAsFactors = FALSE)
  do.call(rbind, list(
    g("Exposure", "The 0/1 variable whose effect is being studied. 1 = exposed, 0 = comparison (control) patients."),
    g("Outcome", "The result being compared between the exposed and comparison patients (the dependent variable)."),
    g("Covariate", "A patient characteristic (such as age or a comorbidity) that could differ between the groups and affect the outcome."),
    g("Confounding", "When a covariate affects both who is exposed and the outcome, so a raw comparison mixes its effect with the exposure's. Matching and weighting aim to remove it for the measured covariates."),
    g("Propensity score", "A patient's estimated probability of being exposed, given the matching or weighting variables. Patients with similar scores are comparable."),
    g("Matching", "Pairing each exposed patient with unexposed patients who have a similar propensity score, then comparing the matched patients only."),
    g("Caliper", "The largest propensity score difference allowed for a match. Exposed patients with no comparison patient inside the caliper are left unmatched."),
    g("Weighting (IPTW)", "Giving each patient a weight based on the probability of the exposure they actually had, so the weighted groups resemble each other. No patient is dropped."),
    g("Entropy balancing", "Choosing weights directly so the weighted groups have equal covariate averages (and, optionally, variances) without fitting a propensity score."),
    g("Analysis sample", "The patients the effect is estimated from: the matched patients for matching, all patients with their weights for weighting and entropy balancing."),
    g("ATT, ATE, ATO", "The population the effect describes. ATT: patients who were exposed. ATE: all patients. ATO: patients with the most overlap between the groups (those who could plausibly have been in either group)."),
    g("Standardized mean difference (SMD)", "The difference between the group averages of a covariate divided by its standard deviation, so covariates on different scales can be compared. 0 = identical averages; an absolute value below 0.1 is the usual target for balance."),
    g("Balance", "How alike the exposed and comparison groups are on a covariate. Good balance after matching or weighting makes the comparison fairer."),
    g("Variance ratio", "The spread of a covariate in the exposed group divided by its spread in the comparison group. 1 = same spread."),
    g("Kolmogorov-Smirnov (KS) statistic", "The largest gap between the two groups' cumulative distributions of a covariate. 0 = identical distributions."),
    g("Effective sample size (ESS)", "The number of equally weighted patients that would carry the same information as the weighted group. A much smaller ESS than the real count means a few patients dominate."),
    g("Common support (overlap)", "The range of propensity scores found in both groups. Patients outside it have no comparable patients in the other group."),
    g("Odds ratio", "The odds of the outcome in the exposed group divided by the odds in the comparison group. 1 = no difference."),
    g("Risk ratio / risk difference", "The risk (proportion with the outcome) in the exposed group divided by, or minus, the risk in the comparison group. 1 (ratio) or 0 (difference) = no difference."),
    g("Confidence interval", "The range of effects that is compatible with the data at the stated level (95%). An interval that includes 1 (ratio) or 0 (difference) is compatible with no effect."),
    g("p-value", "How unusual the result would be if the exposure had no effect. Small values (below 0.05 by convention) are evidence against no effect. It does not measure the size or importance of the effect."),
    g("Bonferroni adjustment", "A strict correction for testing many models: each p-value is multiplied by the number of tests."),
    g("False discovery rate (FDR)", "A less strict correction (Benjamini-Hochberg) that controls the expected share of false findings among the results called significant."),
    g("Crude (unadjusted) estimate", "The same comparison with no matching or weighting. Comparing it with the adjusted estimate shows how much the adjustment mattered.")))
}

#' "Start here" page of the README
#'
#' @param prefix \code{"PS"}, \code{"IPTW"} or \code{"EB"}.
#' @param detail Logical; whether the detailed files were written.
#' @return A data frame: Topic, Explanation.
#' @keywords internal
#' @noRd
Matching_StartHere <- function(prefix, detail = FALSE) {
  what <- switch(prefix,
    PS = "Each file asks one question about one outcome: does the exposure change the outcome once exposed and comparison patients are made similar? Here that is done by propensity score matching. Every exposed patient is paired with unexposed patients who had a similar probability (propensity score) of being exposed, based on the matching variables, and the effect is estimated in the matched patients only (the 'analysis sample').",
    IPTW = "Each file asks one question about one outcome: does the exposure change the outcome once exposed and comparison patients are made similar? Here that is done by inverse probability of treatment weighting. Every patient receives a weight based on the probability of the exposure they actually had, so the weighted exposed and comparison groups resemble each other on the weighting variables. All patients stay in the analysis (the 'analysis sample'), but some count more than others.",
    EB = "Each file asks one question about one outcome: does the exposure change the outcome once exposed and comparison patients are made similar? Here that is done by entropy balancing. The weights are chosen so the weighted groups have the same averages (and, with moments = 2, the same variances) on the balancing variables. No propensity score is used. All patients stay in the analysis (the 'analysis sample'), but some count more than others.",
    "Each file compares an outcome between exposed and comparison patients after adjusting for measured covariates.")
  verdict <- switch(prefix,
    PS = "matching_quality is 'poor' when fewer than 80% of exposed patients were matched or any propensity score is exactly 0 or 1, 'fair' when any covariate still has an absolute SMD of 0.1 or more, and otherwise 'SMD/retention thresholds passed'.",
    IPTW = "weighting_quality is 'poor' when either group keeps less than 25% of its effective sample size or any propensity score is exactly 0 or 1, 'fair' when any covariate still has an absolute SMD of 0.1 or more, and otherwise 'SMD/ESS thresholds passed'.",
    EB = "balancing_quality is 'poor' when either group keeps less than 25% of its effective sample size or a balanced covariate still has an absolute SMD of 0.01 or more, 'fair' when any other covariate has an absolute SMD of 0.1 or more, and otherwise 'SMD/ESS thresholds passed'.",
    "The quality column gives an automatic verdict ('poor', 'fair' or 'thresholds passed').")
  quality_col <- switch(prefix, PS = "matching_quality", IPTW = "weighting_quality", EB = "balancing_quality", "the quality column")
  sheets <- paste0(
    "RegressionAnalysis: one row per model (one ModelID). The estimated effect, how precise it is, how many patients were used, and an automatic quality verdict.\n",
    "Balance: one row per covariate per model. Shows whether exposed and comparison patients look alike before and after the adjustment. This is how you check that the adjustment worked.\n",
    if (isTRUE(detail)) "OutcomeBalance: the same balance table recomputed on only the patients who have this outcome recorded.\n" else "",
    "Flags: one row per candidate covariate, screened before fitting for problems such as a covariate that almost is the exposure. Empty when no flags were supplied. Flags describe concerns; nothing is removed automatically.",
    if (isTRUE(detail)) "" else "\n\nThe files were written with detail = FALSE (compact). Rerun with detail = TRUE to add every per-covariate column to RegressionAnalysis and to write the full Balance and OutcomeBalance sheets.")
  data.frame(Topic = c("What these files are", "Files", "Sheets in each workbook", "How to read a result (suggested order)",
                       "Rules of thumb used", "Many models at once", "What the results can and cannot tell you", "Where to find definitions"),
             Explanation = c(
    what,
    paste0("One workbook per outcome, named ", prefix, "_Sub_DataPredictions_<outcome>.xlsx, plus this guide. The exposure is coded 0/1: 1 = exposed, 0 = comparison (control) patients. A ratio above 1, or a positive difference, means the outcome is higher in the exposed group."),
    sheets,
    paste0(
      if (isTRUE(detail)) "1. Check status and error in RegressionAnalysis. status = ok means the model was fitted; failed means no effect was estimated and error says why.\n"
      else "1. Check that Estimate is filled in. A blank Estimate means the model could not be fitted; rerun with detail = TRUE to see the reason (status and error columns).\n",
      "2. Read ", quality_col, ". 'poor': do not rely on the effect until you know why. 'fair': the adjustment worked except for some covariates (see step 3). 'thresholds passed': the automatic checks were met.\n",
      "3. Open the Balance sheet. For every covariate compare smd_before with smd_after: after the adjustment smd_after should be close to 0, and balanced should be TRUE. Use smd_after_outcome_sample and balanced_outcome_sample to judge balance among the patients who actually contribute to this outcome. A covariate with balanced = FALSE still differs between the groups, and the effect may partly reflect that.\n",
      "4. Read the effect. Yes/no outcome: OddsRatio (also RiskRatio and RiskDifference). Count outcome: ExpectedCountRatio. Continuous outcome: MeanDifference. Then look at the confidence limits: if the interval includes 1 (ratios) or 0 (differences) the data are compatible with no effect. For OddsRatio and the count ratios, exp() of '2.5 %' and '97.5 %' gives the limits.\n",
      "5. Judge significance with Pr(>|t|), or with the adjusted columns (FDR, p_fdr_overall) when many models were run.\n",
      if (isTRUE(detail)) "6. Read flagged_covariates and the Flags sheet as cautions about the matching variables."
      else "6. Read the Flags sheet as cautions about the matching variables."),
    paste0(
      "Balance: an absolute standardized mean difference (SMD) below 0.1 counts as balanced. This is applied automatically in the 'balanced' columns.\n",
      "Significance: p below 0.05 is marked significant (Significance, sig_unadjusted); adjusted p-values use the same cut-off.\n",
      "Quality verdict: ", verdict, "\n",
      "These are conventions, not proof. Judge them together with the size of the effect and its confidence limits."),
    "With many models, some will look significant by chance. Use the adjusted p-values: FDR (p_fdr) within one outcome, or p_fdr_overall across all outcomes in the analysis. Bonferroni (p.adjust, p_bonferroni) is stricter.",
    "The effect compares exposed and comparison patients who look alike on the measured variables. It cannot account for differences that were not measured, and it describes an association, not proof of cause. The confidence limits treat the matching or weights as fixed, so they leave out the uncertainty from estimating them and from choosing the variables. A statistically significant effect can be too small to matter, and a non-significant effect can be large but imprecise: read the size of the effect and its confidence limits, not only the p-value.",
    "Columns sheet: every column in the workbooks, what it is and how to read it (filter by Sheet or Group). Glossary sheet: the statistical terms used."),
    stringsAsFactors = FALSE)
}

#' Write <prefix>_README.xlsx
#'
#' A "Start here" page, a "Columns" sheet defining every column present and a
#' "Glossary".
#'
#' @param res,balance,outcome_balance,flags Results and sheets written.
#' @param Fpath Output folder.
#' @param prefix \code{"PS"}, \code{"IPTW"} or \code{"EB"}.
#' @param detail Logical; whether the detailed sheets were written.
#' @return Invisibly, the column definitions.
#' @keywords internal
#' @noRd
Matching_WriteReadme <- function(res, balance, outcome_balance, flags, Fpath, prefix, detail = FALSE) {
  guide <- Matching_ColumnGuide()
  sheet_columns <- list(RegressionAnalysis = names(res))
  if (isTRUE(detail)) {
    sheet_columns$Balance <- if (is.null(balance)) character(0) else names(balance)
    sheet_columns$OutcomeBalance <- if (is.null(outcome_balance) || !ncol(outcome_balance)) character(0) else
      setdiff(names(outcome_balance), c("full_smd_after", "full_mean_control_after", "full_mean_exposed_after", "full_balanced", "in_propensity_model"))
  } else {
    res_df <- as.data.frame(res)
    slim <- unlist(lapply(unique(res_df$dependent), function(dep) names(Matching_BalanceSheet(res_df[res_df$dependent == dep, , drop = FALSE], balance, outcome_balance, dep))))
    sheet_columns$Balance <- if (is.null(slim)) character(0) else unique(slim)
  }
  sheet_columns$Flags <- if (is.null(flags)) character(0) else names(flags)
  balance_metrics <- c("mean_control_before", "mean_exposed_before", "mean_control_after", "mean_exposed_after", "smd_before", "smd_after",
                       "var_ratio_before", "var_ratio_after", "ks_before", "ks_after", "balanced", "in_propensity_model")
  lookup <- function(sheet, column) {
    hit <- guide[guide$Sheet == sheet & guide$Column == column, , drop = FALSE]
    if (nrow(hit)) return(hit[1, c("Group", "Meaning", "HowToRead")])
    if (sheet == "Balance" && grepl("^(exposed|unexposed)_", column)) {
      hit <- guide[guide$Sheet == "Group" & guide$Column == column, , drop = FALSE]
      if (nrow(hit)) return(hit[1, c("Group", "Meaning", "HowToRead")])
    }
    if (sheet == "RegressionAnalysis") {
      for (metric in balance_metrics) {
        outcome_prefix <- paste0("outcome_", metric, "_")
        plain_prefix <- paste0(metric, "_")
        for (p in c(outcome_prefix, plain_prefix)) if (startsWith(column, p)) {
          src <- guide[guide$Sheet == "Balance" & guide$Column == metric, , drop = FALSE]
          if (!nrow(src)) next
          scope <- if (p == outcome_prefix) " On the patients with this outcome recorded." else ""
          return(data.frame(Group = "Per-covariate balance (detail = TRUE)",
                            Meaning = paste0(src$Meaning[1], scope, " Covariate: ", sub(paste0("^", p), "", column), "."),
                            HowToRead = src$HowToRead[1], stringsAsFactors = FALSE))
        }
      }
    }
    data.frame(Group = "", Meaning = "", HowToRead = "", stringsAsFactors = FALSE)
  }
  columns <- do.call(rbind, lapply(names(sheet_columns), function(s) {
    cols <- sheet_columns[[s]]
    if (!length(cols)) return(NULL)
    info <- do.call(rbind, lapply(cols, function(cn) lookup(s, cn)))
    data.frame(Sheet = s, Group = info$Group, Column = cols, `What it is` = info$Meaning, `How to read it` = info$HowToRead,
               check.names = FALSE, stringsAsFactors = FALSE, row.names = NULL)
  }))
  start <- Matching_StartHere(prefix, detail)
  glossary <- Matching_Glossary()
  tryCatch({
    wb <- openxlsx::createWorkbook()
    header <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#DDEBF7", wrapText = TRUE, valign = "top", border = "bottom")
    body <- openxlsx::createStyle(wrapText = TRUE, valign = "top")
    bold <- openxlsx::createStyle(textDecoration = "bold", wrapText = TRUE, valign = "top")
    put <- function(name, df, widths, bold_first = FALSE) {
      openxlsx::addWorksheet(wb, name)
      openxlsx::writeData(wb, name, df, headerStyle = header)
      openxlsx::setColWidths(wb, name, cols = seq_along(widths), widths = widths)
      openxlsx::addStyle(wb, name, body, rows = seq_len(nrow(df)) + 1L, cols = seq_len(ncol(df)), gridExpand = TRUE)
      if (bold_first) openxlsx::addStyle(wb, name, bold, rows = seq_len(nrow(df)) + 1L, cols = 1L, gridExpand = TRUE)
      openxlsx::freezePane(wb, name, firstRow = TRUE)
    }
    put("Start here", start, c(30, 130), bold_first = TRUE)
    put("Columns", columns, c(20, 34, 32, 70, 90))
    openxlsx::addFilter(wb, "Columns", row = 1, cols = seq_len(ncol(columns)))
    put("Glossary", glossary, c(36, 120), bold_first = TRUE)
    openxlsx::saveWorkbook(wb, file.path(Fpath, paste0(prefix, "_README.xlsx")), overwrite = TRUE)
  }, error = function(e) message(sprintf("  Could not write the README: %s", conditionMessage(e))))
  invisible(columns)
}

#' Write one workbook per dependent and the README
#'
#' @param res Compiled results.
#' @param balance,outcome_balance,flags Sheets to write.
#' @param Fpath Output folder.
#' @param prefix \code{"PS"}, \code{"IPTW"} or \code{"EB"}.
#' @param detail Logical; write the full Balance and OutcomeBalance sheets.
#' @return Called for its side effect.
#' @keywords internal
#' @noRd
Matching_WriteResults <- function(res, balance, outcome_balance, flags, Fpath, prefix, detail = FALSE) {
  dir.create(Fpath, recursive = TRUE, showWarnings = FALSE)
  extra_outcome_columns <- c("full_smd_after", "full_mean_control_after", "full_mean_exposed_after", "full_balanced", "in_propensity_model")
  for (dep in unique(res$dependent)) {
    sub <- res[res$dependent == dep]
    sheets <- list(RegressionAnalysis = sub)
    if (isTRUE(detail)) {
      sheets$Balance <- if (is.null(balance)) data.frame() else balance[paste(balance$independent, balance$MatchingVariable, balance$Method) %in%
                                                                          paste(sub$independent, sub$MatchingVariable, sub$Method), ]
      sheets$OutcomeBalance <- if (is.null(outcome_balance) || !nrow(outcome_balance)) data.frame() else {
        ob <- as.data.frame(outcome_balance[outcome_balance$dependent == dep, ])
        ob[, setdiff(names(ob), extra_outcome_columns), drop = FALSE] }
    } else {
      sheets$Balance <- Matching_BalanceSheet(sub, balance, outcome_balance, dep)
    }
    sheets$Flags <- if (is.null(flags)) data.frame() else flags
    tryCatch(openxlsx::write.xlsx(sheets, file = file.path(Fpath, paste0(prefix, "_Sub_DataPredictions_", dep, ".xlsx")), colWidths = "auto"),
             error = function(e) message(sprintf("  Could not write the file for '%s': %s", dep, conditionMessage(e))))
  }
  Matching_WriteReadme(res, balance, outcome_balance, flags, Fpath, prefix, detail)   # <prefix>_README.xlsx: what each column means
}

#' Checkpoint folder
#'
#' @param Fpath Output folder.
#' @param prefix File prefix.
#' @param checkpoint Logical; whether checkpoints are used.
#' @return The folder, or \code{NULL} when checkpoints are off or there is no
#'   \code{Fpath}.
#' @keywords internal
#' @noRd
Matching_CheckpointDir <- function(Fpath, prefix, checkpoint) {
  if (!isTRUE(checkpoint) || is.null(Fpath)) return(NULL)
  dir <- file.path(Fpath, paste0(prefix, "_checkpoints"))
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  dir
}

#' What a saved result depends on
#'
#' The setup rows, the settings and each used column's length, missing count,
#' distinct values and sum.
#'
#' @param DT Data frame.
#' @param rows Setup rows.
#' @param settings List of settings.
#' @param columns Columns used.
#' @return A list.
#' @keywords internal
#' @noRd
Matching_Fingerprint <- function(DT, rows, settings, columns) {
  columns <- unique(columns)
  if (!all(columns %in% names(DT))) stop("Fingerprint columns missing from data.")
  data <- as.data.frame(DT)[, columns, drop=FALSE]
  rownames(data) <- NULL
  rows <- as.data.frame(rows); rownames(rows) <- NULL
  if (!requireNamespace("digest", quietly=TRUE)) stop("Package digest is required for checkpoint content hashing.")
  list(version=2L, rows=rows, settings=settings,
       data=digest::digest(data, algo="sha256", serialize=TRUE))
}

#' Load a checkpoint when its key and fingerprint match
#'
#' @param dir Checkpoint folder.
#' @param n Chunk number.
#' @param key Chunk key.
#' @param fingerprint \code{Matching_Fingerprint()} result.
#' @return The saved list, or \code{NULL}.
#' @keywords internal
#' @noRd
Matching_CheckpointLoad <- function(dir, n, key, fingerprint) {
  if (is.null(dir)) return(NULL)
  file <- file.path(dir, sprintf("chunk_%05d.rds", n))
  if (!file.exists(file)) return(NULL)
  saved <- tryCatch(readRDS(file), error = function(e) NULL)
  if (is.null(saved) || !identical(saved$key, key) || !identical(saved$fingerprint, fingerprint)) return(NULL)
  saved
}

#' Save a checkpoint
#'
#' Writes to a temporary file and renames it, so a stopped run never leaves a
#' broken checkpoint.
#'
#' @param dir Checkpoint folder.
#' @param n Chunk number.
#' @param key Chunk key.
#' @param fingerprint \code{Matching_Fingerprint()} result.
#' @param rows,balance,outcome_balance Results to save.
#' @return Invisibly, the file path.
#' @keywords internal
#' @noRd
Matching_CheckpointSave <- function(dir, n, key, fingerprint, rows, balance, outcome_balance) {
  if (is.null(dir)) return(invisible(NULL))
  file <- file.path(dir, sprintf("chunk_%05d.rds", n))
  tmp <- paste0(file, ".tmp")   # write then rename, so a run stopped mid-save never leaves a broken checkpoint
  saveRDS(list(key = key, fingerprint = fingerprint, rows = rows, balance = balance, outcome_balance = outcome_balance), tmp)
  if (file.exists(file)) unlink(file)   # file.rename does not overwrite on Windows
  if (!file.rename(tmp, file)) stop("Could not finalize checkpoint: ", file)
  invisible(file)
}

#' Remove the checkpoint folder
#'
#' @param dir Checkpoint folder, or \code{NULL}.
#' @return \code{NULL}, invisibly.
#' @keywords internal
#' @noRd
Matching_CheckpointClear <- function(dir) {
  if (!is.null(dir)) unlink(dir, recursive = TRUE)
  invisible(NULL)
}
