#' Propensity score matched regression driven by a setup table
#'
#' The matched counterpart of \code{RegressionCompile()}. For each exposure x
#' matching specification x method in \code{setup}: a propensity score from the
#' matching covariates (\code{PSFit()}), 1:\code{ratio} nearest-neighbour
#' matching on its logit within the caliper (\code{PSMatchWeight()}), then each
#' row's outcome model on the matched sample (\code{PSOutcomeCompile()}).
#'
#' \code{setup} columns: \code{DependentVariable}; \code{IndependentVariable}
#' (the exposure, plus optional outcome-model terms); \code{MatchingVariable}
#' (covariates joined with \code{+}); optional \code{Method} (\code{"svyglm"},
#' \code{"glm"} or \code{"lm"}), \code{Type} (\code{"Logistic"},
#' \code{"Linear"} or \code{"Poisson"}), \code{ModelID},
#' \code{ExposureVariable}, \code{Ratio}, \code{PSCaliperSD} and
#' \code{CovariateCalipers}.
#'
#' Results follow \code{RegressionCompile()}'s columns (Estimate, Std. Error, t
#' value, Pr(>|t|), odds or count ratio, confidence limits, Bonferroni and FDR
#' within each dependent and overall) followed by the effect measures on their
#' natural scale and a short set of quality checks (sample size, effective
#' sample size, overlap, largest weight, largest SMD after matching, and the
#' matching quality verdict). Compact results retain \code{status},
#' \code{error}, and the setup fields needed to identify failed or repeated
#' configurations. With \code{Fpath}, one
#' \code{PS_Sub_DataPredictions_<dependent>.xlsx} per dependent (sheets
#' RegressionAnalysis, Balance and Flags) and \code{PS_README.xlsx} are
#' written.
#'
#' @param DT Data frame.
#' @param setup Setup table, for example read from
#'   \code{Matching_RegressionSetup.xlsx}.
#' @param weights,cluster,strata Design columns, or \code{NULL}.
#' @param method Method for rows without one.
#' @param type Type for rows without one (\code{NULL}: from the data).
#' @param ratio Controls per exposed record for rows without \code{Ratio}.
#' @param Fpath Optional output folder.
#' @param balance_covariates Covariates whose balance is assessed in addition to
#'   the matching covariates.
#' @param flags Optional \code{MatchingAssumptionCheck()$table}.
#' @param caliper Caliper in SDs of the logit score for rows without
#'   \code{PSCaliperSD}.
#' @param checkpoint Logical; with \code{Fpath}, save each matching as it
#'   finishes and resume a stopped run.
#' @param detail Logical; \code{TRUE} returns every column and writes the full
#'   Balance and OutcomeBalance sheets.
#' @return A data.table with one row per setup row, with attributes
#'   \code{balance}, \code{outcome_balance} and \code{flags}.
#' @seealso \code{\link{IPTWRegression}}, \code{\link{EBRegression}},
#'   \code{\link{MatchedComparisonOverviewBubblePlot}}
#' @examples
#' \dontrun{
#' setup <- MatchingRegressionSetup(model_dt, c("Death_01", "LOS"), "sud_any_01",
#'   matching = c("AGE", "sex_01"), type = c(LOS = "Linear"))
#' res <- PSRegression(model_dt, setup, weights = "DISCWT",
#'   cluster = "cluster_id", strata = "strata_id", Fpath = tempdir())
#' }
#' @export
PSRegression <- function(DT, setup, weights = NULL, cluster = NULL, strata = NULL, method = "glm", type = NULL, ratio = 3, Fpath = NULL, balance_covariates = NULL, flags = NULL, caliper = 0.2,
                         checkpoint = TRUE, detail = FALSE) {
  if (!is.numeric(caliper) || length(caliper) != 1L || !is.finite(caliper) || caliper <= 0)
    stop("caliper must be one positive finite number, in SD units of the logit propensity score.")
  if (!is.numeric(ratio) || length(ratio) != 1L || !is.finite(ratio) || ratio < 1 || ratio != floor(ratio))
    stop("ratio must be a positive whole number of controls per exposed patient.")
  DT <- as.data.frame(DT)
  prep <- PS_PrepareSetup(DT, setup, weights, cluster, strata, method, type, balance_covariates, flags, verb = "matched")
  setup <- prep$setup; balance_covariates <- prep$balance_covariates
  for (nm in c("Ratio", "PSCaliperSD", "CovariateCalipers")) if (!nm %in% names(setup)) setup[[nm]] <- NA
  resolve <- function(x, default, integer=FALSE) {
    x[PS_Optional(x)] <- default
    x <- suppressWarnings(as.numeric(x))
    if (any(!is.finite(x) | x <= 0) || (integer && any(x != floor(x)))) stop("Invalid Ratio or PSCaliperSD row override.")
    x
  }
  setup$Ratio <- resolve(setup$Ratio, ratio, TRUE)
  setup$PSCaliperSD <- resolve(setup$PSCaliperSD, caliper)
  setup$CovariateCalipers[PS_Optional(setup$CovariateCalipers)] <- "NA"
  parsed_calipers <- lapply(setup$CovariateCalipers, PS_ParseCalipers, DT=DT)
  balance_covariates <- unique(c(balance_covariates, unlist(lapply(parsed_calipers,names)),
    unlist(lapply(setup$.OutcomeAdjust,PS_BaseVars))))
  crude <- PS_CrudeEffects(DT, setup, weights, cluster, strata, estimand = "ATT")
  key <- paste(setup$IndependentVariable, setup$MatchingVariable, setup$Method, setup$.OutcomeAdjust, setup$Ratio, setup$PSCaliperSD, setup$CovariateCalipers, sep = " | ")
  cp <- Matching_CheckpointDir(Fpath, "PS", checkpoint)
  settings <- list(analysis = "PS", weights = weights, cluster = cluster, strata = strata, balance_covariates = balance_covariates)
  rows <- list(); balance <- list(); outcome_balance <- list()
  for (k in unique(key)) {
    idx <- which(key == k); ind <- setup$IndependentVariable[idx[1]]; covs <- setdiff(PS_SplitVars(setup$MatchingVariable[idx[1]]), ind)
    ratio <- setup$Ratio[idx[1]]; caliper <- setup$PSCaliperSD[idx[1]]
    row_calipers <- parsed_calipers[[idx[1]]]
    adjust <- if (nzchar(setup$.OutcomeAdjust[idx[1]])) PS_SplitVars(setup$.OutcomeAdjust[idx[1]]) else NULL
    meth <- setup$Method[idx[1]]; svy <- meth == "svyglm"; ps_method <- if (svy) "svyglm" else "glm"
    n_k <- match(k, unique(key))
    fp <- Matching_Fingerprint(DT, setup[idx, ], settings, c(setup$DependentVariable[idx], ind, PS_BaseVars(covs),
                                                            PS_BaseVars(setup$.OutcomeAdjust[idx]), weights, cluster, strata, balance_covariates))
    saved <- Matching_CheckpointLoad(cp, n_k, k, fp)
    if (!is.null(saved)) {
      message(sprintf("[%d / %d] %s (from checkpoint)", n_k, length(unique(key)), k))
      rows[[k]] <- saved$rows; balance[[k]] <- saved$balance; outcome_balance[names(saved$outcome_balance)] <- saved$outcome_balance
      next
    }
    message(sprintf("[%d / %d] %s", n_k, length(unique(key)), k))
    out <- tryCatch({
      fit <- PSFit(DT, ind, covs, method = ps_method, weights = weights, cluster = if (svy) cluster else NULL, strata = if (svy) strata else NULL)
      m <- PSMatchWeight(fit, approach = "match", estimand = "ATT", ratio = ratio, caliper = caliper,
                         balance_covariates = setdiff(unique(c(balance_covariates, PS_BaseVars(covs))), ind), covariate_calipers=row_calipers)
      eff <- lapply(setNames(nm = unique(setup$Type[idx])), function(tp)
        PSOutcomeCompile(list(matched = m), unique(setup$DependentVariable[idx][setup$Type[idx] == tp]), family = PS_TypeFamily[[tp]], adjust=adjust, PrintSummary = FALSE))
      list(fit = fit, m = m, eff = eff)
    }, error = function(err) err)
    if (inherits(out, "error")) { rows[[k]] <- PS_FailedRows(setup, idx, conditionMessage(out)); next }
    m <- out$m; fit <- out$fit; dg <- m$diagnostics; ex <- dg$group == "exposed"; un <- dg$group == "unexposed"
    balance[[k]] <- data.frame(independent = ind, MatchingVariable = setup$MatchingVariable[idx[1]], Method = meth, m$balance, row.names = NULL)
    retained <- 1 - m$n_unmatched / sum(DT[[ind]] == 1)
    info <- PS_PropensitySummary(DT, fit, m, ind, ps_method, quality_label = "matching", matching = list(
      covariate_calipers = setup$CovariateCalipers[idx[1]], match_ratio = ratio, caliper_sd = caliper, caliper_logit = caliper * sd(fit$logit_ps),
      n_matched_sets = length(unique(m$data$subclass)), n_controls_unused = dg$n_unmatched[un],
      controls_per_exposed = sum(m$data[[ind]] == 0) / sum(m$data[[ind]] == 1),
      exposed_retained_weighted = dg$weighted_retained_fraction[ex],
      controls_retained = dg$retained_fraction[un], controls_retained_weighted = dg$weighted_retained_fraction[un]))
    smd <- m$balance
    quality <- if (retained < 0.8 || fit$n_clipped > 0) {
      sprintf("poor: %.0f%% of exposed matched, %d scores at bounds; check the matching variable(s)", 100 * retained, fit$n_clipped)
    } else if (!all(smd$balanced)) {
      sprintf("fair: %d covariate(s) with |SMD| >= 0.1 after matching", sum(!smd$balanced))
    } else {
      sprintf("SMD/retention thresholds passed: %.0f%% of exposed matched, all assessed |SMD| < 0.1", 100 * retained)
    }
    ob_before <- names(outcome_balance)
    rows[[k]] <- as.data.frame(data.table::rbindlist(lapply(idx, function(i) {
      er <- PS_EffectRow(DT, setup, i, out$eff, m, fit, crude, covs, balance_covariates, info,
                         retention = data.frame(n_exposed_unmatched = m$n_unmatched, exposed_retained = retained),
                         quality = data.frame(matching_quality = quality, stringsAsFactors = FALSE))
      outcome_balance[[paste(k, setup$ModelID[i], setup$DependentVariable[i], setup$Type[i], sep = " | ")]] <<- er$outcome_balance
      er$row
    }), fill = TRUE))
    Matching_CheckpointSave(cp, n_k, k, fp, rows[[k]], balance[[k]], outcome_balance[setdiff(names(outcome_balance), ob_before)])
  }
  res <- Matching_CompileResults(rows, setup, extra_numeric = c("n_exposed_unmatched", "exposed_retained", "ps_at_bounds"))
  res$ExposureVariable <- setup$ExposureVariable[match(res$ModelID,setup$ModelID)]
  res$IndependentVariable <- setup$.IndependentExpression[match(res$ModelID,setup$ModelID)]
  res$Ratio <- setup$Ratio[match(res$ModelID,setup$ModelID)]
  res$PSCaliperSD <- setup$PSCaliperSD[match(res$ModelID,setup$ModelID)]
  res$CovariateCalipers <- setup$CovariateCalipers[match(res$ModelID,setup$ModelID)]
  balance <- do.call(rbind, balance)
  outcome_balance <- data.table::rbindlist(outcome_balance, fill = TRUE)
  res <- Matching_SelectColumns(res, detail)   # per-covariate balance stays in the Balance sheet (and OutcomeBalance with detail = TRUE)
  if (!is.null(Fpath)) Matching_WriteResults(res, balance, outcome_balance, flags, Fpath, prefix = "PS", detail = detail)
  Matching_CheckpointClear(cp)   # stage 2 is written; the checkpoints are no longer needed
  message(sprintf("%d models: %d significant unadjusted, %d after FDR within dependent, %d after overall FDR.", nrow(res),
                  sum(res$sig_unadjusted, na.rm = TRUE), sum(res$sig_fdr, na.rm = TRUE), sum(res$sig_fdr_overall, na.rm = TRUE)))
  attr(res, "balance") <- balance
  attr(res, "outcome_balance") <- outcome_balance
  attr(res, "flags") <- flags
  res
}

#' Inverse probability of treatment weighted regression driven by a setup table
#'
#' The weighting counterpart of \code{PSRegression()}, using the same setup
#' table and result layout. For each exposure x matching specification x
#' method: a propensity score (\code{PSFit()}), inverse probability weights
#' (\pkg{WeightIt}) times the survey weights, then each row's outcome model on
#' every record with those weights. \code{"svyglm"} uses the survey-design
#' covariance; \code{"glm"} and \code{"lm"} use the HC1 sandwich covariance.
#' Both treat the weights as fixed. With \code{Fpath}, one
#' \code{IPTW_Sub_DataPredictions_<dependent>.xlsx} per dependent and
#' \code{IPTW_README.xlsx} are written.
#'
#' @inheritParams PSRegression
#' @param estimand \code{"ATE"} (whole population), \code{"ATT"} (exposed) or
#'   \code{"ATO"} (overlap population).
#' @param trim \code{NULL}, or a quantile (for example 0.99) at which the
#'   weights are capped.
#' @return A data.table with one row per setup row, with attributes
#'   \code{balance}, \code{outcome_balance} and \code{flags}.
#' @examples
#' \dontrun{
#' res <- IPTWRegression(model_dt, setup, weights = "DISCWT",
#'   cluster = "cluster_id", strata = "strata_id", estimand = "ATE")
#' }
#' @export
IPTWRegression <- function(DT, setup, weights = NULL, cluster = NULL, strata = NULL, method = "glm", type = NULL,
                           estimand = "ATE", trim = NULL, Fpath = NULL, balance_covariates = NULL, flags = NULL, checkpoint = TRUE, detail = FALSE) {
  if (length(estimand) != 1L || !estimand %in% c("ATE", "ATT", "ATO")) stop("estimand must be ATE, ATT or ATO.")
  if (!is.null(trim) && (!is.numeric(trim) || length(trim) != 1L || !is.finite(trim) || trim <= 0.5 || trim >= 1))
    stop("trim must be NULL or one quantile between 0.5 and 1, e.g. 0.99.")
  DT <- as.data.frame(DT)
  prep <- PS_PrepareSetup(DT, setup, weights, cluster, strata, method, type, balance_covariates, flags, verb = "weighted")
  setup <- prep$setup; balance_covariates <- unique(c(prep$balance_covariates, unlist(lapply(setup$.OutcomeAdjust, PS_BaseVars))))
  crude <- PS_CrudeEffects(DT, setup, weights, cluster, strata, estimand = estimand)
  key <- paste(setup$IndependentVariable, setup$MatchingVariable, setup$Method, setup$.OutcomeAdjust, sep = " | ")
  cp <- Matching_CheckpointDir(Fpath, "IPTW", checkpoint)
  settings <- list(analysis = "IPTW", weights = weights, cluster = cluster, strata = strata, balance_covariates = balance_covariates,
                   estimand = estimand, trim = trim)
  rows <- list(); balance <- list(); outcome_balance <- list()
  for (k in unique(key)) {
    idx <- which(key == k); ind <- setup$IndependentVariable[idx[1]]; covs <- setdiff(PS_SplitVars(setup$MatchingVariable[idx[1]]), ind)
    adjust <- if (nzchar(setup$.OutcomeAdjust[idx[1]])) PS_SplitVars(setup$.OutcomeAdjust[idx[1]]) else NULL
    meth <- setup$Method[idx[1]]; svy <- meth == "svyglm"; ps_method <- if (svy) "svyglm" else "glm"
    n_k <- match(k, unique(key))
    fp <- Matching_Fingerprint(DT, setup[idx, ], settings, c(setup$DependentVariable[idx], ind, PS_BaseVars(covs),
                                                            PS_BaseVars(setup$.OutcomeAdjust[idx]), weights, cluster, strata, balance_covariates))
    saved <- Matching_CheckpointLoad(cp, n_k, k, fp)
    if (!is.null(saved)) {
      message(sprintf("[%d / %d] %s (from checkpoint)", n_k, length(unique(key)), k))
      rows[[k]] <- saved$rows; balance[[k]] <- saved$balance; outcome_balance[names(saved$outcome_balance)] <- saved$outcome_balance
      next
    }
    message(sprintf("[%d / %d] %s", n_k, length(unique(key)), k))
    out <- tryCatch({
      fit <- PSFit(DT, ind, covs, method = ps_method, weights = weights, cluster = if (svy) cluster else NULL, strata = if (svy) strata else NULL)
      m <- PSMatchWeight(fit, approach = "weight", estimand = estimand, trim = trim,
                         balance_covariates = setdiff(unique(c(balance_covariates, PS_BaseVars(covs))), ind))
      eff <- lapply(setNames(nm = unique(setup$Type[idx])), function(tp)
        PSOutcomeCompile(list(weighted = m), unique(setup$DependentVariable[idx][setup$Type[idx] == tp]), family = PS_TypeFamily[[tp]], adjust = adjust, PrintSummary = FALSE))
      list(fit = fit, m = m, eff = eff)
    }, error = function(err) err)
    if (inherits(out, "error")) { rows[[k]] <- PS_FailedRows(setup, idx, conditionMessage(out)); next }
    m <- out$m; fit <- out$fit; dg <- m$diagnostics; ex <- dg$group == "exposed"; un <- dg$group == "unexposed"
    balance[[k]] <- data.frame(independent = ind, MatchingVariable = setup$MatchingVariable[idx[1]], Method = meth, m$balance, row.names = NULL)
    kept <- c(exposed = dg$ess_after[ex] / dg$ess_before[ex], controls = dg$ess_after[un] / dg$ess_before[un])
    info <- PS_PropensitySummary(DT, fit, m, ind, ps_method, quality_label = "weighting", weighting = list(
      estimand = estimand, trim = if (is.null(trim)) NA_real_ else trim,
      ess_kept_exposed = kept[["exposed"]], ess_kept_controls = kept[["controls"]]))
    smd <- m$balance
    quality <- if (min(kept) < 0.25 || fit$n_clipped > 0) {
      sprintf("poor: ESS kept %.0f%% (exposed) / %.0f%% (controls), %d scores at bounds; check the weighting variable(s) or trim",
              100 * kept[["exposed"]], 100 * kept[["controls"]], fit$n_clipped)
    } else if (!all(smd$balanced)) {
      sprintf("fair: %d covariate(s) with |SMD| >= 0.1 after weighting", sum(!smd$balanced))
    } else {
      sprintf("SMD/ESS thresholds passed: ESS kept %.0f%% (exposed) / %.0f%% (controls), all assessed |SMD| < 0.1",
              100 * kept[["exposed"]], 100 * kept[["controls"]])
    }
    ob_before <- names(outcome_balance)
    rows[[k]] <- as.data.frame(data.table::rbindlist(lapply(idx, function(i) {
      er <- PS_EffectRow(DT, setup, i, out$eff, m, fit, crude, covs, balance_covariates, info, retention = NULL,
                         quality = data.frame(weighting_quality = quality, stringsAsFactors = FALSE))
      outcome_balance[[paste(k, setup$ModelID[i], setup$DependentVariable[i], setup$Type[i], sep = " | ")]] <<- er$outcome_balance
      er$row
    }), fill = TRUE))
    Matching_CheckpointSave(cp, n_k, k, fp, rows[[k]], balance[[k]], outcome_balance[setdiff(names(outcome_balance), ob_before)])
  }
  res <- Matching_CompileResults(rows, setup, extra_numeric = "ps_at_bounds")
  res$ExposureVariable <- setup$ExposureVariable[match(res$ModelID, setup$ModelID)]
  res$IndependentVariable <- setup$.IndependentExpression[match(res$ModelID, setup$ModelID)]
  for (nm in intersect(c("Ratio", "PSCaliperSD", "CovariateCalipers"), names(setup)))
    res[[nm]] <- setup[[nm]][match(res$ModelID, setup$ModelID)]
  res$matching_settings_note <- "Ratio, PSCaliperSD and CovariateCalipers are matching-only settings; not applied to weighting."
  balance <- do.call(rbind, balance)
  outcome_balance <- data.table::rbindlist(outcome_balance, fill = TRUE)
  res <- Matching_SelectColumns(res, detail)   # per-covariate balance stays in the Balance sheet (and OutcomeBalance with detail = TRUE)
  if (!is.null(Fpath)) Matching_WriteResults(res, balance, outcome_balance, flags, Fpath, prefix = "IPTW", detail = detail)
  Matching_CheckpointClear(cp)   # stage 2 is written; the checkpoints are no longer needed
  message(sprintf("%d models: %d significant unadjusted, %d after FDR within dependent, %d after overall FDR.", nrow(res),
                  sum(res$sig_unadjusted, na.rm = TRUE), sum(res$sig_fdr, na.rm = TRUE), sum(res$sig_fdr_overall, na.rm = TRUE)))
  attr(res, "balance") <- balance
  attr(res, "outcome_balance") <- outcome_balance
  attr(res, "flags") <- flags
  res
}

#' Entropy balancing weighted regression driven by a setup table
#'
#' The entropy balancing counterpart of \code{PSRegression()}, using the same
#' setup table and result layout. For each exposure x matching specification x
#' method: entropy balancing weights (\code{EBWeight()}; no propensity score)
#' times the survey weights, then each row's outcome model on every record with
#' those weights. With \code{Fpath}, one
#' \code{EB_Sub_DataPredictions_<dependent>.xlsx} per dependent and
#' \code{EB_README.xlsx} are written.
#'
#' @inheritParams PSRegression
#' @param estimand \code{"ATT"} (controls reweighted to the exposed) or
#'   \code{"ATE"} (both groups reweighted to the whole population).
#' @param moments 1 balances means; 2 also balances variances.
#' @return A data.table with one row per setup row, with attributes
#'   \code{balance}, \code{outcome_balance} and \code{flags}.
#' @examples
#' \dontrun{
#' res <- EBRegression(model_dt, setup, weights = "DISCWT",
#'   cluster = "cluster_id", strata = "strata_id", estimand = "ATT")
#' }
#' @export
EBRegression <- function(DT, setup, weights = NULL, cluster = NULL, strata = NULL, method = "glm", type = NULL,
                         estimand = "ATT", moments = 1, Fpath = NULL, balance_covariates = NULL, flags = NULL, checkpoint = TRUE, detail = FALSE) {
  if (length(estimand) != 1L || !estimand %in% c("ATT", "ATE")) stop("estimand must be ATT or ATE.")
  if (length(moments) != 1L || !moments %in% c(1, 2)) stop("moments must be 1 (means) or 2 (means and variances).")
  DT <- as.data.frame(DT)
  prep <- PS_PrepareSetup(DT, setup, weights, cluster, strata, method, type, balance_covariates, flags, verb = "balanced")
  setup <- prep$setup; balance_covariates <- unique(c(prep$balance_covariates, unlist(lapply(setup$.OutcomeAdjust, PS_BaseVars))))
  crude <- PS_CrudeEffects(DT, setup, weights, cluster, strata, estimand = estimand)
  key <- paste(setup$IndependentVariable, setup$MatchingVariable, setup$Method, setup$.OutcomeAdjust, sep = " | ")
  cp <- Matching_CheckpointDir(Fpath, "EB", checkpoint)
  settings <- list(analysis = "EB", weights = weights, cluster = cluster, strata = strata, balance_covariates = balance_covariates,
                   estimand = estimand, moments = moments)
  rows <- list(); balance <- list(); outcome_balance <- list()
  for (k in unique(key)) {
    idx <- which(key == k); ind <- setup$IndependentVariable[idx[1]]; covs <- setdiff(PS_SplitVars(setup$MatchingVariable[idx[1]]), ind)
    adjust <- if (nzchar(setup$.OutcomeAdjust[idx[1]])) PS_SplitVars(setup$.OutcomeAdjust[idx[1]]) else NULL
    meth <- setup$Method[idx[1]]; svy <- meth == "svyglm"
    n_k <- match(k, unique(key))
    fp <- Matching_Fingerprint(DT, setup[idx, ], settings, c(setup$DependentVariable[idx], ind, PS_BaseVars(covs),
                                                            PS_BaseVars(setup$.OutcomeAdjust[idx]), weights, cluster, strata, balance_covariates))
    saved <- Matching_CheckpointLoad(cp, n_k, k, fp)
    if (!is.null(saved)) {
      message(sprintf("[%d / %d] %s (from checkpoint)", n_k, length(unique(key)), k))
      rows[[k]] <- saved$rows; balance[[k]] <- saved$balance; outcome_balance[names(saved$outcome_balance)] <- saved$outcome_balance
      next
    }
    message(sprintf("[%d / %d] %s", n_k, length(unique(key)), k))
    out <- tryCatch({
      m <- EBWeight(DT, ind, covs, weights = weights, cluster = if (svy) cluster else NULL, strata = if (svy) strata else NULL,
                    estimand = estimand, moments = moments, balance_covariates = setdiff(unique(c(balance_covariates, PS_BaseVars(covs))), ind))
      eff <- lapply(setNames(nm = unique(setup$Type[idx])), function(tp)
        PSOutcomeCompile(list(balanced = m), unique(setup$DependentVariable[idx][setup$Type[idx] == tp]), family = PS_TypeFamily[[tp]], adjust = adjust, PrintSummary = FALSE))
      list(m = m, eff = eff)
    }, error = function(err) err)
    if (inherits(out, "error")) { rows[[k]] <- PS_FailedRows(setup, idx, conditionMessage(out)); next }
    m <- out$m; dg <- m$diagnostics; ex <- dg$group == "exposed"; un <- dg$group == "unexposed"
    balance[[k]] <- data.frame(independent = ind, MatchingVariable = setup$MatchingVariable[idx[1]], Method = meth, m$balance, row.names = NULL)
    kept <- c(exposed = dg$ess_after[ex] / dg$ess_before[ex], controls = dg$ess_after[un] / dg$ess_before[un])
    info <- PS_PropensitySummary(DT, fit = NULL, m, ind, ps_method = NULL, quality_label = "balancing", weighting = list(
      weights_model = sprintf("ebal (moments = %d): %s", moments, deparse1(m$formula)), estimand = estimand, moments = moments,
      ess_kept_exposed = kept[["exposed"]], ess_kept_controls = kept[["controls"]]))
    smd <- m$balance
    target_smd <- abs(smd$smd_after[smd$in_propensity_model])
    worst_target <- if (any(is.finite(target_smd))) max(target_smd, na.rm = TRUE) else NA_real_
    quality <- if (min(kept) < 0.25 || !is.finite(worst_target) || worst_target >= 0.01) {
      sprintf("poor: ESS kept %.0f%% (exposed) / %.0f%% (controls), largest |SMD| on balanced covariates %.3f; check ps_fit_warnings and the balancing variable(s)",
              100 * kept[["exposed"]], 100 * kept[["controls"]], worst_target)
    } else if (!all(smd$balanced)) {
      sprintf("fair: %d other assessed covariate(s) with |SMD| >= 0.1 after balancing", sum(!smd$balanced))
    } else {
      sprintf("SMD/ESS thresholds passed: ESS kept %.0f%% (exposed) / %.0f%% (controls), all assessed |SMD| < 0.1",
              100 * kept[["exposed"]], 100 * kept[["controls"]])
    }
    ob_before <- names(outcome_balance)
    rows[[k]] <- as.data.frame(data.table::rbindlist(lapply(idx, function(i) {
      er <- PS_EffectRow(DT, setup, i, out$eff, m, fit = NULL, crude, covs, balance_covariates, info, retention = NULL,
                         quality = data.frame(balancing_quality = quality, stringsAsFactors = FALSE))
      outcome_balance[[paste(k, setup$ModelID[i], setup$DependentVariable[i], setup$Type[i], sep = " | ")]] <<- er$outcome_balance
      er$row
    }), fill = TRUE))
    Matching_CheckpointSave(cp, n_k, k, fp, rows[[k]], balance[[k]], outcome_balance[setdiff(names(outcome_balance), ob_before)])
  }
  res <- Matching_CompileResults(rows, setup)
  res$ExposureVariable <- setup$ExposureVariable[match(res$ModelID, setup$ModelID)]
  res$IndependentVariable <- setup$.IndependentExpression[match(res$ModelID, setup$ModelID)]
  for (nm in intersect(c("Ratio", "PSCaliperSD", "CovariateCalipers"), names(setup)))
    res[[nm]] <- setup[[nm]][match(res$ModelID, setup$ModelID)]
  res$matching_settings_note <- "Ratio, PSCaliperSD and CovariateCalipers are matching-only settings; not applied to weighting."
  balance <- do.call(rbind, balance)
  outcome_balance <- data.table::rbindlist(outcome_balance, fill = TRUE)
  res <- Matching_SelectColumns(res, detail)   # per-covariate balance stays in the Balance sheet (and OutcomeBalance with detail = TRUE)
  if (!is.null(Fpath)) Matching_WriteResults(res, balance, outcome_balance, flags, Fpath, prefix = "EB", detail = detail)
  Matching_CheckpointClear(cp)   # stage 2 is written; the checkpoints are no longer needed
  message(sprintf("%d models: %d significant unadjusted, %d after FDR within dependent, %d after overall FDR.", nrow(res),
                  sum(res$sig_unadjusted, na.rm = TRUE), sum(res$sig_fdr, na.rm = TRUE), sum(res$sig_fdr_overall, na.rm = TRUE)))
  attr(res, "balance") <- balance
  attr(res, "outcome_balance") <- outcome_balance
  attr(res, "flags") <- flags
  res
}
