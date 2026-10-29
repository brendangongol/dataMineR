#' Outcome family from the data
#'
#' @param y Outcome values.
#' @param family Optional family that overrides the data.
#' @return \code{"binomial"} for 0/1 values and \code{"gaussian"} otherwise.
#'   Count outcomes must request \code{"poisson"} explicitly because integer
#'   storage alone cannot distinguish counts from ordinal or rounded continuous
#'   outcomes.
#' @keywords internal
#' @noRd
PS_family <- function(y, family = NULL) {
  if (!is.null(family)) { return(match.arg(family, c("binomial", "poisson", "gaussian"))) }
  y <- y[!is.na(y)]
  if (length(y) && all(y %in% c(0, 1))) "binomial" else "gaussian"
}

#' Validate complete model inputs
#'
#' @param data Data frame.
#' @param columns Columns required by a model.
#' @param context Short model description used in the error.
#' @return \code{NULL}, invisibly.
#' @keywords internal
#' @noRd
PS_ValidateComplete <- function(data, columns, context) {
  columns <- unique(columns[!is.na(columns) & nzchar(columns)])
  missing_columns <- setdiff(columns, names(data))
  if (length(missing_columns))
    stop(context, " column(s) not found: ", paste(missing_columns, collapse = ", "), ".", call. = FALSE)
  bad <- vapply(columns, function(v) {
    x <- data[[v]]
    sum(is.na(x) | (is.numeric(x) & !is.finite(x)))
  }, integer(1))
  bad <- bad[bad > 0]
  if (length(bad))
    stop(context, " requires complete finite inputs. Problem columns: ",
         paste(sprintf("%s (%d)", names(bad), bad), collapse = ", "),
         ". Resolve or impute these values before fitting.", call. = FALSE)
  invisible(NULL)
}

#' Fit one weighted model
#'
#' Survey-weighted (\code{survey::svyglm}) when cluster or strata are given,
#' otherwise \code{glm} with optional weights rescaled to mean 1. Quasi families
#' allow non-integer weights. Stops when the fit does not converge or has
#' non-finite coefficients.
#'
#' @param formula Model formula.
#' @param data Data frame.
#' @param family \code{"binomial"}, \code{"poisson"} or \code{"gaussian"}.
#' @param weights,cluster,strata Design columns, or \code{NULL}.
#' @param lonely_psu Value of \code{options(survey.lonely.psu)} for the fit.
#' @return A list: \code{model}, \code{coefs} and \code{vcov}.
#' @keywords internal
#' @noRd
PS_glm_fit <- function(formula, data, family = "binomial", weights = NULL, cluster = NULL, strata = NULL, lonely_psu = "adjust") {
  fam <- switch(family, binomial = quasibinomial(), poisson = quasipoisson(), gaussian = gaussian())
  data <- as.data.frame(data)
  y_check <- model.response(model.frame(formula, data = data, na.action = na.fail))
  if (family == "binomial" && (!all(y_check %in% c(0, 1)) || length(unique(y_check)) != 2))
    stop("Binary model requires both 0 and 1 outcomes; no finite contrast is available otherwise.")
  if (!is.null(weights) && weights %in% names(data) &&
      any(!is.finite(data[[weights]]) | data[[weights]] <= 0)) stop("Analysis weights must be finite and positive.")
  missing_cols <- setdiff(c(weights, cluster, strata), names(data))
  if (length(missing_cols)) { stop("Survey design column(s) not found in the data: ", paste(missing_cols, collapse = ", "),
                                   ". Create them before this step (e.g. strata_id / cluster_id).", call. = FALSE) }
  old_options <- options(survey.lonely.psu = lonely_psu)
  on.exit(options(old_options))
  fit_warnings <- character()
  capture_warning <- function(w) { fit_warnings <<- c(fit_warnings, conditionMessage(w)) }
  if (!is.null(cluster) || !is.null(strata)) {
    if (is.null(weights)) { stop("A survey design needs weights.") }
    design <- survey::svydesign(ids = if (is.null(cluster)) ~1 else as.formula(paste0("~", cluster)),
                                strata = if (is.null(strata)) NULL else as.formula(paste0("~", strata)),
                                weights = as.formula(paste0("~", weights)), nest = TRUE, data = data)
    model <- withCallingHandlers(survey::svyglm(formula, design = design, family = fam, na.action = na.fail), warning = capture_warning)
  } else {
    # Weights rescaled to mean 1, as svyglm does: estimates and robust SEs are unchanged, but glm's IRLS can diverge
    # with large weights (e.g. inverse probability x discharge weights in the hundreds).
    data$.ps_fit_weight <- if (is.null(weights)) rep(1, nrow(data)) else data[[weights]] / mean(data[[weights]])
    model <- withCallingHandlers(glm(formula, data = data, family = fam, weights = .ps_fit_weight, na.action = na.fail), warning = capture_warning)
  }
  if (!isTRUE(model$converged) || any(!is.finite(model$coefficients)) || any(!is.finite(vcov(model))))
    stop("Model failed convergence or has non-finite/aliased coefficients or covariance; review predictors and overlap.")
  model$ps_warnings <- unique(fit_warnings)
  list(model = model, coefs = summary(model)$coefficients, vcov = vcov(model))
}

#' Positivity and overlap check before matching or weighting
#'
#' For each candidate covariate and exposure: its weighted correlation with the
#' exposure (phi for 0/1 covariates) and, for 0/1 covariates, the exposed share
#' among records with and without it. Flags are: \code{near_duplicate}
#' (\eqn{|r| \ge} \code{phi_threshold}: the covariate nearly is the exposure);
#' \code{joint_duplicate} (0/1 covariates that together nearly are the
#' exposure); \code{extreme_level} (a level whose exposed share leaves no
#' comparable patients for the estimand); and \code{post_exposure} (user
#' supplied). Flags annotate; nothing is removed. With \code{ps}, the share of
#' each group near the score bounds and the common support are also reported.
#'
#' @param data Data frame.
#' @param exposure One or more 0/1 exposure columns.
#' @param covariates Candidate covariates.
#' @param weights Name of the weight column, or \code{NULL}.
#' @param post_exposure Covariates known to be measured after the exposure.
#' @param estimand \code{"ATT"} and/or \code{"ATE"}; with several, the checks
#'   are run for each and stacked.
#' @param ps Optional fitted propensity scores for the overlap summary.
#' @param phi_threshold Absolute correlation at which a covariate is a
#'   near-duplicate.
#' @param joint_min_cor Minimum correlation for a 0/1 covariate to join a
#'   joint near-duplicate set.
#' @param share_threshold Exposed share defining an extreme level.
#' @param min_level_n Minimum records in a level before it can be extreme.
#' @param ps_bounds Score bounds for the overlap summary.
#' @param PrintSummary Logical; print the flagged covariates.
#' @param IPTWSampleSizeMinimum Minimum effective sample size per group, as a
#'   share of its records, for the weighting check.
#' @return Invisibly, a list: \code{table} (one row per covariate, exposure and
#'   estimand with the flags and suggested action), \code{keep} (all
#'   candidates), \code{checks}, \code{ps_summary} and \code{common_support}.
#' @examples
#' \dontrun{
#' chk <- MatchingAssumptionCheck(model_dt, exposure = "sud_any_01",
#'   covariates = c("AGE", "sex_01", "elix_drug", "elix_alcohol"),
#'   weights = "DISCWT", estimand = c("ATT", "ATE"))
#' chk$table
#' }
#' @export
MatchingAssumptionCheck <- function(data, exposure, covariates, weights = NULL, post_exposure = NULL, estimand = "ATT", ps = NULL,
                              phi_threshold = 0.8, joint_min_cor = 0.3, share_threshold = 0.01, min_level_n = 20,
                              ps_bounds = c(0.01, 0.99), PrintSummary = TRUE, IPTWSampleSizeMinimum = 0.10) {
  # Minimum group ESS as a fraction of its original record count (0.10 = 10%).
  if (!is.numeric(IPTWSampleSizeMinimum) || length(IPTWSampleSizeMinimum) != 1L ||
      !is.finite(IPTWSampleSizeMinimum) || IPTWSampleSizeMinimum < 0 || IPTWSampleSizeMinimum > 1)
    stop("IPTWSampleSizeMinimum must be a number between 0 and 1 (0.10 means 10%).")
  estimand <- unique(match.arg(estimand, c("ATT", "ATE"), several.ok = TRUE))
  if (length(estimand) > 1L) {
    runs <- setNames(lapply(estimand, function(es) MatchingAssumptionCheck(data, exposure, covariates, weights, post_exposure, es, ps,
      phi_threshold, joint_min_cor, share_threshold, min_level_n, ps_bounds, PrintSummary, IPTWSampleSizeMinimum)), estimand)
    stack <- function(nm) { x <- do.call(rbind, lapply(runs, `[[`, nm)); rownames(x) <- NULL; x }
    return(invisible(list(table = stack("table"), keep = unique(unlist(lapply(runs, `[[`, "keep"))), checks = stack("checks"),
      by_estimand = runs, ps_summary = lapply(runs, `[[`, "ps_summary"), common_support = lapply(runs, `[[`, "common_support"))))
  }
  data <- as.data.frame(data)
  absent <- setdiff(c(exposure, covariates, weights), names(data))
  if (length(absent)) stop("Columns not found: ", paste(absent, collapse=", "))
  if (!length(exposure) || anyNA(exposure) || any(!exposure %in% names(data)))
    stop("exposure must name one or more columns in data.")
  exposure <- unique(exposure)
  covariates <- unique(covariates)
  if (length(exposure) > 1L) {
    if (!is.null(ps) && (!is.list(ps) || !all(exposure %in% names(ps))))
      stop("For multiple exposures, ps must be NULL or a named list of score vectors.")
    checks <- setNames(lapply(exposure, function(a) MatchingAssumptionCheck(data, a, setdiff(covariates, a),
      weights, post_exposure, estimand, if (is.null(ps)) NULL else ps[[a]],
      phi_threshold, joint_min_cor, share_threshold, min_level_n, ps_bounds, PrintSummary, IPTWSampleSizeMinimum)), exposure)
    return(invisible(list(table = do.call(rbind, lapply(checks, `[[`, "table")),
      keep = unique(unlist(lapply(checks, `[[`, "keep"))), by_exposure = checks,
      checks = do.call(rbind, lapply(checks, `[[`, "checks")),
      ps_summary = lapply(checks, `[[`, "ps_summary"), common_support = lapply(checks, `[[`, "common_support"))))
  }
  covariates <- setdiff(covariates, exposure)
  if (!length(covariates)) stop("Supply at least one covariate to assess.")
  a <- data[[exposure]]
  if (!all(a %in% c(0, 1))) { stop("The exposure must be coded 0/1.") }
  w <- if (is.null(weights)) rep(1, nrow(data)) else data[[weights]]
  if (length(unique(a)) != 2L) stop("Both exposure groups must be present.")
  if (!is.numeric(w) || any(!is.finite(w) | w <= 0)) stop("Sampling weights must be finite and positive.")
  if (!is.null(ps) && (!is.numeric(ps) || length(ps) != nrow(data) || any(!is.finite(ps) | ps < 0 | ps > 1)))
    stop("ps must contain one finite probability per data row.")
  table <- do.call(rbind, lapply(covariates, function(v) {
    x <- data[[v]]; binary <- all(na.omit(x) %in% c(0, 1)) && any(!is.na(x))
    ok <- !is.na(x)
    if (is.numeric(x)) ok <- ok & is.finite(x)
    r <- if (is.numeric(x) && sum(ok)>1 && length(unique(x[ok]))>1 && length(unique(a[ok]))>1)
      stats::cov.wt(cbind(x[ok], a[ok]), wt=w[ok], cor=TRUE)$cor[1,2] else NA_real_
    x[!ok] <- NA
    share1 <- if (binary && any(x == 1, na.rm=TRUE)) weighted.mean(a[which(x == 1)], w[which(x == 1)]) else NA_real_
    share0 <- if (binary && any(x == 0, na.rm=TRUE)) weighted.mean(a[which(x == 0)], w[which(x == 0)]) else NA_real_
    n1 <- if (binary) sum(x == 1, na.rm=TRUE) else NA; n0 <- if (binary) sum(x == 0, na.rm=TRUE) else NA
    high <- function(s, n) !is.na(s) && n >= min_level_n && s >= 1 - share_threshold
    low <- function(s, n) !is.na(s) && n >= min_level_n && s <= share_threshold && estimand != "ATT"
    extreme <- high(share1, n1) || high(share0, n0) || low(share1, n1) || low(share0, n0)
    near_dup <- !is.na(r) && abs(r) >= phi_threshold
    post <- v %in% post_exposure
    action <- if (post) "drop: post-exposure (user-specified)" else if (near_dup) sprintf("drop: near-duplicate of the exposure (r = %.2f)", r) else
      if (extreme) "drop or combine levels: almost all one exposure group in a level" else "keep"
    data.frame(variable = v, missing_n = sum(!ok), constant = length(unique(x[ok])) < 2,
               binary = binary, correlation = r, n_with = n1, exposed_share_with = share1, exposed_share_without = share0,
               near_duplicate = near_dup, extreme_level = extreme, post_exposure = post, action = action, row.names = NULL)
  }))
  table$joint_duplicate <- FALSE
  pool <- table[table$binary & is.finite(table$correlation) & table$correlation >= joint_min_cor & table$missing_n == 0 & !table$post_exposure, ]
  pool <- pool[order(-pool$correlation), ]
  if (nrow(pool) >= 2) {
    phi_any <- function(vars) stats::cov.wt(cbind(as.numeric(rowSums(data[, vars, drop = FALSE]) > 0), a), wt = w, cor = TRUE)$cor[1, 2]
    chosen <- pool$variable[1]; best <- phi_any(chosen)
    for (v in pool$variable[-1]) { p <- phi_any(c(chosen, v)); if (is.finite(p) && is.finite(best) && p > best) { chosen <- c(chosen, v); best <- p } }
    if (length(chosen) >= 2 && best >= phi_threshold) {
      hit <- table$variable %in% chosen
      table$joint_duplicate[hit] <- TRUE
      table$action[hit & table$action == "keep"] <- sprintf("drop: together with %s nearly the exposure (r of any = %.2f)",
                                                            paste(chosen, collapse = ", "), best)
      table$action[hit] <- ifelse(grepl("^drop: near-duplicate", table$action[hit]),
                                  sprintf("drop: near-duplicate of the exposure, alone and with %s (r of any = %.2f)", paste(chosen, collapse = ", "), best),
                                  table$action[hit])
    }
  }
  ps_summary <- NULL; common_support <- NULL
  if (!is.null(ps)) {
    ps_summary <- do.call(rbind, lapply(c(1, 0), function(g) {
      p <- ps[a == g]
      data.frame(group = if (g == 1) "exposed" else "unexposed", records = length(p), ps_min = min(p), ps_max = max(p),
                 below_lower = mean(p < ps_bounds[1]), above_upper = mean(p > ps_bounds[2]))
    }))
    support <- c(lower = max(min(ps[a == 1]), min(ps[a == 0])), upper = min(max(ps[a == 1]), max(ps[a == 0])))
    outside <- ps < support[["lower"]] | ps > support[["upper"]]
    common_support <- list(range = support, outside_exposed = sum(outside & a == 1), outside_unexposed = sum(outside & a == 0))
  }
  # Necessary marginal support checks apply to all three approaches.
  table$unsupported_levels <- ""
  table$outside_control_range <- NA_integer_
  table$eb_mean_outside_support <- FALSE
  for (i in seq_len(nrow(table))) {
    x <- data[[table$variable[i]]]
    ok <- !is.na(x); if (is.numeric(x)) ok <- ok & is.finite(x)
    if (!is.numeric(x) || table$binary[i]) {
      lev1 <- unique(as.character(x[ok & a==1])); lev0 <- unique(as.character(x[ok & a==0]))
      bad <- if (estimand=="ATT") setdiff(lev1,lev0) else union(setdiff(lev1,lev0),setdiff(lev0,lev1))
      table$unsupported_levels[i] <- paste(bad,collapse="; ")
    } else if (any(ok & a==0) && any(ok & a==1)) {
      r0 <- range(x[ok & a==0]); r1 <- range(x[ok & a==1])
      table$outside_control_range[i] <- sum(ok & a==1 & (x<r0[1] | x>r0[2]),na.rm=TRUE)
      target <- weighted.mean(x[ok & (a==1 | estimand=="ATE")],w[ok & (a==1 | estimand=="ATE")])
      table$eb_mean_outside_support[i] <- target<r0[1] || target>r0[2] ||
        (estimand=="ATE" && (target<r1[1] || target>r1[2]))
    }
    issues <- c(if(table$missing_n[i]>0) "missing/nonfinite values",
      if(table$constant[i]) "constant or entirely missing covariate",
      if(nzchar(table$unsupported_levels[i])) "levels lack comparison-group support",
      if(isTRUE(table$outside_control_range[i]>0)) "exposed observations outside control range",
      if(table$eb_mean_outside_support[i]) "entropy target mean outside available range")
    if(length(issues)) table$action[i] <- paste(c(if(table$action[i]!="keep") table$action[i],paste("review:",issues)),collapse="; ")
  }
  checks <- data.frame(exposure=character(),algorithm=character(),check=character(),status=character(),detail=character())
  add <- function(algorithm,check,status,detail) {
    checks <<- rbind(checks,data.frame(exposure=exposure,algorithm=algorithm,check=check,status=status,detail=detail))
  }
  add("All", "Input groups and sampling weights", "checked", "Both groups present; sampling weights finite and positive.")
  add("All", "Measured covariate support", if(any(table$action!="keep")) "review" else "no marginal flags",
      "See covariate table. Marginal overlap does not establish joint overlap.")
  add("All", "Causal identification", "not testable from data",
      "Requires clinically justified baseline confounders, no unmeasured confounding, consistency and no interference; timing flags require post_exposure input.")
  add("PSM", "Matched balance and retention", "requires fitted results", "Inspect all covariate balance, unmatched counts and retained target population in PSResults.")
  add("IPTW", "Weighted balance and effective sample size", "requires fitted results", "Inspect actual final weights, group ESS and all covariate balance in IPTW results.")
  add("EB", "Moment feasibility", if(any(table$eb_mean_outside_support | nzchar(table$unsupported_levels))) "review" else "necessary marginal checks only",
      "Target means must be attainable; marginal ranges do not prove joint feasibility or feasibility of higher moments.")
  add("EB", "Optimization and achieved balance", "requires fitted results", "Inspect solver warnings, achieved moment balance, extreme weights and group ESS in EB results.")
  if (!is.null(ps)) {
    add("PSM/IPTW", "Propensity-score support", if(any(outside)) "review" else "no range flags",
      sprintf("Outside common score range: %d exposed; %d controls.",sum(outside & a==1),sum(outside & a==0)))
    iw <- if(estimand=="ATT") ifelse(a==1,1,ps/(1-ps)) else ifelse(a==1,1/ps,1/(1-ps))
    for(g in c(1,0)) {
      z <- (iw*w)[a==g]
      valid <- all(is.finite(z)) && sum(z)>0
      ess <- if(valid) sum(z)^2/sum(z^2) else NA_real_
      add("IPTW", paste("Untrimmed implied weights, group",g),if(!valid || ess<sum(a==g)*IPTWSampleSizeMinimum) "review" else "diagnostic",
          sprintf("ESS=%s of %d; maximum combined weight=%s. Uses supplied scores and sampling weights; no clipping or trimming.",signif(ess,4),sum(a==g),signif(max(z),4)))
    }
  } else add("PSM/IPTW", "Propensity-score support", "not assessed", "Supply ps from the fitted propensity model to assess score support and implied IPTW weights.")
  # Flags describe concerns; no candidate is automatically removed.
  table$action <- sub("^drop[^:]*: ", "review: ", table$action)
  table$flagged <- table$action != "keep"
  table$exposure <- exposure
  table$estimand <- estimand; checks$estimand <- estimand
  keep <- table$variable
  if (PrintSummary) {
    print(checks, row.names=FALSE)
    flagged <- table[table$action != "keep", c("variable", "correlation", "exposed_share_with", "action")]
    if (nrow(flagged)) { print(flagged, row.names = FALSE, digits = 3) } else { message("No covariates flagged.") }
    if (!is.null(ps_summary)) {
      print(ps_summary, row.names = FALSE, digits = 3)
      message(sprintf("Common support %.3g to %.3g: %d exposed and %d unexposed records outside it.", support[["lower"]], support[["upper"]],
                      common_support$outside_exposed, common_support$outside_unexposed))
    }
  }
  invisible(list(table = table, keep = keep, checks = checks, ps_summary = ps_summary, common_support = common_support))
}

#' Fit a propensity score model
#'
#' Scores are bounded to \code{[1e-6, 1 - 1e-6]} so their logit, used for the
#' matching caliper, is finite.
#'
#' @param data Data frame.
#' @param exposure The 0/1 exposure column.
#' @param covariates Covariates (formula terms such as \code{"AGE + I(AGE^2)"}
#'   are allowed).
#' @param method \code{"svyglm"} (survey-weighted logistic; needs weights and
#'   cluster or strata), \code{"glm"} (logistic, optional weights) or
#'   \code{"lasso"} (\code{glmnet::cv.glmnet} at lambda.min; ordinary
#'   cross-validation).
#' @param weights,cluster,strata Design columns, or \code{NULL}.
#' @param lonely_psu Value of \code{options(survey.lonely.psu)} for the fit.
#' @param nfolds Folds for \code{method = "lasso"}.
#' @param seed Random seed for \code{method = "lasso"}.
#' @return A \code{"PSFit"} list: \code{ps}, \code{logit_ps}, \code{model},
#'   \code{method}, \code{exposure}, \code{covariates}, \code{data} (with
#'   \code{.ps}), the design columns, \code{formula}, \code{overlap_plot} and
#'   \code{n_clipped} (scores at the bounds).
#' @examples
#' \dontrun{
#' fit <- PSFit(model_dt, "sud_any_01", c("AGE", "sex_01"), method = "glm",
#'   weights = "DISCWT")
#' }
#' @export
PSFit <- function(data, exposure, covariates, method = c("svyglm", "glm", "lasso"), weights = NULL, cluster = NULL, strata = NULL,
                  lonely_psu = "adjust", nfolds = 5, seed = 1621) {
  method <- match.arg(method)
  data <- as.data.frame(data)
  f <- if (length(covariates)) reformulate(covariates, exposure, env=environment(PS_Formula(covariates))) else as.formula(paste(exposure, "~ 1"))
  PS_ValidateComplete(data, c(all.vars(f), weights, cluster, strata), "Propensity model")
  if (!all(data[[exposure]] %in% c(0, 1)) || length(unique(data[[exposure]])) != 2L)
    stop("Propensity exposure must contain both 0 and 1 and no other values.", call. = FALSE)
  if (!is.null(weights) && any(!is.finite(data[[weights]]) | data[[weights]] <= 0))
    stop("Propensity-model weights must be finite and positive.", call. = FALSE)
  if (method == "lasso") {
    if (!is.null(seed)) { set.seed(seed) }
    x <- model.matrix(f, data)[, -1, drop=FALSE]
    model <- glmnet::cv.glmnet(x, data[[exposure]], family = "binomial", nfolds = nfolds,
                               weights = if (is.null(weights)) NULL else data[[weights]])
    ps <- as.numeric(predict(model, newx = x, s = "lambda.min", type = "response"))
  } else {
    if (method == "svyglm" && is.null(cluster) && is.null(strata)) { stop("method = 'svyglm' needs cluster and/or strata (or use 'glm').") }
    fit <- PS_glm_fit(f, data, "binomial", weights, if (method == "svyglm") cluster else NULL, if (method == "svyglm") strata else NULL, lonely_psu)
    model <- fit$model
    ps <- as.numeric(fitted(model))
  }
  if (length(ps) != nrow(data) || any(!is.finite(ps))) stop("Invalid propensity predictions.")
  n_clipped <- sum(ps < 1e-6 | ps > 1 - 1e-6)
  if (n_clipped) warning(n_clipped, " propensity scores at numerical bounds; inspect separation and overlap.")
  ps <- pmin(pmax(ps, 1e-6), 1 - 1e-6)
  data$.ps <- ps
  plot_df <- data.frame(logit_ps = qlogis(ps), group = factor(ifelse(data[[exposure]] == 1, "exposed", "unexposed")),
                        w = if (is.null(weights)) 1 else data[[weights]])
  overlap_plot <- ggplot2::ggplot(plot_df, ggplot2::aes(x = logit_ps, fill = group, weight = w)) +
    ggplot2::geom_density(alpha = 0.4) + ggplot2::labs(x = "logit(propensity score)", y = "Density", fill = NULL) + ggplot2::theme_minimal()
  structure(list(ps = ps, logit_ps = qlogis(ps), model = model, method = method, exposure = exposure, covariates = PS_BaseVars(covariates),
                 data = data, weights = weights, cluster = cluster, strata = strata, formula = f, overlap_plot = overlap_plot, n_clipped = n_clipped),
            class = "PSFit")
}

#' Covariate balance and group diagnostics
#'
#' Standardized mean differences, variance ratios and Kolmogorov-Smirnov
#' statistics before and after matching or weighting (\pkg{cobalt}), and per
#' group: records, effective sample sizes, weight range and, when a score
#' exists, score range and records outside the common support.
#'
#' @param d Data with the exposure and covariates (and \code{.ps} when
#'   available).
#' @param a Exposure column.
#' @param sw Survey weights.
#' @param after Analysis weights after matching or weighting.
#' @param balance_covariates Covariates whose balance is assessed.
#' @param model_covariates Covariates in the propensity or balancing model.
#' @param estimand Target estimand.
#' @param balance_threshold Absolute SMD below which a covariate is balanced.
#' @param exact Optional exact-matching covariates.
#' @param n_clipped Optional number of scores at the bounds.
#' @param fit_warnings Optional warnings to report.
#' @return A list: \code{bal_tab}, \code{balance} and \code{diagnostics}.
#' @keywords internal
#' @noRd
PS_BalanceDiagnostics <- function(d, a, sw, after, balance_covariates, model_covariates, estimand, balance_threshold = 0.1,
                                  exact = NULL, n_clipped = NULL, fit_warnings = NULL) {
  has_ps <- ".ps" %in% names(d)
  bt <- cobalt::bal.tab(x = d[, balance_covariates, drop = FALSE], treat = d[[a]],
                       weights = after / sw, s.weights = sw, estimand = estimand,
                       stats = c("m", "v", "ks"), binary = "std", un = TRUE, quick = FALSE,
                       thresholds = c(m = balance_threshold))
  b <- bt$Balance
  col <- function(nm) if (nm %in% names(b)) b[[nm]] else rep(NA_real_, nrow(b))
  balance <- data.frame(variable = rownames(b), mean_control_before = col("M.0.Un"),
      mean_exposed_before = col("M.1.Un"), mean_control_after = col("M.0.Adj"),
      mean_exposed_after = col("M.1.Adj"), smd_before = col("Diff.Un"), smd_after = col("Diff.Adj"),
      var_ratio_before = col("V.Ratio.Un"), var_ratio_after = col("V.Ratio.Adj"),
      ks_before = col("KS.Un"), ks_after = col("KS.Adj"), row.names = NULL)
  balance$balanced <- is.finite(balance$smd_after) & abs(balance$smd_after) < balance_threshold
  # An SMD this large means the exposed group has almost no variance in the covariate, i.e. the covariate is a
  # near-duplicate of the exposure. Report it here, named, instead of through cobalt's anonymous console warning.
  extreme <- balance$variable[(is.finite(balance$smd_before) & abs(balance$smd_before) > 5) |
                              (is.finite(balance$smd_after) & abs(balance$smd_after) > 5)]
  if (length(extreme)) fit_warnings <- c(fit_warnings, sprintf(
    "|SMD| above 5 for %s: the exposed group has almost no variance in it (usually a near-duplicate of the exposure); consider leaving it out of the balance covariates.",
    paste(extreme, collapse = ", ")))
  # Ask cobalt for its actual expanded row names; prefix matching can confuse
  # similarly named columns and does not reliably handle factor-level labels.
  row_map <- setNames(lapply(balance_covariates, function(v) {
    if (is.numeric(d[[v]])) return(v)
    one <- cobalt::bal.tab(x=d[, v, drop=FALSE], treat=d[[a]],
      weights=after/sw, s.weights=sw, estimand=estimand,
      stats="m", binary="std", un=FALSE, quick=TRUE)
    rownames(one$Balance)
  }), balance_covariates)
  rows_for <- function(vars) unique(unlist(row_map[intersect(vars, names(row_map))], use.names=FALSE))
  balance$in_propensity_model <- balance$variable %in% rows_for(model_covariates)
  balance$exact_match <- balance$variable %in% rows_for(exact)
  ess_fun <- function(w) if (sum(w^2) > 0) sum(w)^2 / sum(w^2) else 0
  diagnostics <- do.call(rbind, lapply(c(1, 0), function(g) {
    id <- d[[a]] == g; retained <- id & after > 0
    out <- data.frame(group = if (g == 1) "exposed" else "unexposed", n_before = sum(id),
       n_after = sum(retained), n_unmatched = sum(id & !retained), retained_fraction = sum(retained)/sum(id),
       weighted_retained_fraction = sum(sw[retained])/sum(sw[id]),
       ess_before = ess_fun(sw[id]), ess_after = ess_fun(after[id]),
       weight_min = if (any(retained)) min(after[retained]) else NA,
       weight_max = if (any(retained)) max(after[retained]) else NA)
    if (has_ps) {
      out$ps_min <- min(d$.ps[id]); out$ps_max <- max(d$.ps[id])
      out$fraction_ps_extreme <- mean(d$.ps[id] < .01 | d$.ps[id] > .99)
    }
    out
  }))
  if (has_ps) {
    support <- c(max(tapply(d$.ps, d[[a]], min)), min(tapply(d$.ps, d[[a]], max)))
    diagnostics$outside_common_support <- vapply(c(1, 0), function(g)
      sum(d[[a]] == g & (d$.ps < support[1] | d$.ps > support[2])), integer(1))
  }
  diagnostics$max_abs_smd_after <- if (any(is.finite(balance$smd_after))) max(abs(balance$smd_after), na.rm = TRUE) else NA_real_
  diagnostics$n_balance_flags <- sum(!balance$balanced)
  if (!is.null(n_clipped)) diagnostics$n_clipped_scores <- n_clipped
  diagnostics$fit_warnings <- paste(fit_warnings, collapse = "; ")
  diagnostics$review_status <- if (all(balance$balanced)) "SMD threshold passed; review distributions and overlap" else "Review imbalance before interpreting outcomes"
  # Include group-level diagnostics in the main returned balance table as well.
  for (i in seq_len(nrow(diagnostics))) for (nm in setdiff(names(diagnostics), "group"))
    balance[[paste0(diagnostics$group[i], "_", nm)]] <- diagnostics[[nm]][i]
  list(bal_tab = bt, balance = balance, diagnostics = diagnostics)
}

#' Love plot of absolute standardized mean differences
#'
#' Wraps \code{cobalt::love.plot()} and drops its generic "Large mean
#' differences detected" warning, which names neither the model nor the
#' covariate; the same condition is reported per covariate in
#' \code{PS_BalanceDiagnostics()}'s \code{fit_warnings}.
#'
#' @param bt A \code{cobalt::bal.tab()} result.
#' @param balance_threshold Absolute SMD below which a covariate is balanced.
#' @return A ggplot.
#' @keywords internal
#' @noRd
PS_LovePlot <- function(bt, balance_threshold) {
  withCallingHandlers(
    cobalt::love.plot(bt, stats = "m", abs = TRUE, thresholds = c(m = balance_threshold)),
    warning = function(w) if (grepl("Large mean differences detected", conditionMessage(w), fixed = TRUE)) invokeRestart("muffleWarning"))
}

#' Match or weight on a propensity score and check balance
#'
#' \code{"match"}: \pkg{MatchIt} nearest-neighbour matching on the logit score,
#' \code{ratio} controls per exposed record, \code{caliper} in standard
#' deviations of the logit score, optional exact matching and covariate
#' calipers; estimates the ATT. If MatchIt finds no match with covariate
#' calipers, the match is retried on the negated score. \code{"weight"}:
#' \pkg{WeightIt} weights from the fitted score for the ATT, ATE or ATO, with
#' optional trimming. \code{"none"}: survey weights only. The analysis weight
#' \code{.ps_weight} is the matching or propensity weight times the survey
#' weight.
#'
#' @param psfit A \code{PSFit()} result.
#' @param approach \code{"match"}, \code{"weight"} or \code{"none"}.
#' @param estimand \code{"ATT"}, \code{"ATE"} or \code{"ATO"} (matching is ATT).
#' @param ratio Controls per exposed record.
#' @param caliper Caliper in standard deviations of the logit score.
#' @param exact Optional exact-matching covariates (for example a year).
#' @param replace Logical; match with replacement.
#' @param trim Optional weight quantile to cap weights at.
#' @param balance_threshold Absolute SMD below which a covariate is balanced.
#' @param balance_covariates Covariates whose balance is assessed.
#' @param covariate_calipers Optional named vector of maximum differences per
#'   covariate.
#' @return A \code{"PSMatchWeight"} list: \code{approach}, \code{estimand},
#'   \code{data} (analysis records with \code{.ps_weight}), \code{balance},
#'   \code{diagnostics}, \code{max_smd_after}, \code{love_plot}, \code{ess},
#'   \code{n_unmatched}, \code{object} and the design fields.
#' @examples
#' \dontrun{
#' m <- PSMatchWeight(fit, approach = "match", ratio = 3, caliper = 0.2)
#' m$balance
#' }
#' @export
PSMatchWeight <- function(psfit, approach = c("match", "weight", "none"), estimand = c("ATT", "ATE", "ATO"), ratio = 3, caliper = 0.2,
                          exact = NULL, replace = FALSE, trim = NULL, balance_threshold = 0.1, balance_covariates = psfit$covariates, covariate_calipers = NULL) {
  approach <- match.arg(approach); estimand <- match.arg(estimand)
  if (!inherits(psfit, "PSFit")) { stop("psfit must come from PSFit().") }
  d <- psfit$data; a <- psfit$exposure; f <- psfit$formula
  d$.ps_row_id <- seq_len(nrow(d))
  balance_covariates <- unique(c(balance_covariates, exact, names(covariate_calipers)))
  if (replace && approach == "match") stop("This workflow supports matching without replacement; replacement needs repeated-unit inference.")
  sw <- if (is.null(psfit$weights)) rep(1, nrow(d)) else d[[psfit$weights]]
  object <- NULL; n_unmatched <- NA; match_note <- NULL
  if (approach == "match") {
    if (estimand != "ATT") { stop("Nearest-neighbour matching estimates ATT; use approach = 'weight' for ATE or ATO.") }
    widths <- caliper; standardized <- TRUE
    if (length(covariate_calipers)) {
      if (is.null(names(covariate_calipers)) || any(!names(covariate_calipers) %in% names(d)) || any(!is.finite(covariate_calipers) | covariate_calipers <= 0))
        stop("Covariate calipers must be positive named widths for existing variables.")
      widths <- c(setNames(caliper, ""), covariate_calipers)
      standardized <- setNames(c(TRUE, rep(FALSE,length(covariate_calipers))), names(widths))
    }
    run_match <- function(distance) MatchIt::matchit(f, data = d, distance = distance, method = "nearest", estimand = "ATT", ratio = ratio,
                               caliper = widths, std.caliper = standardized, replace = replace,
                               exact = if (is.null(exact)) NULL else reformulate(exact), s.weights = sw)
    # MatchIt 4.7.2 returns "No units were matched" with a covariate caliper when the supplied score decreases linearly
    # in that covariate (e.g. a score from AGE alone with a negative coefficient). The negated score has the same
    # distances and SD-unit caliper, so on that error the match is retried with it (ties may then break differently).
    object <- tryCatch(run_match(psfit$logit_ps), error = function(e) {
      if (!length(covariate_calipers) || !grepl("No units were matched", conditionMessage(e), fixed = TRUE)) stop(e)
      match_note <<- "Matched on the negated logit score: MatchIt found no matches with the covariate caliper on the original score."
      run_match(-psfit$logit_ps)
    })
    out <- MatchIt::match.data(object, data = d, include.s.weights = TRUE)
    out$.ps_weight <- out$weights  # already matching weight times survey weight
    n_unmatched <- sum(d[[a]] == 1) - sum(out[[a]] == 1)
  } else if (approach == "weight") {
    object <- WeightIt::weightit(f, data = d, ps = psfit$ps, estimand = estimand, s.weights = sw)
    if (!is.null(trim)) { object <- WeightIt::trim(object, at = trim) }
    out <- d
    out$.ps_weight <- object$weights * sw
  } else {
    out <- d
    out$.ps_weight <- sw
  }
  # Use the same combined weights for diagnostics and outcome analysis.
  after <- numeric(nrow(d)); after[out$.ps_row_id] <- out$.ps_weight
  bd <- PS_BalanceDiagnostics(d, a, sw, after, balance_covariates, psfit$covariates, estimand, balance_threshold, exact = exact,
                              n_clipped = psfit$n_clipped, fit_warnings = c(psfit$model$ps_warnings, match_note))
  bt <- bd$bal_tab; balance <- bd$balance; diagnostics <- bd$diagnostics
  ess <- setNames(diagnostics$ess_after, diagnostics$group)
  target_population <- if (approach == "match" && n_unmatched > 0)
    "Retained matched exposed population (not all exposed patients)" else if (estimand == "ATT") "Exposed population" else estimand
  structure(list(approach = approach, estimand = estimand, target_population = target_population,
      data = out, balance = balance, diagnostics = diagnostics,
      max_smd_after = diagnostics$max_abs_smd_after[1],
      love_plot = PS_LovePlot(bt, balance_threshold),
      ess = ess, n_unmatched = n_unmatched, object = object, exposure = a, covariates = psfit$covariates,
      weights = psfit$weights, cluster = psfit$cluster, strata = psfit$strata), class = "PSMatchWeight")
}

#' Exposure effects on many dependents in matched or weighted designs
#'
#' For each design and dependent, fits \code{dependent ~ exposure (+ adjust)}
#' with the analysis weight (survey-weighted when the design has cluster or
#' strata) and reports marginal effects by g-computation over the estimand's
#' target population: risk difference, risk ratio and odds ratio (binomial);
#' mean difference and ratio of expected counts (poisson); mean difference
#' (gaussian), with delta-method intervals. The intervals treat the matching or
#' weights as fixed. Bonferroni and FDR p-values are within each design and
#' measure, across dependents.
#'
#' @param designs Named list of \code{PSMatchWeight()} results.
#' @param dependents Outcome columns.
#' @param family \code{NULL} (0/1 outcomes use binomial; all others use
#'   gaussian), one family, or a named vector by dependent. Request poisson
#'   explicitly for count outcomes.
#' @param adjust Outcome-model covariates: a character vector for every design,
#'   or a named list by design.
#' @param level Confidence level.
#' @param Fpath Optional folder; one \code{PS_Effects_<dependent>.xlsx} per
#'   dependent is written.
#' @param lonely_psu Value of \code{options(survey.lonely.psu)} for the fits.
#' @param PrintSummary Logical; print the results.
#' @return A list: \code{table} (one row per design x dependent x measure) and
#'   \code{models}.
#' @examples
#' \dontrun{
#' eff <- PSOutcomeCompile(list(matched = m), dependents = c("Death_01", "LOS"))
#' eff$table
#' }
#' @export
PSOutcomeCompile <- function(designs, dependents, family = NULL, adjust = NULL, level = 0.95, Fpath = NULL,
                             lonely_psu = "adjust", PrintSummary = TRUE) {
  if (inherits(designs, "PSMatchWeight")) { designs <- list(design = designs) }
  z <- qnorm(1 - (1 - level) / 2)
  rows <- list(); models <- list()
  for (dn in names(designs)) {
    des <- designs[[dn]]; a <- des$exposure
    adj <- if (is.list(adjust)) adjust[[dn]] else adjust
    for (dep in dependents) {
      tryCatch({
      outcome_formula <- reformulate(c(a, adj), dep, env=environment(PS_Formula(c(a,adj))))
      needed <- all.vars(outcome_formula)
      d <- des$data[complete.cases(des$data[, needed, drop=FALSE]), , drop = FALSE]
      fam <- PS_family(d[[dep]], if (is.null(family)) NULL else if (length(family) == 1) family else family[[dep]])
      fit <- PS_glm_fit(outcome_formula, d, fam, ".ps_weight", des$cluster, des$strata, lonely_psu)
      survey_vcov <- fit$vcov
      weighted <- des$approach %in% c("weight", "ebal")
      inference <- if (weighted) "Survey-design covariance; weights treated as fixed" else "Survey-design covariance; matching treated as fixed"
      inference_df <- if (inherits(fit$model, "svyglm")) fit$model$df.residual else Inf
      if (weighted && !inherits(fit$model, "svyglm")) {
        # Weighted glm: the model-based covariance ignores the weighting; use the robust (sandwich) covariance.
        fit$vcov <- sandwich::vcovHC(fit$model, type = "HC1")
        inference <- sprintf("Robust (HC1 sandwich); %s weights treated as fixed",
                             if (des$approach == "ebal") "entropy balancing" else "inverse probability")
      }
      if (des$approach == "match") {
        # Weighted GLM has the same point estimates as svyglm. vcovCL accounts for
        # crossed hospital and matched-set dependence; it does NOT use survey stratification.
        robust <- PS_glm_fit(outcome_formula, d, fam, ".ps_weight")
        clusters <- data.frame(matched_set = d$subclass)
        if (!is.null(des$cluster)) clusters$hospital <- d[[des$cluster]]
        counts <- vapply(clusters, function(x) length(unique(x)), integer(1))
        if (any(counts < 2)) stop("At least two hospitals/matched sets are required for cluster inference.")
        fit$vcov <- sandwich::vcovCL(robust$model, cluster = clusters, type = "HC1", fix = FALSE)
        if (any(!is.finite(fit$vcov)) || min(eigen(fit$vcov, symmetric = TRUE, only.values = TRUE)$values) < -1e-10)
          stop("Multiway covariance is invalid; do not interpret this confidence interval.")
        inference_df <- min(counts) - 1
        inference <- if (is.null(des$cluster)) "Matched-set robust; fixed matching" else
          "Hospital + matched-set robust; survey weights, no stratification variance; fixed matching"
      }
      z <- qt(1 - (1 - level)/2, df = inference_df)
      models[[dn]][[dep]] <- fit$model
      b <- coef(fit$model)
      if (anyNA(b)) { stop(sprintf("Aliased coefficients for %s in design %s; remove collinear adjust variables.", dep, dn)) }
      target <- if (des$estimand == "ATT") d[[a]] == 1 else rep(TRUE, nrow(d))
      tt <- delete.response(terms(fit$model))
      d1 <- d[target, , drop = FALSE]; d1[[a]] <- 1
      d0 <- d1; d0[[a]] <- 0
      X1 <- model.matrix(tt, d1); X0 <- model.matrix(tt, d0)
      ww <- d$.ps_weight[target]
      linkinv <- family(fit$model)$linkinv; mu_eta <- family(fit$model)$mu.eta
      e1 <- drop(X1 %*% b); e0 <- drop(X0 %*% b)
      mu <- c(sum(ww * linkinv(e1)), sum(ww * linkinv(e0))) / sum(ww)
      J <- rbind(colSums(ww * mu_eta(e1) * X1), colSums(ww * mu_eta(e0) * X0)) / sum(ww)
      Vm <- J %*% fit$vcov[colnames(X1), colnames(X1)] %*% t(J)
      survey_Vm <- J %*% survey_vcov[colnames(X1), colnames(X1)] %*% t(J)
      measure <- function(name, est, grad, log_scale) {
        se <- sqrt(max(0, drop(t(grad) %*% Vm %*% grad)))
        survey_se <- sqrt(max(0, drop(t(grad) %*% survey_Vm %*% grad)))
        survey_z <- if (inherits(fit$model, "svyglm")) qt(1 - (1-level)/2, fit$model$df.residual) else qnorm(1 - (1-level)/2)
        transform <- if (log_scale) exp else identity
        data.frame(design = dn, approach = des$approach, estimand = des$estimand, dependent = dep, family = fam,
                   records = nrow(d), events = if (fam == "binomial") sum(d[[dep]]) else NA,
                   adjust = paste(adj, collapse = " + "), measure = name,
                   estimate = if (log_scale) exp(est) else est, lower = if (log_scale) exp(est - z * se) else est - z * se,
                   upper = if (log_scale) exp(est + z * se) else est + z * se, p_value = 2 * pt(-abs(est / se), df = inference_df),
                   survey_lower = transform(est - survey_z * survey_se), survey_upper = transform(est + survey_z * survey_se),
                   inference = inference, inference_df = inference_df,
                   inference_note = "Conditional on selected covariates, fitted propensity scores and matching; selection uncertainty omitted",
                   target_population = des$target_population,
                   n_exposed = sum(d[[a]] == 1), n_control = sum(d[[a]] == 0),
                   n_missing_outcome = nrow(des$data) - nrow(d),
                   review_status = paste(c(unique(des$diagnostics$review_status),
                     if (nrow(d) < nrow(des$data)) "Outcome missingness changes the analysis sample; review balance for this outcome"), collapse = "; "),
                   mean_exposed = mu[1], mean_unexposed = mu[2], status = "ok", error = "",
                   fit_warnings = paste(fit$model$ps_warnings, collapse = "; "), row.names = NULL)
      }
      out <- if (fam == "binomial") {
        rbind(measure("Risk difference", mu[1] - mu[2], c(1, -1), FALSE),
              measure("Risk ratio", log(mu[1] / mu[2]), c(1 / mu[1], -1 / mu[2]), TRUE),
              measure("Odds ratio", qlogis(mu[1]) - qlogis(mu[2]), c(1 / (mu[1] * (1 - mu[1])), -1 / (mu[2] * (1 - mu[2]))), TRUE))
      } else if (fam == "poisson") {
        rbind(measure("Mean difference", mu[1] - mu[2], c(1, -1), FALSE),
              measure(PS_MeasureFor("poisson"), log(mu[1] / mu[2]), c(1 / mu[1], -1 / mu[2]), TRUE))
      } else {
        measure("Mean difference", mu[1] - mu[2], c(1, -1), FALSE)
      }
      rows[[length(rows) + 1]] <- out
      }, error = function(e) {
        warning("Outcome ", dep, " failed: ", conditionMessage(e))
        rows[[length(rows) + 1]] <<- data.frame(design = dn, approach = des$approach,
          estimand = des$estimand, dependent = dep, family = NA_character_, measure = "Not estimated",
          estimate = NA_real_, lower = NA_real_, upper = NA_real_, p_value = NA_real_,
          status = "failed", error = conditionMessage(e), target_population = des$target_population,
          review_status = "Failed outcome fit; review error", stringsAsFactors = FALSE)
      })
    }
  }
  columns <- unique(unlist(lapply(rows, names)))
  table <- do.call(rbind, lapply(rows, function(r) {
    for (nm in setdiff(columns, names(r))) r[[nm]] <- NA
    r[, columns, drop = FALSE]
  }))
  key <- paste(table$design, table$measure)
  table$p_bonferroni <- ave(table$p_value, key, FUN = function(p) p.adjust(p, "bonferroni"))
  table$p_fdr <- ave(table$p_value, key, FUN = function(p) p.adjust(p, "BH"))
  if (!is.null(Fpath)) {
    dir.create(Fpath, recursive = TRUE, showWarnings = FALSE)
    for (dep in dependents) {
      sheets <- list(Effects = table[table$dependent == dep, ])
      for (dn in names(designs)) {
        sheets[[paste0(substr(dn, 1, 18), "_balance")]] <- designs[[dn]]$balance
        sheets[[paste0(substr(dn, 1, 18), "_diagnostics")]] <- designs[[dn]]$diagnostics
      }
      openxlsx::write.xlsx(sheets, file = file.path(Fpath, paste0("PS_Effects_", dep, ".xlsx")), colWidths = "auto")
    }
  }
  if (PrintSummary) { print(table[, c("design", "dependent", "measure", "estimate", "lower", "upper", "p_value")], row.names = FALSE, digits = 3) }
  invisible(list(table = table, models = models))
}

#' Balance on the records used for one outcome
#'
#' Recomputes covariate balance on the matched or weighted records that have
#' the outcome (and outcome-model predictors) recorded, standardizing by the
#' original population, because records with a missing outcome are dropped from
#' the effect estimate.
#'
#' @param DT The original data (the SMD standardization reference).
#' @param matched A \code{PSMatchWeight()} result.
#' @param dep Outcome column.
#' @param covariates Covariates to assess.
#' @param outcome_predictors Optional outcome-model predictors that must also be
#'   complete.
#' @return A list: \code{balance} (one row per covariate) and \code{summary}
#'   (one row).
#' @export
PSOutcomeBalance <- function(DT, matched, dep, covariates, outcome_predictors = NULL) {
  d <- matched$data; a <- matched$exposure
  keep <- complete.cases(d[, unique(c(dep, outcome_predictors)), drop=FALSE])
  sw <- if (is.null(matched$weights)) rep(1, nrow(DT)) else DT[[matched$weights]]
  aw <- numeric(nrow(DT)); aw[d$.ps_row_id[keep]] <- d$.ps_weight[keep]
  template <- matched$balance[, c("variable", "smd_before", "smd_after", "var_ratio_before", "var_ratio_after",
    "ks_before", "ks_after", "mean_control_before", "mean_exposed_before", "mean_control_after", "mean_exposed_after"), drop = FALSE]
  note <- ""
  b <- tryCatch({
    if (!all(c(0, 1) %in% d[[a]][keep])) stop("Both exposure groups are required for outcome-sample balance.")
    cobalt::bal.tab(x = DT[, covariates, drop = FALSE], treat = DT[[a]], weights = aw/sw,
      s.weights = sw, estimand = matched$estimand, stats = c("m", "v", "ks"), binary = "std", un = TRUE, quick = FALSE)$Balance
  }, error = function(e) { note <<- conditionMessage(e); NULL })
  columns <- c(smd_after="Diff.Adj", var_ratio_after="V.Ratio.Adj", ks_after="KS.Adj",
               mean_control_after="M.0.Adj", mean_exposed_after="M.1.Adj")
  for (nm in names(columns)) template[[nm]] <- if (is.null(b) || !columns[[nm]] %in% names(b)) NA_real_ else b[match(template$variable, rownames(b)), columns[[nm]]]
  template$balanced <- is.finite(template$smd_after) & abs(template$smd_after) < .1
  template$diagnostic_error <- note
  ess <- function(g) { w <- aw[DT[[a]] == g]; if (sum(w^2) > 0) sum(w)^2/sum(w^2) else 0 }
  max_smd <- if (any(is.finite(template$smd_after))) max(abs(template$smd_after), na.rm=TRUE) else NA_real_
  summary <- data.frame(outcome_max_smd_after=max_smd, outcome_balanced_all_assessed=all(template$balanced),
    outcome_n_imbalanced=sum(!template$balanced), outcome_ess_exposed=ess(1), outcome_ess_controls=ess(0),
    outcome_balance_note=if (nzchar(note)) note else if (all(template$balanced)) "SMD threshold passed for outcome sample" else "Outcome sample has residual imbalance")
  for (j in seq_len(nrow(template))) for (nm in c(names(columns), "balanced"))
    summary[[paste0("outcome_", nm, "_", template$variable[j])]] <- template[[nm]][j]
  list(balance=template, summary=summary)
}

#' Entropy balancing weights and balance
#'
#' Weights chosen with \code{WeightIt::weightit(method = "ebal")} so the
#' weighted covariate moments of the controls (ATT) or both groups (ATE) match
#' the target, while staying close to the survey weights. No propensity score is
#' fitted.
#'
#' @param data Data frame.
#' @param exposure The 0/1 exposure column.
#' @param covariates Covariates to balance.
#' @param weights,cluster,strata Design columns, or \code{NULL}.
#' @param estimand \code{"ATT"} or \code{"ATE"}.
#' @param moments 1 balances means; 2 also balances variances.
#' @param balance_covariates Covariates whose balance is assessed.
#' @param balance_threshold Absolute SMD below which a covariate is balanced.
#' @return A \code{"PSMatchWeight"}-style list with \code{approach = "ebal"}
#'   that \code{PSOutcomeCompile()} accepts.
#' @examples
#' \dontrun{
#' eb <- EBWeight(model_dt, "sud_any_01", c("AGE", "sex_01"), weights = "DISCWT")
#' }
#' @export
EBWeight <- function(data, exposure, covariates, weights = NULL, cluster = NULL, strata = NULL, estimand = c("ATT", "ATE"),
                     moments = 1, balance_covariates = covariates, balance_threshold = 0.1) {
  estimand <- match.arg(estimand)
  if (!length(covariates)) stop("Entropy balancing needs at least one covariate.")
  if (!moments %in% c(1, 2)) stop("moments must be 1 (means) or 2 (means and variances).")
  d <- as.data.frame(data); d$.ps_row_id <- seq_len(nrow(d))
  PS_ValidateComplete(d, c(exposure, PS_BaseVars(covariates), weights, cluster, strata,
                           PS_BaseVars(balance_covariates)), "Entropy-balancing model")
  if (!all(d[[exposure]] %in% c(0, 1)) || length(unique(d[[exposure]])) != 2L)
    stop("Entropy-balancing exposure must contain both 0 and 1 and no other values.", call. = FALSE)
  sw <- if (is.null(weights)) rep(1, nrow(d)) else d[[weights]]
  if (any(!is.finite(sw) | sw <= 0)) stop("Entropy-balancing weights must be finite and positive.", call. = FALSE)
  f <- reformulate(covariates, exposure, env=environment(PS_Formula(covariates)))
  fit_warnings <- character()
  object <- withCallingHandlers(WeightIt::weightit(f, data = d, method = "ebal", estimand = estimand, s.weights = sw, moments = moments),
                                warning = function(w) { fit_warnings <<- c(fit_warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
  if (any(!is.finite(object$weights))) stop("Entropy balancing returned non-finite weights; check the covariates for separation.")
  out <- d
  out$.ps_weight <- object$weights * sw
  bd <- PS_BalanceDiagnostics(d, exposure, sw, out$.ps_weight, unique(c(PS_BaseVars(balance_covariates), PS_BaseVars(covariates))), PS_BaseVars(covariates), estimand,
                              balance_threshold, fit_warnings = fit_warnings)
  structure(list(approach = "ebal", estimand = estimand, target_population = if (estimand == "ATT") "Exposed population" else estimand,
      data = out, balance = bd$balance, diagnostics = bd$diagnostics, max_smd_after = bd$diagnostics$max_abs_smd_after[1],
      love_plot = PS_LovePlot(bd$bal_tab, balance_threshold),
      ess = setNames(bd$diagnostics$ess_after, bd$diagnostics$group), n_unmatched = NA, object = object, exposure = exposure,
      covariates = PS_BaseVars(covariates), moments = moments, formula = f, fit_warnings = fit_warnings,
      weights = weights, cluster = cluster, strata = strata), class = "PSMatchWeight")
}
