#' Extract stable columns from a survey design
#'
#' @param design A \code{survey.design} object.
#' @return A list containing \code{data}, \code{cluster}, \code{strata} and
#'   \code{weights}. Design values are copied into collision-safe internal
#'   columns, so the original data symbol need not remain in scope.
#' @keywords internal
#' @noRd
svyVarSel_design_data <- function(design) {
  if (!inherits(design, c("survey.design", "survey.design2")))
    stop("design must be created by survey::svydesign().", call. = FALSE)
  data <- as.data.frame(design$variables)
  internal_name <- function(stem) {
    candidate <- stem
    while (candidate %in% names(data)) candidate <- paste0(candidate, "_")
    candidate
  }
  weights <- internal_name(".svyVarSel_weight")
  data[[weights]] <- as.numeric(stats::weights(design))
  cluster <- NULL
  id_vars <- tryCatch(all.vars(design$call$id), error = function(e) character())
  if (length(id_vars)) {
    cluster <- internal_name(".svyVarSel_cluster")
    cl <- as.data.frame(design$cluster)
    data[[cluster]] <- if (ncol(cl) == 1L) cl[[1]] else
      do.call(interaction, c(unname(as.list(cl)), list(drop = TRUE, lex.order = TRUE)))
  }
  strata <- NULL
  strata_vars <- tryCatch(all.vars(design$call$strata), error = function(e) character())
  if (length(strata_vars)) {
    strata <- internal_name(".svyVarSel_strata")
    st <- as.data.frame(design$strata)
    data[[strata]] <- if (ncol(st) == 1L) st[[1]] else
      do.call(interaction, c(unname(as.list(st)), list(drop = TRUE, lex.order = TRUE)))
  }
  list(data = data, cluster = cluster, strata = strata, weights = weights)
}

#' Survey-weighted LASSO with replicate-weight cross-validation
#'
#' A copy of \code{svyVarSel::wlasso} (svyVarSel 1.0.1; Iparragirre, Lumley,
#' Barrio and Arostegui; GPL (>= 3)) with two changes. (1) The final
#' single-lambda \pkg{glmnet} refit at lambda.min is replaced by the matching
#' column of the full lambda path; that refit cold-starts at a small lambda and
#' can fail to converge, leaving an empty \code{$model$min}. (2)
#' \code{family = "poisson"} is supported: the cross-validated error is computed
#' by svyVarSel's internal error function with a loss that adds the Poisson
#' deviance (gaussian and binomial losses are unchanged).
#'
#' The internal helpers \code{replicate.weights} and \code{error.f} of
#' \pkg{svyVarSel} are used, so \pkg{svyVarSel} must be installed;
#' \code{svyVarSel::wlasso.plot()} plots the result.
#'
#' @param data Data frame with the outcome, predictors and design columns.
#'   Ignored when \code{design} is given.
#' @param col.y Name of the outcome column.
#' @param col.x Character vector of candidate predictor columns.
#' @param cluster Name of the cluster (primary sampling unit) column, or
#'   \code{NULL}.
#' @param strata Name of the strata column, or \code{NULL}.
#' @param weights Name of the sampling-weight column.
#' @param design Optional \code{survey::svydesign} object used instead of
#'   \code{data}, \code{cluster}, \code{strata} and \code{weights}.
#' @param family \code{"gaussian"}, \code{"binomial"} or \code{"poisson"}.
#' @param lambda.grid Optional decreasing lambda grid; \code{NULL} uses
#'   glmnet's own path.
#' @param method Replicate-weight method: \code{"dCV"}, \code{"JKn"},
#'   \code{"bootstrap"}, \code{"subbootstrap"}, \code{"BRR"}, \code{"split"} or
#'   \code{"extrapolation"}.
#' @param k Number of folds for \code{method = "dCV"}.
#' @param R Number of repeats for \code{"dCV"}, \code{"split"} and
#'   \code{"extrapolation"}.
#' @param B Number of replicates for the bootstrap methods.
#' @param dCV.sw.test Logical; for \code{"dCV"}, use the sampling weights of the
#'   test set to compute the error.
#' @param train.prob Share of clusters in the training set for \code{"split"}
#'   and \code{"extrapolation"}.
#' @param method.split Replicate method used within \code{"split"}.
#' @param print.rw Logical; add the data with replicate-weight columns to the
#'   result as \code{data.rw}.
#' @return An object of class \code{"wlasso"}: \code{lambda} (\code{grid},
#'   \code{min}), \code{error} (\code{average}, \code{all}), \code{model}
#'   (\code{grid}: intercepts, coefficients and degrees of freedom for every
#'   lambda; \code{min}: the lambda.min column of the path) and \code{call}.
#' @seealso \code{\link{svyVarSel_wlasso_tune}}, \code{\link{LassoCompile}}
#' @examples
#' \dontrun{
#' load(system.file("extdata", "egs_misc_1621.rda", package = "dataMineR"))
#' d <- egs_misc_1621
#' d$strata_id <- paste(d$YEAR, d$NIS_STRATUM)
#' d$cluster_id <- paste(d$YEAR, d$HOSP_NIS)
#' fit <- svyVarSel_wlasso_adjusted(data = d, col.y = "Death_01",
#'   col.x = c("AGE", "sex_01", "elix_chf", "elix_fed"), family = "binomial",
#'   strata = "strata_id", weights = "DISCWT", cluster = "cluster_id",
#'   method = "dCV", k = 5, R = 1)
#' }
#' @export
svyVarSel_wlasso_adjusted <- function(data = NULL, col.y = NULL, col.x = NULL, cluster = NULL, strata = NULL, weights = NULL,
                            design = NULL, family = c("gaussian", "binomial", "poisson"), lambda.grid = NULL,
                            method = c("dCV", "JKn", "bootstrap", "subbootstrap", "BRR", "split", "extrapolation"),
                            k = 10, R = 1, B = 200, dCV.sw.test = FALSE, train.prob = 0.7,
                            method.split = c("dCV", "bootstrap", "subbootstrap"), print.rw = FALSE) {
  family <- match.arg(family)
  method <- match.arg(method)
  method.split <- match.arg(method.split)
  if (is.null(data) & is.null(design)) {
    stop("Information about either the data set ('data') or the sampling design ('design') needed.")
  }
  if (method == "split") {
    if (is.null(train.prob)) stop("Selected replicate weights method: 'split'.\nPlease, set a value between 0 and 1 for the argument 'train.prob'.")
    if (train.prob < 0 | train.prob > 1) stop("Selected replicate weights method: 'split'.\nPlease, set a value between 0 and 1 for the argument 'train.prob'.")
    if (length(method.split) != 1) stop("Selected replicate weights method: 'split'.\nPlease, set a valid method for the argument 'method.split'. Choose between: 'dCV', 'bootstrap' or 'subbootstrap'.")
  }
  if (method == "extrapolation") {
    if (is.null(train.prob)) stop("Selected replicate weights method: 'extrapolation'.\nPlease, set a value between 0 and 1 for the argument 'train.prob'.")
    if (train.prob < 0 | train.prob > 1) stop("Selected replicate weights method: 'extrapolation'.\nPlease, set a value between 0 and 1 for the argument 'train.prob'.")
  }
  if (method %in% c("JKn", "bootstrap", "subbootstrap", "BRR")) {
    if (R != 1) message("Selected method:", method, ". For this method, R = 1. Thus, the argument R =", R, "has been ignored.")
  }
  if (method %in% c("dCV", "split", "extrapolation")) {
    if (R != round(R)) stop("The argument 'R' must be an integer greater or equal to 1. R=", R, " is not an integer.\nPlease, set a valid value for 'R' or skip the argument to select the default option R=1.")
    if (R < 1) stop("The argument 'R' must be an integer greater or equal to 1. R=", R, " lower than 1.\nPlease, set a valid value for 'R' or skip the argument to select the default option R=1.")
  }
  if (method != "dCV") {
    if (!is.null(k) & k != 10) message("Selected method:", method, ". The argument k =", k, "is not needed and, hence, has been ignored.")
  }
  if (method == "dCV") {
    if (k != round(k)) stop("The argument 'k' must be an integer. k=", k, " is not an integer.\nPlease, set a valid value for 'k' or skip the argument to select the default option k=10.")
    if (k < 1) stop("The argument 'k' must be a positive integer. k=", k, " is not a positive integer.\nPlease, set a valid value for 'k' or skip the argument to select the default option k=10.")
  }
  if (!(method %in% c("bootstrap", "subbootstrap"))) {
    if (!is.null(B) & B != 200) message("Selected method:", method, ". The argument B=", B, " is not needed and, hence, has been ignored.")
  }
  if (method %in% c("bootstrap", "subbootstrap")) {
    if (B != round(B)) stop("The argument 'B' must be an integer. B=", B, " is not an integer.\nPlease, set a valid value for 'B' or skip the argument to select the default option B=200.")
    if (B < 1) stop("The argument 'B' must be a positive integer. B=", B, " is not a positive integer.\nPlease, set a valid value for 'B' or skip the argument to select the default option B=200.")
  }
  # Copy design values into stable internal columns rather than evaluating the
  # original data expression, which may no longer exist in this environment.
  if (!is.null(design)) {
    design_data <- svyVarSel_design_data(design)
    data <- design_data$data
    cluster <- design_data$cluster
    strata <- design_data$strata
    weights <- design_data$weights
  }

  newdata <- utils::getFromNamespace("replicate.weights", "svyVarSel")(data = data, method = method, cluster = cluster, strata = strata,
                                           weights = weights, k = k, R = R, B = B, train.prob = train.prob,
                                           method.split = method.split, rw.test = TRUE, dCV.sw.test = dCV.sw.test)

  if (is.null(lambda.grid)) {
    model.orig <- glmnet::glmnet(y = as.numeric(newdata[, col.y]), x = as.matrix(newdata[, col.x]),
                                 weights = as.numeric(newdata[, weights]), family = family)
    lambda.grid <- model.orig$lambda
  } else {
    model.orig <- glmnet::glmnet(y = as.numeric(newdata[, col.y]), x = as.matrix(newdata[, col.x]),
                                 weights = as.numeric(newdata[, weights]), family = family, lambda = lambda.grid)
  }

  rwtraincols <- grep("_train", colnames(newdata))
  l.yhat <- list()
  for (col.w in rwtraincols) {
    model <- glmnet::glmnet(y = as.numeric(newdata[, col.y]), x = as.matrix(newdata[, col.x]),
                            weights = as.numeric(newdata[, col.w]), lambda = lambda.grid, family = family)
    yhat <- predict(model, newx = as.matrix(newdata[, col.x]), type = "response")
    l.yhat[[length(l.yhat) + 1]] <- yhat
    names(l.yhat)[[length(l.yhat)]] <- paste0("yhat_", colnames(newdata)[col.w])
  }

  error_f <- utils::getFromNamespace("error.f", "svyVarSel")
  environment(error_f) <- list2env(list(loss.f = svyVarSel_loss), parent = asNamespace("svyVarSel"))
  error <- error_f(data = newdata, l.yhat = l.yhat, method = method, cv.error.ind = dCV.sw.test,
                   R = R, k = k, B = B, col.y = col.y, family = family, weights = weights)
  mean.error <- apply(error, 2, mean)
  lambda.min <- lambda.grid[which.min(mean.error)]

  j <- match(lambda.min, model.orig$lambda)
  if (is.na(j)) { stop("Selected lambda is missing from the fitted path.") }
  model <- list(a0 = model.orig$a0[j], beta = model.orig$beta[, j, drop = FALSE], df = model.orig$df[j])
  result <- list()
  result$lambda <- list(grid = lambda.grid, min = lambda.min)
  result$error <- list(average = mean.error, all = error)
  result$model <- list()
  result$model$grid <- list(a0 = model.orig$a0, beta = model.orig$beta, df = model.orig$df)
  result$model$min <- list(a0 = model$a0, beta = model$beta, df = model$df)
  result$call <- match.call()
  if (print.rw == TRUE) { result$data.rw <- newdata }   # data with the replicate-weight columns
  class(result) <- "wlasso"
  return(result)
}

#' Cross-validation loss for one prediction per record
#'
#' Used by \code{svyVarSel_wlasso_adjusted} in place of svyVarSel's internal
#' loss: squared error (gaussian), minus the log of the predicted probability of
#' the observed class (binomial) and the Poisson deviance
#' \code{2 * (y log(y / mu) - (y - mu))} (poisson).
#'
#' @param y.est Predicted means (\code{type = "response"}).
#' @param y Observed outcomes.
#' @param family \code{"gaussian"}, \code{"binomial"} or \code{"poisson"}.
#' @return A numeric vector of losses, one per record.
#' @keywords internal
#' @noRd
svyVarSel_loss <- function(y.est, y, family = c("gaussian", "binomial", "poisson")) {
  family <- match.arg(family)
  if (family == "gaussian") { return((y.est - y)^2) }
  if (family == "binomial") { return(ifelse(y == 1, -log(y.est), -log(1 - y.est))) }
  2 * (ifelse(y > 0, y * log(y / y.est), 0) - (y - y.est))
}

#' Choose lambda from a wlasso fit
#'
#' Chooses the lambda with the lowest average cross-validated error
#' (\code{"min"}), or the largest lambda whose error is within one standard error
#' of that minimum (\code{"1se"}). The standard error is the spread of the
#' replicate errors divided by the square root of the number of replicates; with
#' a single replicate the 1-SE choice falls back to lambda.min.
#'
#' @param fit A \code{"wlasso"} object.
#' @param lambda_rule \code{"min"} or \code{"1se"}.
#' @param PrintSummary Logical; print the summary, and a message when lambda.min
#'   is at the end of the grid.
#' @return Invisibly, a list: \code{lambda_rule}, \code{grid}, \code{cv_mean},
#'   \code{cv_se}, \code{j_min}, \code{j_1se}, \code{j} (position of the chosen
#'   lambda), \code{lambda}, \code{at_boundary}, \code{summary} (one row per
#'   rule) and \code{checks} (lambda and number of variables for both rules, and
#'   whether lambda.min is at the end of the grid).
#' @examples
#' \dontrun{
#' sel <- svyVarSel_lambda_select(fit, lambda_rule = "min")
#' sel$summary
#' }
#' @export
svyVarSel_lambda_select <- function(fit, lambda_rule = c("min", "1se"), PrintSummary = TRUE) {
  if (!inherits(fit, "wlasso")) { stop("Please, insert an object of class 'wlasso'.") }
  lambda_rule <- match.arg(lambda_rule)
  grid <- fit$lambda$grid
  cv_mean <- fit$error$average
  n_replicates <- nrow(fit$error$all)
  cv_se <- if (n_replicates > 1) { apply(fit$error$all, 2, sd) / sqrt(n_replicates) } else { rep(NA_real_, length(cv_mean)) }
  j_min <- which.min(cv_mean)
  j_1se <- if (is.na(cv_se[j_min])) { j_min } else { min(which(cv_mean <= cv_mean[j_min] + cv_se[j_min])) }
  if (n_replicates < 2) { warning("Only one replicate row in $error$all: the 1-SE rule is undefined, so it equals lambda.min.") }
  j <- if (lambda_rule == "1se") { j_1se } else { j_min }
  at_boundary <- j_min == length(grid)
  summary <- data.frame(rule = c("min", "1se"), position = c(j_min, j_1se), lambda = grid[c(j_min, j_1se)],
                        n_vars = fit$model$grid$df[c(j_min, j_1se)], cv_error = cv_mean[c(j_min, j_1se)],
                        at_grid_end = c(j_min, j_1se) == length(grid), selected = c("min", "1se") == lambda_rule)
  checks <- list(lambda_min = grid[j_min], n_vars_min = fit$model$grid$df[j_min],
                 lambda_1se = grid[j_1se], n_vars_1se = fit$model$grid$df[j_1se],
                 lambda_min_at_grid_end = at_boundary)
  if (PrintSummary) {
    print(summary, row.names = FALSE)
    if (at_boundary) { message("lambda.min is the smallest lambda in the grid: CV error is still falling (or flat) at the end of the path.") }
  }
  invisible(list(lambda_rule = lambda_rule, grid = grid, cv_mean = cv_mean, cv_se = cv_se, j_min = j_min, j_1se = j_1se,
                 j = j, lambda = grid[j], at_boundary = at_boundary, summary = summary, checks = checks))
}

#' LASSO coefficient paths
#'
#' Plots the coefficient path of an existing fit (no refitting), highlighting
#' the predictors selected at position \code{j}. Coefficients are standardized
#' (coefficient x weighted predictor SD) unless \code{standardize = FALSE}.
#' Supply the exact training rows and predictor coding used for the fit.
#'
#' @param fit A \code{"wlasso"} object.
#' @param data The training data used for the fit.
#' @param weights Name of the weight column, or \code{NULL} for equal weights.
#' @param j Position on the lambda path; defaults to lambda.min.
#' @param standardize Logical; plot standardized coefficients.
#' @param Fpath Optional folder to save the plot to.
#' @param FileName File name used with \code{Fpath}.
#' @param PrintPlot Logical; print the plot.
#' @return Invisibly, a list: \code{plot}, \code{paths} (long table of the
#'   path), \code{selected}, \code{correlation} (weighted predictor
#'   correlations), \code{lambda} and \code{omitted_constants}.
#' @examples
#' \dontrun{
#' LassoCoefficientPathPlot(fit, d, weights = "DISCWT")
#' }
#' @export
LassoCoefficientPathPlot <- function(fit, data, weights = NULL,
    j = match(fit$lambda$min, fit$lambda$grid), standardize = TRUE,
    Fpath = NULL, FileName = "LASSO_CoefficientPaths.pdf", PrintPlot = is.null(Fpath)) {
  beta <- as.matrix(fit$model$grid$beta); grid <- fit$lambda$grid
  if (!is.numeric(beta) || !nrow(beta) || ncol(beta) != length(grid) ||
      any(!is.finite(beta)) || any(!is.finite(grid) | grid <= 0) || anyDuplicated(grid))
    stop("A finite coefficient path and matching positive lambda grid are required.")
  if (length(j)!=1L || is.na(j) || j!=as.integer(j) || j<1 || j>ncol(beta)) stop("Invalid lambda position j.")
  vars <- rownames(beta); d <- as.data.frame(data)
  if (is.null(vars) || anyDuplicated(vars) || !all(vars %in% names(d))) stop("Training data must contain each named path predictor.")
  if (!all(vapply(d[,vars,drop=FALSE],is.numeric,logical(1)))) stop("Predictors must retain the numeric coding used in the LASSO fit.")
  x <- as.matrix(d[,vars,drop=FALSE])
  if (nrow(x)<2 || any(!is.finite(x))) stop("Supply complete, finite training predictors; rows are not automatically removed.")
  if (!is.null(weights) && (length(weights)!=1L || !weights %in% names(d))) stop("weights must name a training-data column.")
  w <- if(is.null(weights)) rep(1,nrow(d)) else d[[weights]]
  if(!is.numeric(w) || any(!is.finite(w) | w<=0)) stop("Weights must be finite and positive.")
  cv <- stats::cov.wt(x,wt=w,method="ML")$cov
  sd <- sqrt(diag(cv)); varying <- sd>0
  cor <- matrix(numeric(0),0,0)
  if(any(varying)) cor <- cv[varying,varying,drop=FALSE]/outer(sd[varying],sd[varying])
  cor[] <- pmax(-1,pmin(1,cor))
  selected <- vars[beta[,j]!=0]
  plotted <- if(standardize) sweep(beta,1,sd,"*") else beta
  paths <- data.frame(variable=rep(vars,times=length(grid)),lambda=rep(grid,each=length(vars)),
    coefficient=as.vector(plotted),original_coefficient=as.vector(beta))
  paths$selected <- paths$variable %in% selected
  p <- ggplot2::ggplot(paths,ggplot2::aes(x=log(lambda),y=coefficient,group=variable)) +
    ggplot2::geom_hline(yintercept=0,colour="grey70",linewidth=.3) +
    ggplot2::geom_line(data=paths[!paths$selected,,drop=FALSE],colour="grey75",linewidth=.4) +
    ggplot2::geom_line(data=paths[paths$selected,,drop=FALSE],ggplot2::aes(colour=variable),linewidth=.65) +
    ggplot2::geom_vline(xintercept=log(grid[j]),linetype=2) +
    ggplot2::labs(title="LASSO coefficient paths",subtitle=sprintf("Selected lambda = %.4g; %d predictors selected",grid[j],length(selected)),
      x="log(lambda): stronger penalization to the right",
      y=if(standardize) "Coefficient per weighted SD of predictor" else "Coefficient (original units)",
      caption="Gray: excluded at the selected lambda. Paths describe penalization, not proof of masking.") +
    ggplot2::theme_minimal() + ggplot2::theme(legend.position="none")
  labels <- paths[paths$selected & paths$lambda==grid[j],,drop=FALSE]
  if(nrow(labels)) {
    if(requireNamespace("ggrepel",quietly=TRUE)) p <- p + ggrepel::geom_text_repel(data=labels,
      ggplot2::aes(label=variable,colour=variable),seed=1621,max.overlaps=Inf,size=3)
    else p <- p + ggplot2::geom_text(data=labels,ggplot2::aes(label=variable,colour=variable),hjust=0,size=3)
  }
  if(!is.null(Fpath)) {
    dir.create(Fpath,recursive=TRUE,showWarnings=FALSE)
    ggplot2::ggsave(file.path(Fpath,FileName),plot=p,width=12,height=8)
  }
  # Preview failures must not abort an already completed file export.
  if(PrintPlot) tryCatch(print(p), error=function(e) {
    if (!grepl("Viewport has zero dimension", conditionMessage(e), fixed=TRUE)) stop(e)
    message("Plot preview skipped: the graphics pane is too small. Enlarge it and print the returned $plot; any requested PDF has already been saved.")
  })
  invisible(list(plot=p,paths=paths,selected=selected,correlation=cor,lambda=grid[j],
    omitted_constants=vars[!varying]))
}

#' Selected LASSO predictors and their most correlated alternatives
#'
#' A heatmap of weighted Pearson correlations between the selected predictors
#' and, for each, up to \code{partners_per_variable} excluded predictors with
#' \eqn{|r| \ge} \code{partner_cor}. Rows are clustered on \code{1 - |r|}.
#' Correlation shows overlapping information; it does not prove why a predictor
#' was excluded.
#'
#' @inheritParams LassoCoefficientPathPlot
#' @param partner_cor Minimum absolute correlation for an excluded alternative.
#' @param partners_per_variable Maximum number of alternatives per selected
#'   predictor.
#' @return Invisibly, a list: \code{plot}, \code{correlation}, \code{partners}
#'   (selected predictor, alternative, correlation), \code{selected},
#'   \code{data} and \code{omitted_constants}.
#' @examples
#' \dontrun{
#' LassoCorrelationPlot(fit, d, weights = "DISCWT", partner_cor = 0.5)
#' }
#' @export
LassoCorrelationPlot <- function(fit, data, weights = NULL,
    j = match(fit$lambda$min, fit$lambda$grid), partner_cor = 0.5,
    partners_per_variable = 2, Fpath = NULL,
    FileName = "LASSO_PredictorCorrelations.pdf", PrintPlot = is.null(Fpath)) {
  if(length(partner_cor)!=1L || !is.finite(partner_cor) || partner_cor<0 || partner_cor>1)
    stop("partner_cor must be between 0 and 1.")
  if(length(partners_per_variable)!=1L || !is.finite(partners_per_variable) ||
     partners_per_variable<0 || partners_per_variable!=floor(partners_per_variable))
    stop("partners_per_variable must be a nonnegative integer.")
  info <- LassoCoefficientPathPlot(fit,data,weights,j,PrintPlot=FALSE)
  cor <- info$correlation; selected <- intersect(info$selected,rownames(cor))
  partners <- data.frame(selected=character(),alternative=character(),correlation=numeric())
  for(v in selected) {
    candidates <- setdiff(colnames(cor),info$selected)
    candidates <- candidates[abs(cor[v,candidates])>=partner_cor]
    candidates <- head(candidates[order(-abs(cor[v,candidates]))],partners_per_variable)
    if(length(candidates)) partners <- rbind(partners,data.frame(selected=v,alternative=candidates,correlation=as.numeric(cor[v,candidates])))
  }
  vars <- unique(c(selected,partners$alternative))
  matrix <- cor[vars,vars,drop=FALSE]
  if(length(vars)>1) {
    distance <- 1-abs(matrix);diag(distance)<-0
    vars <- vars[stats::hclust(stats::as.dist(distance),method="average")$order]
    matrix <- matrix[vars,vars,drop=FALSE]
  }
  if(length(vars)) {
    tab <- expand.grid(row=vars,column=vars,stringsAsFactors=FALSE)
    tab$correlation <- as.vector(matrix)
    tab$row <- factor(tab$row,levels=rev(vars));tab$column <- factor(tab$column,levels=vars)
    labels <- setNames(paste0(vars,ifelse(vars %in% info$selected," *","")),vars)
    p <- ggplot2::ggplot(tab,ggplot2::aes(column,row,fill=correlation)) + ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(low="#2166AC",mid="white",high="#B2182B",midpoint=0,limits=c(-1,1),name="r") +
      ggplot2::scale_x_discrete(labels=labels) + ggplot2::scale_y_discrete(labels=labels) +
      ggplot2::coord_fixed() + ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x=ggplot2::element_text(angle=60,hjust=1),panel.grid=ggplot2::element_blank())
  } else {
    tab <- data.frame()
    p <- ggplot2::ggplot() + ggplot2::annotate("text",x=0,y=0,label="No varying predictors selected at this lambda") + ggplot2::theme_void()
  }
  p <- p + ggplot2::labs(title="LASSO predictor correlations",subtitle=sprintf("Selected lambda = %.4g; * = selected predictor",info$lambda),x=NULL,y=NULL,
    caption=sprintf("Weighted Pearson r; up to %d excluded alternatives per selected predictor with |r| >= %.2f.\nCorrelation indicates overlapping information, not proof of masking. Constant predictors omitted: %d.",partners_per_variable,partner_cor,length(info$omitted_constants)))
  if(!is.null(Fpath)) {
    dir.create(Fpath,recursive=TRUE,showWarnings=FALSE)
    ggplot2::ggsave(file.path(Fpath,FileName),plot=p,width=max(9,.28*length(vars)+4),height=max(8,.28*length(vars)+4),limitsize=FALSE)
  }
  # Preview failures must not abort an already completed file export.
  if(PrintPlot) tryCatch(print(p), error=function(e) {
    if (!grepl("Viewport has zero dimension", conditionMessage(e), fixed=TRUE)) stop(e)
    message("Plot preview skipped: the graphics pane is too small. Enlarge it and print the returned $plot; any requested PDF has already been saved.")
  })
  invisible(list(plot=p,correlation=matrix,partners=partners,selected=selected,data=tab,omitted_constants=info$omitted_constants))
}

#' Coefficients of a wlasso fit at one position on the path
#'
#' @param fit A \code{"wlasso"} object.
#' @param j Position on the lambda path (for example
#'   \code{svyVarSel_lambda_select()$j}); defaults to lambda.min.
#' @return A list: \code{coefs} (all predictors), \code{selected} (nonzero
#'   coefficients), \code{excluded} (names of the zero ones) and \code{lambda}.
#'   Coefficients are on the model's link scale (log odds for binomial).
#' @examples
#' \dontrun{
#' svyVarSel_selected_coefs(fit)$selected
#' }
#' @export
svyVarSel_selected_coefs <- function(fit, j = match(fit$lambda$min, fit$lambda$grid)) {
  if (!inherits(fit, "wlasso")) { stop("Please, insert an object of class 'wlasso'.") }
  if (length(j) != 1L || is.na(j) || j < 1 || j > ncol(fit$model$grid$beta)) { stop("j must be one position on the lambda path.") }
  coefs <- as.matrix(fit$model$grid$beta[, j, drop = FALSE])
  colnames(coefs) <- "coef"
  list(coefs = coefs, selected = coefs[coefs[, "coef"] != 0, , drop = FALSE],
       excluded = rownames(coefs)[coefs[, "coef"] == 0], lambda = fit$lambda$grid[j])
}

#' Why were predictors excluded from a wlasso model?
#'
#' For each predictor with a zero coefficient at position \code{j}, describes
#' its group size and events (0/1 predictors), where it is nonzero on the path,
#' its unadjusted association with the outcome, how well the other predictors
#' reproduce it (R-squared), its correlated partners, and its coefficient when
#' the same path is refitted without those partners. The first reason that
#' applies is reported: it enters once its partners are removed (overlap); a
#' small group (binomial outcomes, fewer than \code{min_events} events); largely
#' reproduced by the other predictors (R-squared at least 0.5); no association
#' even unadjusted; otherwise its association is explained by the other
#' predictors.
#'
#' @param fit A \code{"wlasso"} object.
#' @param data The training data used for the fit.
#' @param outcome Name of the outcome column.
#' @param weights Name of the weight column.
#' @param j Position on the lambda path; defaults to lambda.min.
#' @param variables Predictors to explain; defaults to all excluded ones.
#' @param family \code{"binomial"}, \code{"poisson"} or \code{"gaussian"}.
#' @param partner_cor Minimum absolute correlation for a partner.
#' @param min_events Minimum events for a 0/1 predictor before it counts as a
#'   small group.
#' @param PrintSummary Logical; print the summary table.
#' @return Invisibly, a list: \code{details} (one list per predictor) and
#'   \code{summary} (one row per predictor with the reason).
#' @examples
#' \dontrun{
#' svyVarSel_excluded_reasons(fit, d, outcome = "Death_01", weights = "DISCWT")$summary
#' }
#' @export
svyVarSel_excluded_reasons <- function(fit, data, outcome, weights, j = match(fit$lambda$min, fit$lambda$grid),
                                       variables = NULL, family = "binomial", partner_cor = 0.3, min_events = 100,
                                       PrintSummary = TRUE) {
  if (!inherits(fit, "wlasso")) { stop("Please, insert an object of class 'wlasso'.") }
  beta <- fit$model$grid$beta; grid <- fit$lambda$grid; predictors <- rownames(beta)
  if (length(j) != 1L || is.na(j) || j < 1 || j > ncol(beta)) { stop("j must be one position on the lambda path.") }
  if (is.null(variables)) { variables <- predictors[beta[, j] == 0] }
  if (!length(variables)) {
    message("No predictors are excluded at this lambda.")
    return(invisible(list(details = list(), summary = data.frame())))
  }
  data <- as.data.frame(data)
  x <- as.matrix(data[, predictors, drop = FALSE]); y <- data[[outcome]]; w <- data[[weights]]
  cor_x <- stats::cov.wt(x, wt = w, cor = TRUE)$cor
  is_binary <- function(v) all(v %in% c(0, 1))
  details <- list()
  for (v in variables) {
    xv <- x[, v]; binary <- is_binary(xv)
    group <- if (binary) c(patients = sum(xv == 1), events = sum(y[xv == 1]), event_rate = mean(y[xv == 1])) else
      c(patients = length(xv), events = NA, event_rate = NA)
    nz <- which(beta[v, ] != 0)
    path <- data.frame(position = nz, lambda = grid[nz], coef = as.numeric(beta[v, nz]), row.names = NULL)
    glm_family <- switch(family, binomial = quasibinomial(), poisson = quasipoisson(), gaussian())
    uni <- summary(suppressWarnings(glm(y ~ xv, weights = w, family = glm_family)))$coefficients
    unadjusted <- c(estimate = uni["xv", 1], OR = if (family == "binomial") exp(uni["xv", 1]) else NA, p = uni["xv", 4])
    others <- setdiff(predictors, v)
    ls_fit <- lm.wfit(cbind(1, x[, others, drop = FALSE]), xv, w)
    r2 <- 1 - sum(w * ls_fit$residuals^2) / sum(w * (xv - weighted.mean(xv, w))^2)
    r <- cor_x[v, others]; partners <- names(sort(abs(r[abs(r) >= partner_cor]), decreasing = TRUE))
    partner_table <- data.frame(partner = partners, cor = as.numeric(r[partners]),
                                selected_at_j = as.numeric(beta[partners, j]) != 0, row.names = NULL)
    overlap <- NULL
    if (length(partners) && binary) {
      partner_table$share_of_its_patients <- vapply(partners, function(p) if (is_binary(x[, p])) mean(x[xv == 1, p] == 1) else NA_real_, numeric(1))
      binary_partners <- partners[vapply(partners, function(p) is_binary(x[, p]), logical(1))]
      if (length(binary_partners)) {
        overlap <- table(xv, any_partner = as.integer(rowSums(x[, binary_partners, drop = FALSE]) > 0), dnn = c(v, "any_partner"))
      }
    }
    coef_without <- NA_real_
    if (length(partners)) {
      path_without <- glmnet::glmnet(x = x[, setdiff(predictors, partners), drop = FALSE], y = y, weights = w,
                                     family = family, lambda = grid)
      coef_without <- as.numeric(path_without$beta[v, j])
    }
    reason <- if (!is.na(coef_without) && coef_without != 0) {
      paste0("Overlaps with ", paste(partners, collapse = ", "), ": enters once they are removed")
    } else if (binary && family == "binomial" && group[["events"]] < min_events) {
      paste0("Small group (", group[["events"]], " events): too little information to separate its effect")
    } else if (r2 >= 0.5) {
      "Largely reproduced by the other predictors together"
    } else if (unadjusted[["p"]] >= 0.05) {
      "No association with the outcome, even unadjusted"
    } else {
      "Association explained by the other predictors (confounding)"
    }
    details[[v]] <- list(group = group, path = path, unadjusted = unadjusted, r2_other_predictors = r2,
                         partners = partner_table, overlap = overlap, coef_without_partners = coef_without, reason = reason)
  }
  summary <- do.call(rbind, lapply(names(details), function(v) {
    d <- details[[v]]
    data.frame(variable = v, patients = d$group[["patients"]], events = d$group[["events"]],
               lambdas_nonzero = nrow(d$path), first_entry_lambda = if (nrow(d$path)) max(d$path$lambda) else NA,
               unadjusted_estimate = d$unadjusted[["estimate"]], unadjusted_OR = d$unadjusted[["OR"]], unadjusted_p = d$unadjusted[["p"]],
               r2_other_predictors = d$r2_other_predictors,
               top_partner = if (nrow(d$partners)) d$partners$partner[1] else NA, top_partner_cor = if (nrow(d$partners)) d$partners$cor[1] else NA,
               coef_without_partners = d$coef_without_partners, reason = d$reason, row.names = NULL)
  }))
  if (PrintSummary) { print(summary, row.names = FALSE, digits = 3) }
  invisible(list(details = details, summary = summary))
}

#' Build a better lambda grid from a previous wlasso fit
#'
#' Places a dense log-spaced grid over the region that matters: below the old
#' grid when lambda.min was at its smallest lambda, above it when lambda.min was
#' at its largest, and around the 1-SE lambda and lambda.min otherwise. The
#' previous fit's lambdas outside the dense range are kept, so glmnet still
#' warm-starts from the all-zero model.
#'
#' @param fit A \code{"wlasso"} object.
#' @param n_dense Number of lambdas in the dense range.
#' @param extend Factor by which the grid is extended beyond its old end.
#' @param margin Number of previous grid positions added on each side of the
#'   region of interest.
#' @param PrintSummary Logical; print what was done.
#' @return A list: \code{grid} (decreasing; pass as \code{lambda.grid}),
#'   \code{case}, \code{dense_range} and \code{previous} (the previous fit's
#'   lambda choices).
#' @examples
#' \dontrun{
#' grid <- svyVarSel_refine_lambda_grid(fit)$grid
#' }
#' @export
svyVarSel_refine_lambda_grid <- function(fit, n_dense = 40, extend = 10, margin = 3, PrintSummary = TRUE) {
  if (!inherits(fit, "wlasso")) { stop("Please, insert an object of class 'wlasso'.") }
  sel <- svyVarSel_lambda_select(fit, lambda_rule = "min", PrintSummary = FALSE)
  grid <- sel$grid; n <- length(grid)   # decreasing, as glmnet fits it; positions match sel$j_min and sel$j_1se
  if (sel$j_min == n) {
    case <- "lambda.min at the smallest lambda: extended the grid downward"
    upper <- grid[max(1, sel$j_1se - margin)]; lower <- min(grid) / extend
  } else if (sel$j_min == 1) {
    case <- "lambda.min at the largest lambda: extended the grid upward"
    upper <- max(grid) * extend; lower <- grid[min(n, 1 + margin)]
  } else {
    case <- "lambda.min inside the grid: placed a finer grid around the 1-SE lambda and lambda.min"
    upper <- grid[max(1, sel$j_1se - margin)]; lower <- grid[min(n, sel$j_min + margin)]
  }
  dense <- exp(seq(log(upper), log(lower), length.out = n_dense))
  new_grid <- sort(unique(c(grid[grid > upper], dense, grid[grid < lower])), decreasing = TRUE)   # previous lambdas outside the dense range
  if (PrintSummary) {
    message(case, sprintf(" (dense range %.3g to %.3g; %d lambdas).", upper, lower, length(new_grid)))
  }
  list(grid = new_grid, case = case, dense_range = c(upper = upper, lower = lower), previous = sel$summary)
}

#' Fit a wlasso model, refining the lambda grid until lambda.min is inside it
#'
#' Fits \code{svyVarSel_wlasso_adjusted()} and, while lambda.min sits at either
#' end of the grid, refines the grid with \code{svyVarSel_refine_lambda_grid()}
#' and refits, up to \code{max_rounds} fits. \code{seed} is set before every fit
#' so all rounds use the same folds and only the grid changes.
#'
#' @param ... Arguments passed to \code{svyVarSel_wlasso_adjusted()}.
#' @param lambda.grid Starting grid; \code{NULL} uses glmnet's own path.
#' @param max_rounds Maximum number of fits.
#' @param seed Random seed set before every fit, or \code{NULL}.
#' @inheritParams svyVarSel_refine_lambda_grid
#' @param PrintSummary Logical; print the tuning rounds.
#' @return The final \code{"wlasso"} fit with \code{$tuning} added (one row per
#'   round). Warns when lambda.min is still at an end after the last round.
#' @examples
#' \dontrun{
#' tuned <- svyVarSel_wlasso_tune(data = d, col.y = "Death_01", col.x = preds,
#'   family = "binomial", strata = "strata_id", weights = "DISCWT",
#'   cluster = "cluster_id", method = "dCV", k = 5, R = 1, seed = 1621)
#' }
#' @export
svyVarSel_wlasso_tune <- function(..., lambda.grid = NULL, max_rounds = 3, seed = NULL,
                                  n_dense = 40, extend = 10, margin = 3, PrintSummary = TRUE) {
  if (max_rounds < 1 || max_rounds != round(max_rounds)) { stop("max_rounds must be a whole number of at least 1.") }
  tuning <- data.frame()
  for (round in seq_len(max_rounds)) {
    if (!is.null(seed)) { set.seed(seed) }
    fit <- svyVarSel_wlasso_adjusted(..., lambda.grid = lambda.grid)
    sel <- svyVarSel_lambda_select(fit, lambda_rule = "min", PrintSummary = FALSE)
    at_end <- sel$j_min %in% c(1, length(sel$grid))
    refine <- if (at_end && round < max_rounds) svyVarSel_refine_lambda_grid(fit, n_dense = n_dense, extend = extend,
                                                                              margin = margin, PrintSummary = FALSE) else NULL
    tuning <- rbind(tuning, data.frame(round = round, n_lambdas = length(sel$grid), lambda_min = sel$grid[sel$j_min],
                                       position = sel$j_min, at_grid_end = at_end,
                                       next_step = if (!at_end) "stop: lambda.min inside the grid" else
                                         if (is.null(refine)) "stop: max_rounds reached" else refine$case))
    if (PrintSummary) {
      message(sprintf("Round %d: %d lambdas, lambda.min = %.3g at position %d%s", round, length(sel$grid),
                      sel$grid[sel$j_min], sel$j_min, if (at_end) " (end of grid)" else " (inside the grid)"))
    }
    if (is.null(refine)) { break }
    lambda.grid <- refine$grid
  }
  if (at_end) { warning(sprintf("lambda.min is still at the end of the grid after %d round%s; increase max_rounds or extend.",
                                round, if (round == 1) "" else "s")) }
  fit$tuning <- tuning
  if (PrintSummary) { print(tuning, row.names = FALSE, digits = 3) }
  fit
}

#' Stability of a wlasso result across folds (K) and repeats (R)
#'
#' Refits \code{svyVarSel_wlasso_adjusted()} for every K x R x seed with the same
#' lambda grid and reports, for each K x R setting, how much lambda.min and the
#' selected predictors change across seeds. K and R set how the
#' cross-validation is run, so they are not chosen by lowest error. The
#' suggestion is the fastest setting that is stable.
#'
#' @param ... Arguments passed to \code{svyVarSel_wlasso_adjusted()} (not
#'   \code{k}, \code{R} or \code{lambda.grid}).
#' @param K Numbers of folds to try.
#' @param R Numbers of repeats to try.
#' @param seeds At least two random seeds.
#' @param lambda.grid Lambda grid used for every fit (for example a tuned grid).
#' @param lambda_rule \code{"min"} or \code{"1se"}.
#' @param agreement_threshold Minimum share of predictors selected by every seed
#'   among those selected by any seed.
#' @param spread_threshold Maximum \code{log10(max / min)} of lambda.min across
#'   seeds.
#' @param PrintSummary Logical; print progress and the summary.
#' @return Invisibly, a list: \code{fits} (one row per fit), \code{summary} (one
#'   row per K x R), \code{selection_frequency} and \code{suggested} (K, R and
#'   the reason).
#' @examples
#' \dontrun{
#' kr <- svyVarSel_kr_stability(data = d, col.y = "Death_01", col.x = preds,
#'   family = "binomial", strata = "strata_id", weights = "DISCWT",
#'   cluster = "cluster_id", method = "dCV", K = c(5, 10), R = c(1, 3))
#' kr$suggested
#' }
#' @export
svyVarSel_kr_stability <- function(..., K = c(5, 10), R = c(1, 3, 5), seeds = c(1621, 1622), lambda.grid = NULL,
                                   lambda_rule = c("min", "1se"), agreement_threshold = 0.9, spread_threshold = 0.1,
                                   PrintSummary = TRUE) {
  lambda_rule <- match.arg(lambda_rule)
  if (length(seeds) < 2) { stop("Use at least two seeds: stability is measured across seeds.") }
  settings <- expand.grid(K = K, R = R, seed = seeds)
  fits <- data.frame(); selected <- list()
  for (i in seq_len(nrow(settings))) {
    s <- settings[i, ]
    t0 <- Sys.time()
    set.seed(s$seed)
    fit <- svyVarSel_wlasso_adjusted(..., k = s$K, R = s$R, lambda.grid = lambda.grid)
    seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    sel <- svyVarSel_lambda_select(fit, lambda_rule = lambda_rule, PrintSummary = FALSE)
    vars <- svyVarSel_selected_coefs(fit, j = sel$j)
    selected[[i]] <- rownames(vars$selected)
    fits <- rbind(fits, data.frame(K = s$K, R = s$R, seed = s$seed, lambda = sel$lambda, position = sel$j,
                                   n_vars = nrow(vars$selected), at_grid_end = sel$j_min %in% c(1, length(sel$grid)),
                                   seconds = seconds))
    if (PrintSummary) { message(sprintf("Fit %d of %d: K = %d, R = %d, seed = %d (%.0f s)", i, nrow(settings), s$K, s$R, s$seed, seconds)) }
  }
  predictors <- rownames(fit$model$grid$beta)
  key <- paste(fits$K, fits$R)
  summary <- do.call(rbind, lapply(unique(key), function(k) {
    rows <- which(key == k); sets <- selected[rows]
    union_set <- unique(unlist(sets)); common_set <- Reduce(intersect, sets)
    data.frame(K = fits$K[rows[1]], R = fits$R[rows[1]], seeds = length(rows),
               lambda_median = median(fits$lambda[rows]), lambda_min_spread = log10(max(fits$lambda[rows]) / min(fits$lambda[rows])),
               n_vars_low = min(fits$n_vars[rows]), n_vars_high = max(fits$n_vars[rows]),
               selection_agreement = if (length(union_set)) length(common_set) / length(union_set) else 1,
               any_at_grid_end = any(fits$at_grid_end[rows]), seconds = mean(fits$seconds[rows]))
  }))
  summary$stable <- summary$selection_agreement >= agreement_threshold & summary$lambda_min_spread <= spread_threshold
  selection_frequency <- sapply(unique(key), function(k) {
    rows <- which(key == k); vapply(predictors, function(p) mean(vapply(selected[rows], function(x) p %in% x, logical(1))), numeric(1))
  })
  colnames(selection_frequency) <- paste0("K", summary$K, "_R", summary$R)
  if (any(summary$stable)) {
    best <- summary[summary$stable, ][which.min(summary$seconds[summary$stable]), ]
    reason <- sprintf("fastest setting with selection agreement >= %.2f and lambda.min spread <= %.2f", agreement_threshold, spread_threshold)
  } else {
    best <- summary[order(-summary$selection_agreement, summary$lambda_min_spread), ][1, ]
    reason <- "no setting met the thresholds; this is the most stable one tried (consider larger R or more seeds)"
    warning("No K x R setting met the stability thresholds; suggesting the most stable one tried.")
  }
  suggested <- list(K = best$K, R = best$R, reason = reason)
  if (PrintSummary) {
    print(summary, row.names = FALSE, digits = 3)
    message(sprintf("Suggested: K = %d, R = %d (%s).", best$K, best$R, reason))
  }
  invisible(list(fits = fits, summary = summary, selection_frequency = selection_frequency, suggested = suggested))
}

#' Coefficient stability bootstrap for a fitted wlasso model
#'
#' Resamples clusters (hospitals) with replacement within strata and refits the
#' same penalized LASSO at a fixed lambda. Each refit runs a warm-start path down
#' to that lambda. Refits that warn, do not reach the lambda or give nonfinite
#' coefficients are recorded as failures. The percentile ranges describe
#' stability; they are not confidence intervals, and they include zero when a
#' predictor drops out of some refits.
#'
#' @param fit A \code{"wlasso"} object.
#' @param data The training data used for the fit.
#' @param outcome Name of the outcome column.
#' @param weights Name of the weight column.
#' @param cluster Name of the cluster column.
#' @param strata Name of the strata column.
#' @param lambda Lambda to refit at; defaults to the fit's lambda.min.
#' @param B Number of bootstrap refits.
#' @param seed Random seed, or \code{NULL}.
#' @param family \code{"binomial"}, \code{"poisson"} or \code{"gaussian"}.
#' @param probs Percentiles of the stability range.
#' @param min_success Minimum share of successful refits; stops below it.
#' @param progress_every Print progress every this many refits (0 = never).
#' @param PrintSummary Logical; print the results.
#' @param singleton \code{"fixed"} holds single-cluster strata fixed (with a
#'   warning); \code{"fail"} stops.
#' @return Invisibly, a list: \code{table} (selected predictors: LASSO
#'   coefficient, stability limits, percentage of refits keeping it),
#'   \code{coefficients}, \code{failures}, \code{successful}, \code{lambda} and
#'   \code{path}.
#' @examples
#' \dontrun{
#' boot <- svyVarSel_coef_bootstrap(fit, d, outcome = "Death_01",
#'   weights = "DISCWT", cluster = "cluster_id", strata = "strata_id", B = 200)
#' boot$table
#' }
#' @export
svyVarSel_coef_bootstrap <- function(fit, data, outcome, weights, cluster, strata,
                                     lambda = fit$lambda$min, B = 200, seed = 1621, family = "binomial",
                                     probs = c(0.025, 0.975), min_success = 0.9, progress_every = 10, PrintSummary = TRUE,
                                     singleton = c("fixed", "fail")) {
  singleton <- match.arg(singleton)
  if (!inherits(fit, "wlasso")) { stop("Please, insert an object of class 'wlasso'.") }
  if (!is.null(seed)) { set.seed(seed) }
  data <- as.data.frame(data)
  beta_path <- fit$model$grid$beta
  j_fit <- match(lambda, fit$lambda$grid)
  if (is.na(j_fit)) { stop("lambda must be one of the lambdas in fit$lambda$grid.") }
  final_beta <- as.matrix(beta_path[, j_fit, drop = FALSE])
  selected <- rownames(final_beta)[final_beta[, 1] != 0]
  if (!length(selected)) { stop("No predictors are selected at this lambda.") }
  x_boot <- as.matrix(data[, rownames(final_beta), drop = FALSE])
  y_boot <- data[[outcome]]
  hospitals <- unique(data[, c(cluster, strata)])
  if (anyDuplicated(hospitals[[cluster]])) { stop("A hospital (cluster) appears in more than one stratum.") }
  hospital_groups <- split(seq_len(nrow(hospitals)), hospitals[[strata]])
  hospital_row <- match(data[[cluster]], hospitals[[cluster]])
  if (all(lengths(hospital_groups) == 1)) { stop("No strata can be resampled.") }
  if (any(lengths(hospital_groups) == 1) && singleton == "fail") stop("Singleton training strata: bootstrap stopped; choose singleton='fixed' explicitly to condition on these hospitals.")
  if (any(lengths(hospital_groups) == 1)) { warning("Singleton hospitals are held fixed; their between-hospital variation is not represented.") }
  boot_path <- sort(unique(c(exp(seq(log(max(1, lambda * 10)), log(lambda), length.out = 60)), lambda)), decreasing = TRUE)
  boot_coefficients <- matrix(NA_real_, B, length(selected), dimnames = list(NULL, selected))
  boot_failures <- rep(NA_character_, B)
  for (b in seq_len(B)) {
    sampled <- unlist(lapply(hospital_groups, function(h) h[sample.int(length(h), length(h), replace = TRUE)]), use.names = FALSE)
    boot_weights <- data[[weights]] * tabulate(sampled, nbins = nrow(hospitals))[hospital_row]
    tryCatch(withCallingHandlers({
      refit <- glmnet::glmnet(x_boot, y_boot, weights = boot_weights, family = family,
                              alpha = 1, standardize = TRUE, intercept = TRUE, lambda = boot_path,
                              control = list(maxit = 1000000))
      j_boot <- match(lambda, refit$lambda)
      if (refit$jerr != 0 || is.na(j_boot)) { stop("Refit did not reach the selected lambda.") }
      values <- as.numeric(refit$beta[selected, j_boot, drop = FALSE])
      if (any(!is.finite(values))) { stop("Nonfinite coefficients.") }
      boot_coefficients[b, ] <- values
    }, warning = function(e) stop(conditionMessage(e), call. = FALSE)),
    error = function(e) { boot_failures[b] <<- conditionMessage(e) })
    if (progress_every > 0 && (b %% progress_every == 0 || b == B)) { message(b, " of ", B, " complete") }
  }
  successful <- is.na(boot_failures)
  if (PrintSummary) {
    cat("Successful bootstrap refits:", sum(successful), "of", B, "\n")
    if (any(!successful)) { print(table(boot_failures[!successful])) }
  }
  if (sum(successful) < min_success * B) { stop("Too many failed refits to summarize reliably.") }
  stability_limits <- vapply(seq_along(selected), function(j) quantile(boot_coefficients[successful, j], probs), numeric(2))
  table <- data.frame(variable = selected, lasso_coefficient = final_beta[selected, 1],
                      stability_lower = stability_limits[1, ], stability_upper = stability_limits[2, ],
                      selected_percent = 100 * colMeans(boot_coefficients[successful, , drop = FALSE] != 0),
                      row.names = NULL)
  if (PrintSummary) { print(table, digits = 4) }
  invisible(list(table = table, coefficients = boot_coefficients, failures = boot_failures,
                 successful = sum(successful), lambda = lambda, path = boot_path))
}

#' Unpenalized refit of the LASSO-selected predictors
#'
#' Refits the predictors selected at \code{lambda} (or \code{variables}) with
#' \code{survey::svyglm}, \code{stats::glm} or \code{stats::lm} and reports
#' effects with confidence intervals, plus how much the penalty shrank each
#' coefficient. The intervals ignore that the predictors were chosen from the
#' same data. Separation (fitted probabilities of 0 or 1, or a coefficient
#' above 10 in absolute value) is detected and returned, because it makes the
#' refit's effects unreliable.
#'
#' @param fit A \code{"wlasso"} object (optional when \code{variables} is
#'   given).
#' @param data The training data.
#' @param outcome Name of the outcome column.
#' @param weights Name of the weight column. Required for \code{"svyglm"};
#'   optional for \code{"glm"} and \code{"lm"}.
#' @param cluster Name of the cluster column. Required for the usual clustered
#'   \code{"svyglm"} fit; may be \code{NULL} for an unclustered survey design.
#' @param strata Name of the strata column, or \code{NULL}.
#' @param lambda Lambda whose selection is refitted; defaults to lambda.min.
#' @param variables Predictors to refit instead of the LASSO selection.
#' @param add_variables Predictors added to the selection.
#' @param level Confidence level.
#' @param lonely_psu Value of \code{options(survey.lonely.psu)} for this fit.
#' @param family \code{"binomial"} (odds ratios), \code{"poisson"} (rate ratios)
#'   or \code{"gaussian"} (differences).
#' @param method \code{"svyglm"}, \code{"glm"} or \code{"lm"}. The
#'   \code{"lm"} method requires \code{family = "gaussian"}.
#' @param PrintSummary Logical; print the table and comparison.
#' @return Invisibly, a list: \code{table} (effect, limits, p-value, standard
#'   error, LASSO and unpenalized coefficients, shrinkage ratio, shift in
#'   standard errors, sign change), \code{comparison} (linear-predictor slope and
#'   correlation, median shrinkage, largest shift), \code{separation}
#'   (\code{detected}, \code{records_at_bounds}, \code{variables}), \code{model}
#'   and \code{formula}.
#' @examples
#' \dontrun{
#' refit <- svyVarSel_unpenalized_refit(fit, d, outcome = "Death_01",
#'   weights = "DISCWT", cluster = "cluster_id", strata = "strata_id",
#'   method = "svyglm")
#' refit$table
#' }
#' @export
svyVarSel_unpenalized_refit <- function(fit, data, outcome, weights = NULL, cluster = NULL, strata = NULL,
                                        lambda = fit$lambda$min, variables = NULL, add_variables = NULL, level = 0.95,
                                        lonely_psu = "adjust", family = c("binomial", "poisson", "gaussian"),
                                        method = c("svyglm", "glm", "lm"), PrintSummary = TRUE) {
  family <- match.arg(family)
  method <- match.arg(method)
  if (method == "lm" && family != "gaussian")
    stop("method = 'lm' requires family = 'gaussian'.")
  has_fit <- !missing(fit) && inherits(fit, "wlasso")
  if (has_fit) {
    j <- match(lambda, fit$lambda$grid)
    if (is.na(j)) { stop("lambda must be one of the lambdas in fit$lambda$grid.") }
    beta <- as.matrix(fit$model$grid$beta[, j, drop = FALSE])
  }
  if (is.null(variables)) {
    if (!has_fit) { stop("Please, insert an object of class 'wlasso', or pass variables.") }
    variables <- rownames(beta)[beta[, 1] != 0]
  }
  variables <- unique(c(variables, add_variables))
  model_formula <- if (length(variables)) { reformulate(variables, response = outcome) } else { as.formula(paste(outcome, "~ 1")) }
  data <- as.data.frame(data)
  needed <- c(outcome, variables, weights, if (method == "svyglm") c(cluster, strata))
  absent <- setdiff(needed, names(data))
  if (length(absent)) stop("Columns not found for the refit: ", paste(absent, collapse = ", "))
  refit_weights <- if (is.null(weights)) rep(1, nrow(data)) else data[[weights]]
  if (!is.numeric(refit_weights) || any(!is.finite(refit_weights) | refit_weights <= 0))
    stop("Refit weights must be finite and positive.")
  if (method == "svyglm" && is.null(weights))
    stop("method = 'svyglm' requires a weights column.")
  scaled_weights <- refit_weights / mean(refit_weights)
  refit <- local({
    if (method == "svyglm") {
      old_options <- options(survey.lonely.psu = lonely_psu)
      on.exit(options(old_options))
      id_formula <- if (is.null(cluster)) ~1 else as.formula(paste0("~", cluster))
      strata_formula <- if (is.null(strata)) NULL else as.formula(paste0("~", strata))
      design <- survey::svydesign(ids = id_formula, strata = strata_formula,
                                  weights = as.formula(paste0("~", weights)), nest = TRUE, data = data)
      model <- survey::svyglm(model_formula, design = design, na.action = na.fail,
                              family = switch(family, binomial = quasibinomial(), poisson = quasipoisson(), gaussian = gaussian()))
    } else {
      model <- if (method == "glm")
        stats::glm(model_formula, data = data, weights = scaled_weights, na.action = na.fail,
                   family = switch(family, binomial = stats::binomial(), poisson = stats::poisson(), gaussian = stats::gaussian()))
      else
        stats::lm(model_formula, data = data, weights = scaled_weights, na.action = na.fail)
    }
    if (!is.null(model$converged) && !isTRUE(model$converged))
      stop("Confidence-interval refit did not converge.")
    if (any(!is.finite(stats::coef(model))) || any(!is.finite(stats::vcov(model))))
      stop("Confidence-interval refit has non-finite or aliased coefficients or covariance.")
    coefs <- summary(model)$coefficients
    limits <- if (method == "glm") stats::confint.default(model, level = level) else stats::confint(model, level = level)
    list(model = model, limits = limits, se = coefs[, 2], p = coefs[, 4])
  })
  model <- refit$model
  estimates <- coef(model)
  terms <- names(estimates)[names(estimates) != "(Intercept)"]
  pct <- format(100 * level)
  scale_f <- if (family == "gaussian") identity else exp
  table <- data.frame(variable = terms, estimate = scale_f(estimates[terms]),
                      lower = scale_f(refit$limits[terms, 1]), upper = scale_f(refit$limits[terms, 2]),
                      p_value = refit$p[terms], std_error = unname(refit$se[terms]), row.names = NULL)
  names(table)[2:4] <- c(switch(family, binomial = "odds_ratio", poisson = "rate_ratio", gaussian = "difference"),
                         paste0(c("lower_", "upper_"), pct))
  lasso_coef <- if (has_fit) { unname(ifelse(terms %in% rownames(beta), beta[match(terms, rownames(beta)), 1], NA_real_)) } else { rep(NA_real_, length(terms)) }
  unpenalized_coef <- unname(estimates[terms])
  table$lasso_coef <- lasso_coef
  table$unpenalized_coef <- unpenalized_coef
  table$shrinkage_ratio <- ifelse(unpenalized_coef != 0, lasso_coef / unpenalized_coef, NA_real_)
  table$shift_in_SE <- sign(unpenalized_coef) * (unpenalized_coef - lasso_coef) / unname(refit$se[terms])   # > 0 = toward zero
  table$sign_change <- !is.na(lasso_coef) & lasso_coef != 0 & sign(lasso_coef) != sign(unpenalized_coef)
  at_bounds <- if (family == "binomial") sum(model$fitted.values < 1e-8 | model$fitted.values > 1 - 1e-8) else 0L
  extreme <- if (family == "gaussian") character(0) else terms[abs(unpenalized_coef) > 10]
  separation <- list(detected = at_bounds > 0 || length(extreme) > 0, records_at_bounds = at_bounds, variables = extreme)
  comparison <- NULL
  if (has_fit) {
    w <- refit_weights
    lasso_lp <- as.numeric(fit$model$grid$a0[j]) + drop(as.matrix(data[, rownames(beta), drop = FALSE]) %*% beta[, 1])
    refit_lp <- if (method == "lm") as.numeric(stats::fitted(model)) else as.numeric(model$linear.predictors)
    kept <- !is.na(lasso_coef) & lasso_coef != 0
    comparison <- list(lp_slope = unname(coef(lm(refit_lp ~ lasso_lp, weights = w))[2]),
                       lp_correlation = stats::cov.wt(cbind(lasso_lp, refit_lp), wt = w, cor = TRUE)$cor[1, 2],
                       median_shrinkage = if (any(kept)) median(table$shrinkage_ratio[kept], na.rm = TRUE) else NA_real_,
                       max_abs_shift_in_SE = if (length(terms)) max(abs(table$shift_in_SE), na.rm = TRUE) else NA_real_)
  }
  if (PrintSummary) {
    print(table, digits = 4)
    if (!is.null(comparison)) {
      message(sprintf("LASSO vs unpenalized: linear-predictor slope %.3f, correlation %.4f, median shrinkage %.3f, largest shift %.2f SE.",
                      comparison$lp_slope, comparison$lp_correlation, comparison$median_shrinkage, comparison$max_abs_shift_in_SE))
    }
    if (separation$detected) {
      message(sprintf("Separation in the refit (%d records with fitted probability 0 or 1; extreme coefficients: %s): its effects are unreliable.",
                      separation$records_at_bounds, if (length(extreme)) paste(extreme, collapse = ", ") else "none"))
    }
  }
  invisible(list(table = table, comparison = comparison, separation = separation, model = model,
                 formula = model_formula, method = method))
}

#' Hold out whole clusters within strata for validation
#'
#' Splits the data by cluster within each stratum, so no cluster is in both sets
#' and both keep the survey design. In a stratum with n clusters, round(share x
#' n) are held out only if both sides keep at least \code{min_clusters}. The
#' function is model-agnostic and can be used with any data set in which each
#' cluster belongs to exactly one stratum. Run it after data preparation and
#' before fitting, tuning, or exploratory model selection.
#'
#' @param data Data frame.
#' @param cluster Name of the cluster column.
#' @param strata Name of the strata column.
#' @param share Target share of clusters to hold out.
#' @param min_clusters Minimum clusters each side must keep in a stratum.
#' @param outcome Optional numeric outcome column; its nonmissing values are
#'   summed in the summary.
#' @param seed Random seed, or \code{NULL}.
#' @param PrintSummary Logical; print the summary.
#' @return A list: \code{training}, \code{validation}, \code{held_out} (cluster
#'   IDs), \code{summary} (records, clusters, strata and outcome total per set)
#'   and \code{achieved_share}.
#' @examples
#' \dontrun{
#' split <- ClusterStratifiedHoldoutSplit(d, cluster = "site_id",
#'   strata = "region", share = 0.2, outcome = "response")
#' }
#' @export
ClusterStratifiedHoldoutSplit <- function(data, cluster, strata, share = 0.2,
                                          min_clusters = 2, outcome = NULL,
                                          seed = 1621, PrintSummary = TRUE) {
  if (!is.character(cluster) || length(cluster) != 1L || is.na(cluster) || !nzchar(cluster) ||
      !is.character(strata) || length(strata) != 1L || is.na(strata) || !nzchar(strata))
    stop("cluster and strata must each name one column in data.")
  if (identical(cluster, strata)) stop("cluster and strata must name different columns.")
  if (!is.numeric(share) || length(share) != 1L || !is.finite(share) || share <= 0 || share >= 1)
    stop("share must be one finite number between 0 and 1.")
  if (!is.numeric(min_clusters) || length(min_clusters) != 1L || !is.finite(min_clusters) ||
      min_clusters < 1 || min_clusters != floor(min_clusters))
    stop("min_clusters must be one positive integer.")
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)))
    stop("seed must be one finite number or NULL.")
  if (!is.logical(PrintSummary) || length(PrintSummary) != 1L || is.na(PrintSummary))
    stop("PrintSummary must be TRUE or FALSE.")

  df <- as.data.frame(data)
  if (!is.null(outcome) &&
      (!is.character(outcome) || length(outcome) != 1L || is.na(outcome) || !nzchar(outcome)))
    stop("outcome must name one numeric column or be NULL.")
  absent <- setdiff(c(cluster, strata, outcome), names(df))
  if (length(absent)) stop("Columns not found: ", paste(absent, collapse = ", "))
  if (!is.null(outcome) && !is.numeric(df[[outcome]]))
    stop("outcome must name one numeric column or be NULL.")
  if (anyNA(df[[cluster]]) || anyNA(df[[strata]]))
    stop("cluster and strata cannot contain missing values.")

  if (!is.null(seed)) { set.seed(seed) }
  clusters <- unique(df[, c(cluster, strata), drop = FALSE])
  if (anyDuplicated(clusters[[cluster]]))
    stop("A cluster appears in more than one stratum.")
  by_stratum <- split(seq_len(nrow(clusters)), as.character(clusters[[strata]]))
  held_rows <- unlist(lapply(by_stratum, function(rows) {
    n_hold <- round(share * length(rows))
    if (n_hold < min_clusters || length(rows) - n_hold < min_clusters) return(integer(0))
    rows[sample.int(length(rows), n_hold)]
  }), use.names = FALSE)
  if (!length(held_rows))
    stop("No stratum has enough clusters to hold any out; lower min_clusters or change share.")
  held_out <- clusters[[cluster]][held_rows]
  in_validation <- as.character(df[[cluster]]) %in% as.character(held_out)
  describe <- function(rows) c(records = sum(rows), clusters = length(unique(df[[cluster]][rows])),
                               strata = length(unique(df[[strata]][rows])),
                               outcome_total = if (is.null(outcome)) NA_real_ else sum(df[[outcome]][rows], na.rm = TRUE))
  summary <- rbind(training = describe(!in_validation), validation = describe(in_validation))
  achieved <- length(held_out) / nrow(clusters)
  if (PrintSummary) {
    print(summary)
    message(sprintf("Held out %d of %d clusters (%.1f%%; target %.0f%%) from %d of %d strata.", length(held_out), nrow(clusters),
                    100 * achieved, 100 * share, summary["validation", "strata"], length(by_stratum)))
  }
  list(training = data[!in_validation, ], validation = data[in_validation, ], held_out = held_out,
       summary = summary, achieved_share = achieved)
}

#' Deprecated name for cluster-stratified holdout splitting
#'
#' Use \code{ClusterStratifiedHoldoutSplit()} for model-agnostic cluster-level
#' validation splitting. This wrapper is retained for existing scripts.
#'
#' @inheritParams ClusterStratifiedHoldoutSplit
#' @param min_hospitals Deprecated name for \code{min_clusters}.
#' @return The same result as \code{ClusterStratifiedHoldoutSplit()}.
#' @export
svyVarSel_holdout_split <- function(data, cluster, strata, share = 0.2,
                                    min_hospitals = 2, outcome = NULL,
                                    seed = 1621, PrintSummary = TRUE) {
  .Deprecated("ClusterStratifiedHoldoutSplit")
  ClusterStratifiedHoldoutSplit(
    data = data, cluster = cluster, strata = strata, share = share,
    min_clusters = min_hospitals, outcome = outcome, seed = seed,
    PrintSummary = PrintSummary
  )
}

#' Evaluate a fixed binomial wlasso model on held-out records
#'
#' Computes AUC, Brier score, log loss and observed/expected events for a
#' fitted model on validation records that were not used for fitting or tuning,
#' with stratified cluster-bootstrap percentile intervals. The intervals are
#' conditional on the fitted model; they do not include training or tuning
#' uncertainty. Observed/expected checks the overall calibration level, not the
#' calibration slope.
#'
#' @param model A \code{"wlasso"} object; its \code{$model$min} is evaluated.
#' @param validation_data Validation records with the same coding as the
#'   training data.
#' @param training_clusters Cluster IDs of every record used for training or
#'   tuning; validation must not share any.
#' @param outcome Name of the 0/1 outcome column.
#' @param weights Name of the weight column.
#' @param cluster Name of the cluster column.
#' @param strata Name of the strata column.
#' @param B Number of bootstrap draws (at least 20).
#' @param conf Confidence level of the intervals.
#' @param seed Random seed; the caller's random-number stream is restored.
#' @param singleton \code{"fail"} stops on single-cluster strata;
#'   \code{"fixed"} holds them fixed.
#' @return A list: \code{summary} (metric, estimate, lower, upper,
#'   valid_bootstraps), \code{calibration} (observed and predicted event rates
#'   by predicted-risk tenth), \code{predictions}, \code{bootstrap} and settings.
#' @examples
#' \dontrun{
#' ev <- evaluate_prediction_model(fit, validation, unique(training$cluster_id),
#'   outcome = "Death_01", weights = "DISCWT", cluster = "cluster_id",
#'   strata = "strata_id", B = 500)
#' ev$summary
#' }
#' @export
evaluate_prediction_model <- function(model, validation_data, training_clusters,
                                      outcome, weights,
                                      cluster, strata,
                                      B = 500L, conf = 0.95, seed = 1621L,
                                      singleton = c("fail", "fixed")) {
  singleton <- match.arg(singleton)
  if (length(B) != 1L || !is.finite(B) || B < 20 || B != floor(B) ||
      length(conf) != 1L || !is.finite(conf) || conf <= 0 || conf >= 1)
    stop("Use an integer B >= 20 and 0 < conf < 1.")
  if (missing(training_clusters) || !length(training_clusters) || anyNA(training_clusters))
    stop("Supply hospital-year IDs from all records used for training/tuning.")
  beta <- as.matrix(model$model$min$beta)
  intercept <- as.numeric(model$model$min$a0)
  if (ncol(beta) != 1L || is.null(rownames(beta)) || anyDuplicated(rownames(beta)) ||
      length(intercept) != 1L || any(!is.finite(c(beta, intercept))))
    stop("Expected a valid selected binomial wlasso model with named coefficients.")
  predictors <- rownames(beta)
  d <- as.data.frame(validation_data)
  required <- unique(c(predictors, outcome, weights, cluster, strata))
  if (!all(required %in% names(d))) stop("Missing validation columns: ",
                                         paste(setdiff(required, names(d)), collapse = ", "))
  if (!nrow(d) || anyNA(d[, required, drop = FALSE]))
    stop("Validation data must be nonempty and complete for required columns.")
  if (any(as.character(d[[cluster]]) %in% as.character(training_clusters)))
    stop("Training/tuning and validation hospital-year IDs overlap. Do not evaluate on training data.")
  numeric_cols <- unique(c(predictors, outcome, weights))
  if (!all(vapply(d[, numeric_cols, drop = FALSE], is.numeric, logical(1))))
    stop("Predictors, outcome and weights must be numeric, using the training coding.")
  x <- as.matrix(d[, predictors, drop = FALSE])
  y <- d[[outcome]]; w <- d[[weights]]
  if (any(!is.finite(x)) || any(!is.finite(w)) || any(w <= 0) ||
      !all(y %in% c(0, 1)) || length(unique(y)) != 2L)
    stop("Need finite predictors, positive finite weights, and both 0/1 outcome classes.")
  eta <- drop(x %*% beta) + intercept
  if (any(!is.finite(eta))) stop("Nonfinite predictions; check model coefficients and input scaling.")
  p <- plogis(eta)
  # Compute predictions, losses and AUC ordering only once.
  log_loss <- pmax(eta, 0) - y * eta + log1p(exp(-abs(eta)))
  brier <- (y - p)^2
  ord <- order(eta)
  tie_ends <- c(which(diff(eta[ord]) != 0), length(y))
  metrics <- function(wb) {
    pos_cum <- cumsum((wb * y)[ord])[tie_ends]
    neg_cum <- cumsum((wb * (1 - y))[ord])[tie_ends]
    pos <- diff(c(0, pos_cum)); neg <- diff(c(0, neg_cum))
    denom <- sum(wb * y) * sum(wb * (1 - y))
    auc <- if (denom > 0) sum(pos * (neg_cum - neg / 2)) / denom else NA_real_
    expected <- sum(wb * p)
    c(AUC = auc, Brier = sum(wb * brier) / sum(wb),
      Log_loss = sum(wb * log_loss) / sum(wb),
      Observed_expected = if (expected > 0) sum(wb * y) / expected else NA_real_)
  }
  # Resample entire hospitals within strata; retain original discharge weights.
  hospitals <- unique(data.frame(id = as.character(d[[cluster]]),
                                 stratum = as.character(d[[strata]])))
  if (anyDuplicated(hospitals$id)) stop("A hospital-year ID occurs in multiple strata.")
  row_hospital <- match(as.character(d[[cluster]]), hospitals$id)
  groups <- split(seq_len(nrow(hospitals)), hospitals$stratum)
  singletons <- names(groups)[lengths(groups) == 1L]
  if (length(singletons) && singleton == "fail") stop(
    "Validation contains ", length(singletons), " singleton strata. Review their design; ",
    "singleton='fixed' is an explicit sensitivity option, not a variance correction.")
  if (all(lengths(groups) == 1L)) stop("No stratum has multiple hospitals; cluster bootstrap cannot estimate uncertainty.")
  if (length(singletons)) warning("Singleton hospitals held fixed: their between-hospital variance is not estimated.")
  # Keep the caller's random-number stream unchanged.
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)
  point <- metrics(w)
  boot <- matrix(NA_real_, B, length(point), dimnames = list(NULL, names(point)))
  for (b in seq_len(B)) {
    draw <- unlist(lapply(groups, function(ids) ids[sample.int(length(ids), length(ids), replace = TRUE)]),
                   use.names = FALSE)
    multiplicity <- tabulate(draw, nbins = nrow(hospitals))
    boot[b, ] <- metrics(w * multiplicity[row_hospital])
  }
  valid <- colSums(is.finite(boot))
  interval <- vapply(seq_along(point), function(j) {
    if (valid[j] < ceiling(0.9 * B)) return(c(NA_real_, NA_real_))
    unname(quantile(boot[is.finite(boot[, j]), j], c((1 - conf) / 2, (1 + conf) / 2)))
  }, numeric(2))
  if (any(valid < B)) warning("Some bootstrap metrics were undefined; inspect valid_bootstraps. CIs withheld below 90% valid.")
  results <- data.frame(metric = names(point), estimate = unname(point),
                        lower = interval[1, ], upper = interval[2, ], valid_bootstraps = unname(valid), row.names = NULL)
  # Descriptive calibration by predicted-risk bands; no recalibration/refitting.
  cuts <- unique(as.numeric(quantile(p, seq(0, 1, 0.1))))
  bands <- if (length(cuts) > 1L) cut(p, cuts, include.lowest = TRUE) else factor(rep("All", length(p)))
  calibration <- do.call(rbind, lapply(split(seq_along(y), bands, drop = TRUE), function(i) {
    data.frame(records = length(i), predicted = weighted.mean(p[i], w[i]),
               observed = weighted.mean(y[i], w[i]))
  }))
  list(summary = results, calibration = calibration, predictions = p,
       bootstrap = boot, confidence = conf, B = B, seed = seed,
       validation_records = nrow(d), validation_hospitals = nrow(hospitals),
       singleton_strata = singletons, singleton_handling = singleton,
       note = "Conditional performance intervals for a fixed model; no training/tuning uncertainty. O/E ideally equals 1; it does not measure calibration slope.")
}
