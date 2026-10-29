#' Add covariate flags to setup rows
#'
#' @param setup Setup table.
#' @param flags Optional \code{MatchingAssumptionCheck()$table}.
#' @return \code{setup} with \code{flagged_covariates}, \code{covariate_flags}
#'   and \code{has_covariate_flags}.
#' @keywords internal
#' @noRd
PSRegressionFlags <- function(setup, flags = NULL) {
  if (is.null(flags)) {
    if (!"flagged_covariates" %in% names(setup)) setup$flagged_covariates <- ""
    if (!"covariate_flags" %in% names(setup)) setup$covariate_flags <- ""
    if (!"has_covariate_flags" %in% names(setup)) setup$has_covariate_flags <- FALSE
    return(setup)
  }
  setup$flagged_covariates <- ""
  setup$covariate_flags <- ""
  setup$has_covariate_flags <- FALSE
  if (is.null(flags) || !nrow(flags)) return(setup)
  for (i in seq_len(nrow(setup))) {
    vars <- PS_BaseVars(setup$MatchingVariable[i])
    f <- flags[flags$variable %in% vars, , drop = FALSE]
    if ("exposure" %in% names(f)) f <- f[f$exposure == (if ("ExposureVariable" %in% names(setup)) setup$ExposureVariable[i] else setup$IndependentVariable[i]), , drop = FALSE]
    f <- f[f$flagged %in% TRUE, , drop = FALSE]
    setup$has_covariate_flags[i] <- nrow(f) > 0
    setup$flagged_covariates[i] <- paste(unique(f$variable), collapse = "; ")
    # With checks for several estimands, each flag names its estimand, e.g. "AGE [ATE]: review: ..."
    label <- if ("estimand" %in% names(f)) paste0(f$variable, " [", f$estimand, "]") else f$variable
    setup$covariate_flags[i] <- paste(unique(paste(label, f$action, sep = ": ")), collapse = "; ")
  }
  setup
}

#' Resolve each dependent's outcome model type
#'
#' @param DT Data frame.
#' @param dependents Dependent columns.
#' @param type \code{NULL} (0/1 outcomes are Logistic; all others are Linear),
#'   one type for all, or a named vector by dependent. Poisson must be explicit.
#' @return A character vector of \code{"Logistic"}, \code{"Linear"} or
#'   \code{"Poisson"}.
#' @keywords internal
#' @noRd
PSResolveType <- function(DT, dependents, type = NULL) {
  allowed <- c("Logistic", "Linear", "Poisson")
  if (!is.null(type)) {
    if (!is.character(type) || !length(type) || anyNA(type) || any(!type %in% allowed))
      stop("type must contain Logistic, Linear or Poisson.")
    if (is.null(names(type)) && length(type) != 1L)
      stop("Use one type for all outcomes or a named vector, e.g. c(LOS = 'Linear').")
    if (!is.null(names(type)) && (anyNA(names(type)) || any(!nzchar(trimws(names(type)))) ||
        anyDuplicated(names(type)) || any(!names(type) %in% dependents)))
      stop("type names must be unique outcome names present in dependents.")
  }
  vapply(dependents, function(dep) {
    if (!is.null(type) && is.null(names(type))) return(type[[1]])
    if (!is.null(type) && dep %in% names(type)) return(type[[dep]])
    c(binomial = "Logistic", poisson = "Poisson", gaussian = "Linear")[[PS_family(DT[[dep]])]]
  }, character(1), USE.NAMES = FALSE)
}

#' Parse a covariate specification safely
#'
#' Only formula operators, \code{I()}, \code{log()}, \code{sqrt()},
#' \code{factor()}, \code{poly()} and \code{splines::ns()} / \code{bs()} are
#' allowed.
#'
#' @param text Covariate expression(s).
#' @return A one-sided formula.
#' @keywords internal
#' @noRd
PS_Formula <- function(text) {
  parsed <- parse(text = paste("~", paste(text, collapse = " + ")))
  if (length(parsed) != 1L) stop("Supply one formula expression only.")
  expr <- parsed[[1]]
  allowed <- c("~", "+", "-", "*", ":", "^", "/", "(", "I", "log", "sqrt", "factor", "poly", "ns", "bs")
  check <- function(x) {
    if (is.call(x)) {
      if (is.call(x[[1]])) {
        h <- x[[1]]
        if (!identical(as.character(h[[1]]), "::") || as.character(h[[2]]) != "splines" || !as.character(h[[3]]) %in% c("ns", "bs"))
          stop("Only splines::ns and splines::bs namespace calls are supported.")
      } else if (!as.character(x[[1]]) %in% allowed) stop("Unsupported formula function: ", as.character(x[[1]]))
      for (v in as.list(x)[-1]) check(v)
    }
  }
  check(expr)
  env <- new.env(parent = baseenv()); env$ns <- splines::ns; env$bs <- splines::bs
  env$poly <- stats::poly
  as.formula(expr, env = env)
}

#' Base variable names in a covariate specification
#'
#' @param x Covariate expression(s).
#' @return A character vector of column names.
#' @keywords internal
#' @noRd
PS_BaseVars <- function(x) if (!length(x) || all(!nzchar(x))) character(0) else unique(all.vars(PS_Formula(x)))

#' Is a setup cell empty?
#'
#' @param x Values.
#' @return \code{TRUE} for \code{NA}, empty or \code{"NA"} values.
#' @keywords internal
#' @noRd
PS_Optional <- function(x) is.na(x) | trimws(as.character(x)) %in% c("", "NA")

#' Parse a covariate caliper specification
#'
#' @param text Text such as \code{"AGE = 5, BMI = 2"}.
#' @param DT Data frame (the variables must be numeric columns).
#' @return A named numeric vector of caliper widths, or \code{NULL}.
#' @keywords internal
#' @noRd
PS_ParseCalipers <- function(text, DT) {
  if (PS_Optional(text)) return(NULL)
  bits <- trimws(strsplit(text, ";", fixed = TRUE)[[1]])
  pairs <- lapply(bits, function(b) trimws(strsplit(b, "=", fixed = TRUE)[[1]]))
  if (any(lengths(pairs) != 2L)) stop("CovariateCalipers must look like AGE = 5; BMI = 2.")
  vars <- vapply(pairs, `[`, character(1), 1)
  widths <- suppressWarnings(as.numeric(vapply(pairs, `[`, character(1), 2)))
  if (anyDuplicated(vars) || any(!vars %in% names(DT)) || any(!is.finite(widths) | widths <= 0))
    stop("CovariateCalipers require unique existing variables and positive numeric widths.")
  if (any(!vapply(DT[, vars, drop=FALSE], is.numeric, logical(1)))) stop("Caliper variables must be numeric.")
  setNames(widths, vars)
}

#' Example setup rows
#'
#' Rows that show the setup options (ratio, caliper, quadratic and spline
#' terms, covariate calipers, interactions, outcome adjustment), built from the
#' setup's own covariates: the first continuous one and the first 0/1 one other
#' than the exposure.
#'
#' @param DT Data frame.
#' @param setup Setup table.
#' @return A data frame of example rows with a \code{Notes} column.
#' @keywords internal
#' @noRd
PS_SetupExamples <- function(DT, setup) {
  base <- setup[1, , drop=FALSE]
  pool <- intersect(unique(unlist(lapply(setup$MatchingVariable, PS_BaseVars))), names(DT))
  pick <- function(test) { hit <- pool[vapply(pool, function(v) test(DT[[v]]), logical(1))]; if (length(hit)) hit[1] else NA_character_ }
  cv <- pick(function(x) is.numeric(x) && length(unique(x[!is.na(x)])) > 2)
  bv <- setdiff(pool[vapply(pool, function(v) all(DT[[v]] %in% c(0, 1, NA)), logical(1))], base$ExposureVariable)[1]
  exposure <- base$ExposureVariable
  out <- list()
  add <- function(note, match = base$MatchingVariable, independent = base$IndependentVariable,
                  exposure = base$ExposureVariable, ratio = NA_real_, ps = NA_real_, cal = NA_character_) {
    r <- base; r$MatchingVariable <- match; r$IndependentVariable <- independent; r$ExposureVariable <- exposure
    r$Ratio <- ratio; r$PSCaliperSD <- ps; r$CovariateCalipers <- cal; r$Notes <- note
    r$ModelID <- paste0("Example_", length(out)+1L)
    r$flagged_covariates <- ""; r$covariate_flags <- ""; r$has_covariate_flags <- NA
    out[[length(out)+1L]] <<- r
  }
  add("Use function defaults for ratio and propensity caliper; no covariate-specific limit.")
  add("Request two controls and a tighter propensity-score caliper.", ratio=2, ps=.1)
  if (!is.na(cv)) {
    width <- signif(stats::sd(DT[[cv]], na.rm = TRUE) / 4, 2)   # a quarter of its SD
    add(sprintf("Quadratic %s relationship in the propensity model.", cv), match=sprintf("%s + I(%s^2)", cv, cv))
    add(sprintf("Natural spline for %s in the propensity model.", cv), match=sprintf("splines::ns(%s, df = 3)", cv))
    add(sprintf("Both the propensity-score caliper and a maximum %s difference of %g apply.", cv, width),
        match=cv, ratio=3, ps=.2, cal=sprintf("%s = %g", cv, width))
  }
  if (!is.na(cv) && !is.na(bv)) add(sprintf("%s by %s interaction in the propensity model.", cv, bv), match=sprintf("%s * %s", cv, bv))
  if (!is.na(bv)) {
    add(sprintf("Include %s as an additional outcome predictor.", bv), independent=sprintf("%s + %s", exposure, bv))
    add(sprintf("Outcome interaction: report the %s contrast standardized over exposed patients.", exposure),
        independent=sprintf("%s * %s", exposure, bv))
  }
  do.call(rbind, out)
}

#' Build the setup table for PSRegression, IPTWRegression and EBRegression
#'
#' One row per dependent x independent x matching specification: each matching
#' covariate on its own (\code{single}) and all of them together
#' (\code{joint}). Rows from \code{MatchingSetupFromLasso()} can be placed
#' first. The function writes nothing: save the table yourself, for example
#' with \code{openxlsx::write.xlsx(list(Setup = setup), file)}, edit it, and
#' read the \code{Setup} sheet back before fitting.
#'
#' With \code{examples = TRUE}, illustrative rows (ratio, caliper, quadratic and
#' spline terms, covariate calipers, interactions, outcome adjustment) are
#' appended to the same table, with a \code{Notes} column explaining each. They
#' are ordinary rows and are fitted like any other, so delete the ones you do
#' not want before reading the table back.
#'
#' @param DT Data frame.
#' @param dependents Outcome columns.
#' @param independents 0/1 exposure columns (terms after the exposure become
#'   outcome-model covariates).
#' @param matching Matching covariate specifications.
#' @param single Logical; one row per matching covariate.
#' @param joint Logical; one row with all matching covariates.
#' @param method \code{"glm"}, \code{"svyglm"} or \code{"lm"}.
#' @param type \code{NULL} (0/1 outcomes are Logistic; all others are Linear),
#'   one of \code{"Logistic"}, \code{"Linear"} or \code{"Poisson"}, or a
#'   named vector by dependent. Poisson must be explicit.
#' @param flags Optional \code{MatchingAssumptionCheck()$table} to annotate the
#'   rows (nothing is removed).
#' @param examples Logical; append illustrative rows, explained in a
#'   \code{Notes} column, at the end of the table (see Description).
#' @param recommended Optional rows from \code{MatchingSetupFromLasso()}.
#' @return The setup data frame; with \code{examples = TRUE} it also has a
#'   \code{Notes} column and the example rows at the end.
#' @examples
#' \dontrun{
#' setup <- MatchingRegressionSetup(model_dt, dependents = c("Death_01", "LOS"),
#'   independents = "sud_any_01", matching = c("AGE", "sex_01"),
#'   method = "glm", type = c(LOS = "Linear"))
#' }
#' @export
MatchingRegressionSetup <- function(DT, dependents, independents, matching, single = TRUE, joint = TRUE, method = "glm", type = NULL, flags = NULL, examples = FALSE,
                                    recommended = NULL) {
  DT <- as.data.frame(DT)
  for (vars in list(dependents, independents, matching))
    if (!is.character(vars) || !length(vars) || anyNA(vars) || any(!nzchar(trimws(vars))))
      stop("dependents, independents and matching must contain nonblank column names.")
  if (length(method) != 1L || is.na(method)) stop("method must be one nonmissing value.")
  matching <- unique(matching)
  if (!length(matching)) stop("Supply at least one matching covariate.")
  missing_cols <- setdiff(c(dependents, unlist(lapply(independents, PS_BaseVars)), unlist(lapply(matching, PS_BaseVars))), names(DT))
  if (length(missing_cols)) stop("Columns not found in DT: ", paste(missing_cols, collapse = ", "))
  if (!method %in% c("glm", "svyglm", "lm")) stop("method must be glm, svyglm or lm.")
  if (!single && !joint) stop("single and joint cannot both be FALSE.")
  combos <- function(mv) expand.grid(DependentVariable = dependents, IndependentVariable = independents,
                                     MatchingVariable = mv, stringsAsFactors = FALSE)
  setup <- rbind(if (single) combos(matching), if (joint) combos(paste(matching, collapse = " + ")))
  setup <- setup[order(setup$DependentVariable, setup$IndependentVariable), ]
  setup$Method <- method
  setup$Type <- PSResolveType(DT, setup$DependentVariable, type)
  if (!all(setup$Type %in% c("Logistic", "Linear", "Poisson"))) stop("type must be Logistic, Linear or Poisson.")
  if (method == "lm" && any(setup$Type != "Linear")) stop("method lm needs type Linear for every dependent.")
  setup <- unique(setup)
  setup <- data.frame(ModelID = seq_len(nrow(setup)), setup, row.names = NULL)
  setup$ExposureVariable <- vapply(setup$IndependentVariable, function(x) PS_BaseVars(x)[1], character(1))
  setup$Ratio <- NA_real_; setup$PSCaliperSD <- NA_real_; setup$CovariateCalipers <- NA_character_
  if (!is.null(recommended) && nrow(recommended)) {
    setup <- as.data.frame(data.table::rbindlist(list(recommended, setup), fill = TRUE))
    if (anyDuplicated(as.character(setup$ModelID))) stop("recommended ModelIDs clash with the generated ones.")
  }
  setup <- PSRegressionFlags(setup, flags)
  drop_empty_flags <- function(x) {
    flag_columns <- intersect(c("flagged_covariates", "covariate_flags"), names(x))
    if (all(vapply(x[flag_columns], function(v) all(is.na(v) | trimws(as.character(v)) == ""), logical(1)))) x[flag_columns] <- NULL
    x
  }
  setup <- drop_empty_flags(setup)
  if (examples) {
    extra <- PS_SetupExamples(DT, setup)
    extra$ModelID <- tail(make.unique(c(as.character(setup$ModelID), as.character(extra$ModelID))), nrow(extra))
    setup$Notes <- ""
    setup <- drop_empty_flags(as.data.frame(data.table::rbindlist(list(setup, extra), fill = TRUE)))
  }
  setup
}
