#' ANCOVA main-effects and interaction screening
#'
#' Fits an additive multi-covariate linear model
#' (\code{Dependent ~ Covariates}) and reports its ANOVA table, then --
#' for each pair in \code{InteractionPairs}, fits a model with an explicit
#' two-way interaction term and reports that ANOVA table too, screening for
#' effect modification beyond the single-predictor GLM screens in
#' \code{\link{RegressionCompile}}.
#'
#' @param dt A data.frame/data.table containing \code{Dependent} and every
#'   column in \code{Covariates}.
#' @param Dependent Column name of the outcome.
#' @param Covariates Character vector of covariate column names for the
#'   additive model.
#' @param InteractionPairs Optional list of length-2 character vectors, each
#'   giving a pair of covariate names to test for a two-way interaction.
#' @return A list with \code{MainEffects} (the additive model's ANOVA
#'   table), \code{MainModel} (the fitted \code{lm} object), and, when
#'   \code{InteractionPairs} is supplied, \code{InteractionEffects} (the
#'   compiled interaction-model ANOVA tables).
#' @examples
#' \dontrun{
#' data(dfMigrant)
#' AncovaInteractionR(dfMigrant, Dependent = "IndexSOFAscoreattimeofsepsis",
#'                     Covariates = c("Age", "Migrant", "Homeless"),
#'                     InteractionPairs = list(c("Age", "Migrant")))
#' }
#' @export
AncovaInteractionR <- function(dt, Dependent, Covariates, InteractionPairs = NULL){
  dt <- as.data.frame(dt)
  missing_cols <- setdiff(c(Dependent, Covariates), names(dt))
  if (length(missing_cols) > 0L) stop("Columns not found in dt: ", paste(missing_cols, collapse = ", "))
  if (length(Covariates) == 0L) stop("Covariates must contain at least one column.")
  quoted_covariates <- paste0("`", Covariates, "`")
  fmla_main <- stats::as.formula(paste0("`", Dependent, "` ~ ", paste(quoted_covariates, collapse = " + ")))
  model_main <- stats::lm(fmla_main, data = dt)
  MainEffects <- as.data.table(stats::anova(model_main), keep.rownames = "term")
  MainEffects[, model := "additive"]
  Results <- list(MainEffects = MainEffects, MainModel = model_main)
  if(!is.null(InteractionPairs)){
    IntTables <- rbindlist(lapply(InteractionPairs, function(pr){
      if (length(pr) != 2L || !all(pr %in% Covariates))
        stop("Each InteractionPairs entry must contain two names from Covariates.")
      interaction <- paste0("`", pr[1L], "`:`", pr[2L], "`")
      fmla_int <- stats::as.formula(paste0("`", Dependent, "` ~ ", paste(c(quoted_covariates, interaction), collapse = " + ") ))
      m <- stats::lm(fmla_int, data = dt)
      tab <- as.data.table(stats::anova(m), keep.rownames = "term")
      tab[, model := paste(pr, collapse = " x ")]
      tab
    }), fill = TRUE)
    Results$InteractionEffects <- IntTables
  }
  return(Results)
}
