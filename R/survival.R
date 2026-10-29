#' Kaplan-Meier and Cox proportional-hazards survival analysis
#'
#' Derives a survival time and status from admission/death/discharge date
#' columns, fits an overall Kaplan-Meier curve, optionally a curve stratified
#' by \code{StratifyBy} with a log-rank test, and a batch of univariate Cox
#' proportional-hazards models across \code{Covariates}.
#'
#' @param dt A data.table containing \code{AdmitCol}, \code{DeathCol}, and
#'   (if supplied) \code{DischargeCol}, \code{StratifyBy}, and
#'   \code{Covariates}.
#' @param AdmitCol Column name of the admission date/time.
#' @param DeathCol Column name of the death date/time; \code{NA} or blank is
#'   treated as censored (alive/unknown).
#' @param DischargeCol Optional column name of the discharge date/time, used
#'   as the censoring time for patients without a recorded death.
#' @param StratifyBy Optional column name to stratify the Kaplan-Meier curve
#'   and log-rank test by.
#' @param Covariates Optional character vector of column names to screen via
#'   univariate Cox models.
#' @param DateFormat Default \code{\link[base]{strptime}} format string used
#'   for \code{AdmitCol}, \code{DeathCol}, and \code{DischargeCol} alike.
#' @param AdmitFormat,DeathFormat,DischargeFormat Per-column format string
#'   overrides (each defaults to \code{DateFormat}). Real extracts often mix
#'   a date-only death field with datetime admission/discharge stamps. Set
#'   these arguments explicitly when the column formats differ. Values that
#'   fail to parse under \code{DeathFormat} produce a warning and are treated
#'   as censored observations.
#' @param TimeUnits One of \code{"days"} or \code{"hours"} for the derived
#'   survival time.
#' @param Fpath Optional directory to write the univariate Cox results table
#'   to (as \code{SurvivalUnivariateCox.xlsx}).
#' @return A list including \code{KMOverall}/\code{KMOverallPlot}, and (when
#'   applicable) \code{KMStratified}/\code{KMStratifiedPlot},
#'   \code{LogRankTest}/\code{LogRankP}, and \code{UnivariateCox} (a
#'   data.table of HR, CI, Wald statistic, and p-value per covariate).
#' @examples
#' \dontrun{
#' load(system.file("extdata", "dfMigrant.rda", package = "dataMineR"))
#' SurvivalCompileR(dfMigrant, AdmitCol = "DateofAdmission",
#'                   DeathCol = "DateofDeath", DischargeCol = "DateofDischarge",
#'                   StratifyBy = "Migrant", DeathFormat = "%m/%d/%Y",
#'                   Covariates = c("Age", "ICUadmission"), TimeUnits = "days")
#' }
#' @export
SurvivalCompileR <- function(dt, AdmitCol, DeathCol, DischargeCol = NULL, StratifyBy = NULL,
                              Covariates = NULL, DateFormat = "%m/%d/%Y %H:%M",
                              AdmitFormat = DateFormat, DeathFormat = DateFormat, DischargeFormat = DateFormat,
                              TimeUnits = c("days", "hours"), Fpath = NULL){
  if(!requireNamespace("ggfortify", quietly = TRUE)){ stop("Package 'ggfortify' is required for SurvivalCompileR().") }
  TimeUnits <- match.arg(TimeUnits)
  dt <- copy(as.data.table(dt))
  required_cols <- c(AdmitCol, DeathCol, DischargeCol, StratifyBy, Covariates)
  missing_cols <- setdiff(required_cols, names(dt))
  if(length(missing_cols) > 0L){ stop("Columns not found in dt: ", paste(missing_cols, collapse = ", ")) }
  if(!is.null(Fpath) && !dir.exists(Fpath)){ stop("Fpath does not exist: ", Fpath) }
  admit <- as.POSIXct(dt[[AdmitCol]], format = AdmitFormat)
  death_raw <- dt[[DeathCol]]
  status <- as.integer(!is.na(death_raw) & nzchar(trimws(as.character(death_raw))))
  if(any(status == 0L) && is.null(DischargeCol)){
    stop("DischargeCol is required to calculate censoring times for records without a death date.") }
  death_time <- suppressWarnings(as.POSIXct(as.character(death_raw), format = DeathFormat))
  bad_death_parse <- status == 1 & is.na(death_time)
  if(any(bad_death_parse)){
    warning(sprintf(paste0("%d recorded death(s) in '%s' failed to parse with DeathFormat = '%s' ",
                           "and will be dropped from the survival analysis instead of counted as ",
                           "events. Check the date format (for example, DeathCol may lack a time ",
                           "component that DeathFormat expects)."),
                    sum(bad_death_parse), DeathCol, DeathFormat))
  }
  if(!is.null(DischargeCol) && DischargeCol %in% names(dt)){
    censor_time <- as.POSIXct(dt[[DischargeCol]], format = DischargeFormat)
  } else {
    censor_time <- as.POSIXct(rep(NA_character_, nrow(dt)), format = DischargeFormat)
  }
  event_time <- fifelse(status == 1, death_time, censor_time)
  divisor <- if(TimeUnits == "days"){ 86400 } else {3600}
  dt[, SurvTime := as.numeric(difftime(event_time, admit, units = "secs")) / divisor]
  dt[, SurvStatus := status]
  dt <- dt[!is.na(SurvTime) & SurvTime >= 0,]
  if(nrow(dt) == 0L){ stop("No valid non-negative survival times remain after date parsing.") }
  #### Overall KM ####
  km_overall <- survfit(Surv(SurvTime, SurvStatus) ~ 1, data = dt)
  km_overall_plot <- ggplot2::autoplot(km_overall) +
    ggplot2::labs(title = "Overall survival", x = paste0("Time (", TimeUnits, ")"), y = "Survival probability")
  out <- list(KMOverall = km_overall, KMOverallPlot = km_overall_plot)
  if(!is.null(StratifyBy) && data.table::uniqueN(dt[[StratifyBy]], na.rm = TRUE) >= 2L){
    fmla <- stats::as.formula(paste0("Surv(SurvTime, SurvStatus) ~ `", StratifyBy, "`"))
    km_strat <- survfit(fmla, data = dt)
    km_strat_plot <- ggplot2::autoplot(km_strat) +
      ggplot2::labs(title = paste("Survival by", StratifyBy), x = paste0("Time (", TimeUnits, ")"), y = "Survival probability")
    logrank <- survdiff(fmla, data = dt)
    logrank_p <- stats::pchisq(logrank$chisq, length(logrank$n) - 1, lower.tail = FALSE)
    out <- c(out, list(KMStratified = km_strat, KMStratifiedPlot = km_strat_plot,
                       LogRankTest = logrank, LogRankP = logrank_p))
  } else if (!is.null(StratifyBy)) {
    message("Stratified analysis skipped because fewer than two non-missing strata remain.")
  }
  if(!is.null(Covariates)){
    Covariates <- intersect(Covariates, names(dt))
    if(length(Covariates) == 0L){
      univ_results <- data.table(covariate = character(), term = character(),
                                 HR = numeric(), HR_lower = numeric(), HR_upper = numeric(),
                                 wald_test = numeric(), p_value = numeric())
    } else {univ_results <- rbindlist(lapply(Covariates, function(cov){
      fmla <- stats::as.formula(paste0("Surv(SurvTime, SurvStatus) ~ `", cov, "`"))
      tryCatch({
        model <- coxph(fmla, data = dt)
        s <- summary(model)
        data.table(covariate = cov, term = rownames(s$coef),
                   HR = signif(s$coef[, "exp(coef)"], 3),
                   HR_lower = signif(s$conf.int[, "lower .95"], 3),
                   HR_upper = signif(s$conf.int[, "upper .95"], 3),
                   wald_test = signif(s$coef[, "z"], 3),
                   p_value = signif(s$coef[, "Pr(>|z|)"], 4))
      }, error = function(e) data.table(covariate = cov, term = NA_character_,
                                        HR = NA_real_, HR_lower = NA_real_, HR_upper = NA_real_,
                                        wald_test = NA_real_, p_value = NA_real_))
    }), fill = TRUE)
    }
    univ_results <- univ_results[order(p_value)]
    out$UnivariateCox <- univ_results
    if(!is.null(Fpath)) openxlsx::write.xlsx(list(UnivariateCox = univ_results),
                                             file = file.path(Fpath, "SurvivalUnivariateCox.xlsx"), colWidths = "auto")
  }
  return(out)
}
