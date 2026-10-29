#' Batch single-path causal mediation analysis
#'
#' For each row of \code{MediationTable}, fits the independent-on-dependent,
#' independent-on-mediator, and mediator+independent-on-dependent linear
#' models, then runs nonparametric bootstrap causal mediation analysis via
#' \code{\link[mediation]{mediate}} (500 simulations), reporting ACME, ADE,
#' total effect, and proportion mediated. For simultaneous multi-mediator
#' models, see \code{\link{MultiMediationPLSR}}.
#'
#' @param DT A data.table containing every dependent, independent, and
#'   mediator column referenced in \code{MediationTable}.
#' @param MediationTable A data.frame/data.table with columns
#'   \code{DependentVariable}, \code{IndependentVariable}, and
#'   \code{PutativeMediator}, one row per mediation model to fit.
#' @return A list with \code{IndependentModels} (the three component linear
#'   model summaries per row) and \code{MediationAnalysis} (the compiled
#'   ACME/ADE/total-effect/proportion-mediated table).
#' @examples
#' \dontrun{
#' load(system.file("extdata", "dfMigrant.rda", package = "dataMineR"))
#' MediationTable <- data.table::data.table(
#'   DependentVariable = "ICUadmission", IndependentVariable = "Migrant",
#'   PutativeMediator = "hpi2score")
#' MediationTestR(dfMigrant, MediationTable)
#' }
#' @export
MediationTestR <- function(DT, MediationTable){
  if(!requireNamespace("mediation", quietly = TRUE)){ stop("Package 'mediation' is required for MediationTestR().") }
  DT <- copy(as.data.table(DT))
  MediationTable <- as.data.table(MediationTable)
  required_cols <- c("DependentVariable", "IndependentVariable", "PutativeMediator")
  if (!all(required_cols %in% names(MediationTable))){ stop("MediationTable must contain: ", paste(required_cols, collapse = ", ")) }
  extract_mediation_summary <- function (x) {
    clp <- 100 * x$conf.level
    isLinear.y <- ((class(x$model.y)[1] %in% c("lm", "rq")) ||
                     (inherits(x$model.y, "glm") && x$model.y$family$family ==
                        "gaussian" && x$model.y$family$link == "identity") ||
                     (inherits(x$model.y, "survreg") && x$model.y$dist ==
                        "gaussian"))
    printone <- !x$INT && isLinear.y
    if (printone) {
      smat <- c(x$d1, x$d1.ci, x$d1.p)
      smat <- rbind(smat, c(x$z0, x$z0.ci, x$z0.p))
      smat <- rbind(smat, c(x$tau.coef, x$tau.ci, x$tau.p))
      smat <- rbind(smat, c(x$n0, x$n0.ci, x$n0.p))
      rownames(smat) <- c("ACME", "ADE", "Total Effect", "Prop. Mediated")
    } else {
      smat <- c(x$d0, x$d0.ci, x$d0.p)
      smat <- rbind(smat, c(x$d1, x$d1.ci, x$d1.p))
      smat <- rbind(smat, c(x$z0, x$z0.ci, x$z0.p))
      smat <- rbind(smat, c(x$z1, x$z1.ci, x$z1.p))
      smat <- rbind(smat, c(x$tau.coef, x$tau.ci, x$tau.p))
      smat <- rbind(smat, c(x$n0, x$n0.ci, x$n0.p))
      smat <- rbind(smat, c(x$n1, x$n1.ci, x$n1.p))
      smat <- rbind(smat, c(x$d.avg, x$d.avg.ci, x$d.avg.p))
      smat <- rbind(smat, c(x$z.avg, x$z.avg.ci, x$z.avg.p))
      smat <- rbind(smat, c(x$n.avg, x$n.avg.ci, x$n.avg.p))
      rownames(smat) <- c("ACME (control)", "ACME (treated)",
                          "ADE (control)", "ADE (treated)", "Total Effect",
                          "Prop. Mediated (control)", "Prop. Mediated (treated)",
                          "ACME (average)", "ADE (average)", "Prop. Mediated (average)")  }
    colnames(smat) <- c("Estimate", paste(clp, "% CI Lower", sep = ""), paste(clp, "% CI Upper", sep = ""), "p-value")
    smat }
  DependentVariable <- as.character(MediationTable$DependentVariable)
  IndependentVariable <- as.character(MediationTable$IndependentVariable)
  PutativeMediator <- as.character(MediationTable$PutativeMediator)
  model_variables <- unique(c(DependentVariable, IndependentVariable, PutativeMediator))
  missing_variables <- setdiff(model_variables, names(DT))
  if (length(missing_variables) > 0L){ stop("Model variables not found in DT: ", paste(missing_variables, collapse = ", ")) }
  #### Perform analysis ####
  IndependentDT <- data.table()
  MediationDT <- data.table()
  for(i in seq_len(nrow(MediationTable))){
    variables <- c(DependentVariable[i], IndependentVariable[i], PutativeMediator[i])
    if(anyDuplicated(variables)){stop("Dependent, independent, and mediator variables must be distinct.")}
    DTtemp <- copy(DT[, variables, with = FALSE])
    setnames(DTtemp, c(".outcome", ".treatment", ".mediator"))
    DTtemp <- DTtemp[complete.cases(DTtemp),]
    if(nrow(DTtemp) < 5L){stop("Each mediation model requires at least five complete observations.") }
    #### Independent prediction on dependent variable ####
    fmla0 <- ".outcome ~ .treatment"
    model.0 <- lm(as.formula(fmla0), data = DTtemp)
    dt0 <- data.table(coef(summary(model.0)))[2,]
    #### Independent prediction on putative mediator ####
    fmla1 <- ".mediator ~ .treatment"
    model.1 <- lm(as.formula(fmla1), data = DTtemp)
    dt1 <- data.table(coef(summary(model.1)))[2,]
    #### putative mediator prediction on dependent variable ####
    fmla2 <- ".outcome ~ .mediator + .treatment"
    model.2 <- lm(as.formula(fmla2), data = DTtemp)
    dt2 <- data.table(coef(summary(model.2)))[2,]
    #### Create table for independent models ####
    IndependentModels <- rbind(dt0, dt1, dt2)
    IndependentModels$Fmla <- c(
      paste(DependentVariable[i], "~", IndependentVariable[i]),
      paste(PutativeMediator[i], "~", IndependentVariable[i]),
      paste(DependentVariable[i], "~", PutativeMediator[i], "+", IndependentVariable[i])
    )
    IndependentModels$DependMedInd <- paste(DependentVariable[i], PutativeMediator[i], IndependentVariable[i], sep = "-")
    IndependentDT <- rbind(IndependentDT, IndependentModels)
    #### Perform mediation analysis ####
    results1 <- mediation::mediate(model.1, model.2, treat=".treatment", mediator=".mediator", boot=TRUE, sims=500)
    med <- as.data.frame(extract_mediation_summary(summary(results1)))
    med$Effects <- row.names(med)
    med$DependMedInd <- paste(DependentVariable[i], PutativeMediator[i], IndependentVariable[i], sep = "-")
    MediationDT <- rbind(MediationDT, med)
  }
  return(list(IndependentModels = IndependentDT, MediationAnalysis = MediationDT))
}

#' Multi-mediator PLS-SEM mediation analysis
#'
#' Builds a partial-least-squares structural equation model with single-item
#' composite constructs for the treatment, one or more candidate mediators,
#' optional covariates, and the outcome; estimates path coefficients via
#' \pkg{seminr}, bootstraps the model, and tests the specific indirect
#' (mediated) effect through each mediator. A multi-mediator complement to
#' the single-path \code{\link{MediationTestR}}.
#'
#' @param dt A data.frame/data.table containing \code{Treatment},
#'   \code{Mediators}, \code{Outcome}, and \code{Covariates}, with complete
#'   cases on those columns.
#' @param Treatment Column name of the treatment/exposure variable.
#' @param Mediators Character vector of candidate mediator column names.
#' @param Outcome Column name of the outcome variable.
#' @param Covariates Optional character vector of covariate column names
#'   also modeled as direct predictors of \code{Outcome}.
#' @param Nboot Number of bootstrap resamples for \code{seminr::bootstrap_model}.
#' @param Seed Random seed for the bootstrap.
#' @return A list with \code{Model} (the fitted \code{seminr} PLS model),
#'   \code{Bootstrap} (the bootstrapped model), and \code{IndirectEffects} (a
#'   compiled table of specific indirect-effect significance per mediator).
#' @examples
#' \dontrun{
#' load(system.file("extdata", "dfMigrant.rda", package = "dataMineR"))
#' dt <- dfMigrant[!is.na(economic) & !is.na(education) & !is.na(housing), ]
#' MultiMediationPLSR(dt, Treatment = "Migrant",
#'                     Mediators = c("economic", "education", "housing"),
#'                     Outcome = "ICUadmission", Nboot = 100)
#' }
#' @export
MultiMediationPLSR <- function(dt, Treatment, Mediators, Outcome, Covariates = NULL, Nboot = 1000, Seed = 1){
  if(!requireNamespace("seminr", quietly = TRUE)){stop("Package 'seminr' is required for MultiMediationPLSR().")}
  dt <- as.data.frame(dt)
  construct_names <- unique(c(Treatment, Outcome, Mediators, Covariates))
  missing_cols <- setdiff(construct_names, names(dt))
  if(length(missing_cols) > 0L){ stop("Columns not found in dt: ", paste(missing_cols, collapse = ", ")) }
  if(length(Mediators) == 0L){ stop("Mediators must contain at least one column.") }
  if(!is.numeric(Nboot) || length(Nboot) != 1L || Nboot < 1L || Nboot != as.integer(Nboot)){ stop("Nboot must be a positive integer.") }
  if(any(!stats::complete.cases(dt[, construct_names, drop = FALSE]))){ stop("dt must contain complete cases for all model constructs.")}
  constructs_list <- lapply(construct_names, function(v) seminr::composite(v, seminr::single_item(v)))
  mm <- do.call(seminr::constructs, constructs_list)
  paths_list <- c(list(seminr::paths(from = Treatment, to = c(Outcome, Mediators))),
    lapply(Mediators, function(m) seminr::paths(from = m, to = Outcome)),
    if(!is.null(Covariates)){list(seminr::paths(from = Covariates, to = Outcome)) } else {NULL}
  )
  sm <- do.call(seminr::relationships, paths_list)
  model <- seminr::estimate_pls(data = dt, measurement_model = mm, structural_model = sm)
  boot <- seminr::bootstrap_model(model, nboot = Nboot, cores = NULL, seed = Seed)
  IndirectEffects <- rbindlist(lapply(Mediators, function(m){
    eff <- as.data.frame(seminr::specific_effect_significance(boot, from = Treatment, to = Outcome, through = m))
    eff$Mediator <- m
    as.data.table(eff, keep.rownames = "Effect")
  }), fill = TRUE)
  return(list(Model = model, Bootstrap = boot, IndirectEffects = IndirectEffects))
}
