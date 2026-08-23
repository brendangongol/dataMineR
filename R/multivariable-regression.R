#' Multivariable regression sweep on significant univariate predictors
#'
#' For each dependent variable, takes its significant predictors (from
#' \code{IndependentPredictorDT}, typically the output of
#' \code{\link{RegressionCompile}}) and fits every pairwise-combination
#' multivariable model via \code{\link{glmcompileR}}.
#'
#' @param DataTable A data.table containing every dependent and independent
#'   column referenced in \code{IndependentPredictorDT}.
#' @param deps Character vector of dependent variable names to process.
#' @param IndependentPredictorDT A data.table with (at least) columns
#'   \code{dependent}, \code{independent}, and \code{Significance}, as
#'   produced by \code{\link{RegressionCompile}}.
#' @param UseSignificantPredictorsOnly Logical; if \code{TRUE}, restrict to
#'   rows where \code{Significance == "Significant"}.
#' @param parametric Logical; passed through to \code{\link{glmcompileR}}.
#' @param CalculateAUC Logical; passed through to \code{\link{glmcompileR}}
#'   (logistic models only).
#' @param Type One of \code{"Linear"} or \code{"Logistic"}.
#' @param ReturnSignificantOnly Character \code{"T"}/\code{"F"}, passed
#'   through to \code{\link{glmcompileR}}.
#' @return A data.table of multivariable model results across every
#'   dependent with at least two significant predictors.
#' @export
MultiModelRegression <- function(DataTable, deps, IndependentPredictorDT, UseSignificantPredictorsOnly = TRUE, parametric = FALSE, CalculateAUC = FALSE, Type = "Linear", ReturnSignificantOnly = "F"){
  DataTable <- as.data.table(DataTable)
  IndependentPredictorDT <- copy(as.data.table(IndependentPredictorDT))
  required_cols <- c("dependent", "independent", "Significance")
  if(!all(required_cols %in% names(IndependentPredictorDT))){ stop("IndependentPredictorDT must contain: ", paste(required_cols, collapse = ", ")) }
  if(!Type %in% c("Linear", "Logistic")){ stop("Type must be 'Linear' or 'Logistic'.") }
  if(UseSignificantPredictorsOnly){ IndependentPredictorDT = IndependentPredictorDT[Significance == "Significant",] }
  MultiManuscriptModels <- data.table()
  for(i in seq_along(deps)){
    independents <- IndependentPredictorDT[dependent == deps[i],]
    independents <- independents$independent
    if(length(independents) > 1){
      DTManuscript <- DataTable[,unique(c(independents, deps[i])), with = FALSE]
      k <- ComparisonSetup(n=(ncol(DTManuscript)-1), MaxComps = 2)
      multiModels <- glmcompileR(DT = DTManuscript, key = k, significant = ReturnSignificantOnly, AUC = CalculateAUC, Odds = TRUE,
                                 DependentColumn = deps[i], type = Type, Parametric = parametric)
      MultiManuscriptModels <- rbind(MultiManuscriptModels, multiModels, fill = TRUE)
    } else {
      print(paste("There are not enough independent variables for", deps[i]))
    }
  }
  return(MultiManuscriptModels)
}

#' Confounding ("inhibition") screen across paired multivariable models
#'
#' For each dependent in \code{CompiledSignificant}, fits every pairwise
#' multivariable combination of its significant predictors (via
#' \code{\link{glmcompileR}}) and identifies predictors whose significance is
#' lost when co-included with another variable, identifying candidate confounders.
#'
#' @param DT A data.table containing every dependent and independent column
#'   referenced in \code{CompiledSignificant}.
#' @param CompiledSignificant A data.table with (at least) columns
#'   \code{dependent} and \code{independent}, typically the significant
#'   subset of \code{\link{RegressionCompile}}'s output.
#' @param AddBack Optional value; dependent-variable names matching this
#'   value are not excluded from another dependent's predictor set (rarely
#'   needed; defaults to \code{NA}, excluding all dependents from each
#'   other's predictor sets).
#' @return A list with \code{PairedRegressionDT} (all pairwise models),
#'   \code{InhibitionTable} (predictor, co-variate, and inhibitor count for
#'   every lost-significance case), and \code{SignificantPredsInhibitionDT}
#'   (\code{CompiledSignificant} annotated with inhibition info).
#' @export
PairedMultiVariableSummary <- function(DT, CompiledSignificant, AddBack = NA){
  DT <- as.data.table(DT)
  CompiledSignificant <- copy(as.data.table(CompiledSignificant))
  if(!all(c("dependent", "independent") %in% names(CompiledSignificant))){ stop("CompiledSignificant must contain 'dependent' and 'independent' columns.") }
  if(nrow(CompiledSignificant) == 0L){
    return(list(PairedRegressionDT = data.table(), InhibitionTable = data.table(),
                SignificantPredsInhibitionDT = CompiledSignificant))
  }
  depends <- unique(CompiledSignificant$dependent)
  PairedRegressionDT <- data.table()
  DTInhibition <- data.table()
  for(i in seq_along(depends)){
    Tselect <- CompiledSignificant[dependent == depends[i],] #### select significant analysis groups ####
    glmDT <- DT[,(colnames(DT) %in% c(unique(Tselect$dependent), unique(Tselect$independent))), with = FALSE] #### Subset master data table ####
    glmDT <- glmDT[,c(colnames(glmDT)[!colnames(glmDT) %in% depends[i]], depends[i]), with = FALSE] #### re-arrange columns ####
    if((ncol(glmDT) - 1) < 2){
      print(paste("There are not enough significant predictors for", depends[i]))
      next
    }
    k <- ComparisonSetup(n=(ncol(glmDT)-1), MaxComps = 2)
    multiModels <- glmcompileR(DT = glmDT, key = k, significant = "F", AUC = FALSE,
                               DependentColumn = depends[i], type = "Linear", Parametric = TRUE)
    PairedRegressionDT <- rbind(PairedRegressionDT, multiModels)
    Nsig <- multiModels[Significance == "Not-Significant",]
    NinsInd <- unique(Nsig$names)
    for(b in seq_along(NinsInd)){
      temp <- Nsig[names == NinsInd[b],]
      InhibitedBy <- gsub(" + ", "", gsub(NinsInd[b], "",  gsub("^.+ ~ ", "", temp$formula)), fixed = TRUE)
      InhibitedBy <- InhibitedBy[!(InhibitedBy %in% depends[!(depends %in% AddBack)])] #### remove dependent variables
      NumberOfInhibitors <- length(InhibitedBy)
      InhibitedBy <- paste(InhibitedBy, collapse = ", ")
      DTInhibition <- rbind(DTInhibition, data.table(dependent = depends[i], independent =NinsInd[b], InhibitedBy = InhibitedBy, NumberOfInhibitors=NumberOfInhibitors ) ) }
  }
  if(nrow(DTInhibition) == 0L){
    DTInhibition <- data.table(dependent = character(), independent = character(),
                               InhibitedBy = character(), NumberOfInhibitors = integer())
  }
  CompiledSignificantInhib <- merge(CompiledSignificant, DTInhibition, by = c("dependent", "independent"), all.x = TRUE)
  (CompiledSignificantInhib <- CompiledSignificantInhib[!(CompiledSignificantInhib$independent %in% depends[!(depends %in% AddBack)]),]) #### remove dependent variables
  return(list(PairedRegressionDT=PairedRegressionDT, InhibitionTable = DTInhibition, SignificantPredsInhibitionDT = CompiledSignificantInhib))
}

#' Bootstrap "stability" of a regression screen via repeated resampling
#'
#' Repeatedly draws a bootstrap sample of \code{dt} with replacement (size
#' \code{Fraction * nrow(dt)}), re-runs \code{\link{RegressionCompile}} on
#' each draw, and tallies how often each dependent-independent pair remains
#' significant. This provides a resampling-based robustness check for variable
#' selection, complementary to a single full-data model fit.
#'
#' @param dt A data.table containing every column in \code{dependents} and
#'   \code{independents}.
#' @param dependents Character vector of outcome column names.
#' @param independents Character vector of candidate predictor column
#'   names.
#' @param Nboots Number of bootstrap resamples.
#' @param Fraction Fraction of \code{dt}'s rows to draw with replacement on
#'   each resample.
#' @param Parametric Logical; passed through to
#'   \code{\link{RegressionCompile}}.
#' @param type One of \code{"Linear"} or \code{"Logistic"}, passed through
#'   to \code{\link{RegressionCompile}}.
#' @param Seed Optional integer random seed used before bootstrap sampling.
#' @return A data.table with columns \code{dependent}, \code{independent},
#'   \code{TimesSignificant}, and \code{StabilityFraction}, ordered from
#'   most to least stable.
#' @export
RegressionStabilityR <- function(dt, dependents, independents, Nboots = 100, Fraction = 0.5,
                                 Parametric = TRUE, type = "Linear", Seed = NULL){
  dt <- as.data.table(dt)
  if(!is.numeric(Nboots) || length(Nboots) != 1L || Nboots < 1 || Nboots != as.integer(Nboots)){ stop("Nboots must be a positive integer.") }
  if(!is.numeric(Fraction) || length(Fraction) != 1L || Fraction <= 0 || Fraction > 1){ stop("Fraction must be greater than 0 and at most 1.") }
  if(!is.null(Seed)){
    if(!is.numeric(Seed) || length(Seed) != 1L || !is.finite(Seed)){
      stop("Seed must be NULL or a finite integer.")
    }
    seed_integer <- suppressWarnings(as.integer(Seed))
    if(is.na(seed_integer) || Seed != seed_integer){
      stop("Seed must be NULL or a finite integer.")
    }
    set.seed(seed_integer)
  }
  if(nrow(dt) < 3L){ stop("dt must contain at least three rows.") }
  Nsamp <- round(nrow(dt) * Fraction)
  if(Nsamp < 3L){ stop("Fraction produces fewer than three rows per bootstrap sample.") }
  StabilityList <- vector("list", Nboots)
  for(b in seq_len(Nboots)){
    samp <- dt[sample.int(.N, Nsamp, replace = TRUE)]
    sig <- tryCatch(
      RegressionCompile(DT = samp, dependents = dependents, independents = independents,
                        Fpath = tempdir(), ReturnSignificant = TRUE, Parametric = Parametric,
                        type = type, WriteFiles = FALSE),
      error = function(e) data.table())
    if(nrow(sig) > 0) StabilityList[[b]] <- unique(sig[, .(dependent, independent)])
  }
  AllPairs <- rbindlist(StabilityList, fill = TRUE)
  if(nrow(AllPairs) == 0){
    return(data.table(dependent = character(), independent = character(),
                      TimesSignificant = integer(), StabilityFraction = numeric()))
  }
  StabilityDT <- AllPairs[, .N, by = .(dependent, independent)]
  setnames(StabilityDT, "N", "TimesSignificant")
  StabilityDT[, StabilityFraction := round(TimesSignificant / Nboots, 3)]
  setorder(StabilityDT, -StabilityFraction)
  return(StabilityDT)
}
