#' Build a comparison grid of column-index combinations
#'
#' Generates the underscore-joined column-index keys consumed by
#' \code{\link{glmcompileR}} for a multivariable model sweep: either every
#' size-\code{MaxComps} combination of \code{n} columns, or (when
#' \code{MaxComps = "ALL"}) every non-empty subset.
#'
#' @param n Number of independent-variable columns available.
#' @param MaxComps Either an integer combination size, or the string
#'   \code{"ALL"} for every non-empty subset (grows as 2^n; use with
#'   caution for large \code{n}).
#' @return A character vector of underscore-joined column-index
#'   combinations, e.g. \code{c("1_2", "1_3", "2_3")}.
#' @export
ComparisonSetup <- function(n=5, MaxComps = 3){
  if(length(n) != 1L || !is.numeric(n) || is.na(n) || n < 1 || n != as.integer(n)){ stop("n must be a positive integer.") }
  n <- as.integer(n)
  if (identical(MaxComps, "ALL")){
    k <- unlist(lapply(seq_len(n), function(size) {
      apply(utils::combn(seq_len(n), size), 2L, paste, collapse = "_")
    }), use.names = FALSE)
  } else {
    if(length(MaxComps) != 1L || !is.numeric(MaxComps) || is.na(MaxComps) || MaxComps < 1 || MaxComps > n || MaxComps != as.integer(MaxComps)) {
      stop("MaxComps must be 'ALL' or an integer between 1 and n.")
    }
    k <- apply(utils::combn(seq_len(n), as.integer(MaxComps)), 2L,
               paste, collapse = "_")
  }
  return(k)}

#' Fit a GLM for every column-combination key (linear or logistic)
#'
#' For each key in \code{key} (an underscore-joined set of column indices
#' into \code{DT}, as produced by \code{\link{ComparisonSetup}}), fits a
#' linear or logistic regression of \code{DependentColumn} on that subset of
#' predictors, optionally computing AUC/odds ratios, and flags each
#' predictor "Significant"/"Not-Significant" at p < 0.05. Used internally by
#' \code{\link{MultiModelRegression}} and \code{\link{PairedMultiVariableSummary}}
#' to sweep multivariable models.
#'
#' @param DT A data.table containing the numeric predictor columns and
#'   \code{DependentColumn}.
#' @param key Character vector of underscore-joined column-index
#'   combinations, as produced by \code{\link{ComparisonSetup}}.
#' @param significant Character \code{"T"}/\code{"F"}; if \code{"T"}, keep
#'   only rows with p < 0.05 (excluding the intercept).
#' @param AUC Logical; compute model AUC (logistic models only).
#' @param Odds Logical; include odds ratios (logistic models only).
#' @param DependentColumn Name of the outcome column.
#' @param type One of \code{"Logistic"} or \code{"Linear"}.
#' @param Parametric Logical; if \code{FALSE}, drop incomplete cases before
#'   fitting.
#' @return A data.table of model coefficients (one row per predictor per
#'   fitted model), with columns including \code{names}, \code{formula}, and
#'   \code{Significance}.
#' @export
glmcompileR <- function(DT, key, significant, AUC = TRUE, Odds = TRUE, DependentColumn = "outcome", type = "Logistic", Parametric = FALSE){
  if(!type %in% c("Linear", "Logistic")){ stop("type must be 'Linear' or 'Logistic'.") }
  if(type == "Linear" && !Parametric && !requireNamespace("npreg", quietly = TRUE)){
    stop("Package 'npreg' is required for non-parametric linear models.") }
  compiled_dt <- data.table()
  pb <- txtProgressBar(min = 0, max = length(key), style = 3)
  for(m in seq_along(key)){
    #### subset data table to contain only required columns
    spl <- strsplit(key[m], split = "_")[[1]]
    spl <- spl[!spl == "0"]
    newdt <- data.table()
    for(i in seq_along(spl)){
      temp <- DT[,c(as.numeric(spl[i])), with = FALSE]
      newdt <- cbind(newdt, temp) }
    #### find rows with all NA values ####
    allNA <- apply(newdt, 1, function(x){sum(is.na(x))}) < ncol(newdt)
    #### add outcome variable to last column
    newdt[[DependentColumn]] <- DT[[DependentColumn]]
    #### remove rows with all NA values ####
    newdt <- newdt[allNA,]
    #### generate formula
    nam <- colnames(newdt[, -(ncol(newdt)), with = FALSE])
    if(length(nam) > 1){
      plus <- "+"
      p1 <- paste(paste(nam[-length(nam)], plus, sep = " "),collapse = " ")
      p2 <- paste(p1, nam[length(nam)])
      fmla <- as.formula(paste0(colnames(newdt)[ncol(newdt)], " ~ ", p2))
    } else { fmla <- as.formula(paste0(colnames(newdt)[ncol(newdt)], " ~ ", nam)) }
    #### fit glm model
    glm.fit <- tryCatch({
      if(type == "Logistic") {
        glm(fmla, data = newdt, family = binomial)
      } else if(Parametric) {
        lm(fmla, data = newdt)
      } else {
        npreg::gsm(formula = fmla, data = newdt)
      }
    }, error = function(e) {
      message(sprintf("Model %d skipped: %s", m, conditionMessage(e)))
      NULL
    })
    if(is.null(glm.fit)){
      setTxtProgressBar(pb, m)
      next
    }
    #### create data table
    if(type == "Logistic"){
      prob <- suppressMessages(predict(glm.fit, type = c("response")))
      if(AUC == TRUE){
        ROC <- suppressMessages(roc(glm.fit$y, prob)) }
      dtf <- data.table(coef(summary(glm.fit)))
      if(AUC == TRUE){
        dtf$AUC <- as.numeric(ROC$auc) }
      if(Odds){
        dtf$OddsRatio <- exp(dtf$Estimate)
        dtf$logOddsRatio <- dtf$Estimate
        dtf$OddsRatioP <- dtf$`Pr(>|z|)`
      }
      dtf$names <- rownames(summary(glm.fit)$coefficients)
      chr <- as.character(glm.fit$formula)
      dtf$formula <- paste(chr[2], chr[1], chr[3:length(chr)])
    }
    if(type == "Linear"){
      if(Parametric){
        dtf <- data.table(coef(summary(glm.fit)))
        dtf$names <- rownames(summary(glm.fit)$coefficients)
      } else {
        dtf2 <- summary(glm.fit)$s.table
        dtf <- dtf2
        dtf$names <- gsub("\\(", "", gsub("\\)", "", gsub("s\\(", "", rownames(dtf))))
      }
      chr <- as.character(fmla)
      dtf$formula <- paste(chr[2], chr[1], chr[3])
      residuals <- glm.fit$residuals[!is.na(glm.fit$residuals)]
      if(length(residuals) >= 3L && length(residuals) <= 5000L){
        test <- shapiro.test(residuals)
        dtf$Residualdistribution <- if (test$p.value < 0.05)
          "Not normally distributed" else "normally distributed"
      } else {
        dtf$Residualdistribution <- "Not tested"
      }
    }
    if(significant == "T"){
      #### subset data table to retain only significant indepenent variables
      p_name <- names(dtf)[grepl("Pr(>", names(dtf), fixed = TRUE)][1L]
      if (!is.na(p_name)){
        dtf <- dtf[dtf[[p_name]] < 0.05 & !dtf$names %in% c("(Intercept)", "Intercept"),] }
    }
    if(nrow(dtf) > 0){
      dtf$NcountsNoNA <-  nrow(newdt[complete.cases(newdt),])
      dtf <- dtf[(!dtf$names == "(Intercept)"),]
      dtf <- dtf[(!dtf$names == "Intercept"),]
      #### Add confidence intervals ####
      Conf <- suppressMessages(stats::confint.default(glm.fit))
      Conf <- Conf[!grepl("knot|Intercept", rownames(Conf)),]
      if(sum(class(Conf) == "matrix") > 0){
        dtf[[colnames(Conf)[[1]]]] <- c(Conf[,1], rep(NA, (nrow(dtf)-length(Conf[,1]))   ))
        dtf[[colnames(Conf)[[2]]]] <- c(Conf[,2], rep(NA, (nrow(dtf)-length(Conf[,2]))   ))
      } else {
        dtf[[names(Conf)[1]]] <- Conf[1]
        dtf[[names(Conf)[2]]] <- Conf[2]
      }
      if(sum(dtf$names == "Residuals") > 0){
        if(sum(class(Conf) == "matrix") > 0){
          dtf[dtf$names == "Residuals",][[colnames(Conf)[1]]] <- NA
          dtf[dtf$names == "Residuals",][[colnames(Conf)[2]]] <- NA
        } else {
          dtf[dtf$names == "Residuals",][[names(Conf)[1]]] <- NA
          dtf[dtf$names == "Residuals",][[names(Conf)[2]]] <- NA
        }
      }
      compiled_dt <- rbind(compiled_dt, dtf) }
    setTxtProgressBar(pb, m)
  }
  close(pb)
  if(nrow(compiled_dt) == 0L){ return(compiled_dt) }
  if(Parametric){
    compiled_dt <- compiled_dt[!duplicated(compiled_dt[,c(1:6)]),]
    compiled_dt <- compiled_dt[(!compiled_dt$names == "(Intercept)"),]
    compiled_dt <- compiled_dt[(!compiled_dt$names == "Intercept"),]
    compiled_dt$Significance <- "NA"
    if(sum(grepl("Pr(>|t|)", colnames(compiled_dt), fixed = TRUE)) > 0){
      compiled_dt[`Pr(>|t|)` < 0.05, Significance := "Significant"]
      compiled_dt[`Pr(>|t|)` >= 0.05, Significance := "Not-Significant"]
    } else {
      compiled_dt[`Pr(>|z|)` < 0.05, Significance := "Significant"]
      compiled_dt[`Pr(>|z|)` >= 0.05, Significance := "Not-Significant"]
    }
  } else {
    compiled_dt <- compiled_dt[(!compiled_dt$names == "(Intercept)"),]
    compiled_dt <- compiled_dt[(!compiled_dt$names == "Intercept"),]
    compiled_dt$Significance <- "NA"
    if(sum(colnames(compiled_dt) %in% "Pr(>F)") > 0){
      compiled_dt[`Pr(>F)` < 0.05, Significance := "Significant"]
      compiled_dt[`Pr(>F)` >= 0.05, Significance := "Not-Significant"]
    }
    if(sum(colnames(compiled_dt) %in% "Pr(>|z|)") > 0){
      compiled_dt[`Pr(>|z|)` < 0.05, Significance := "Significant"]
      compiled_dt[`Pr(>|z|)` >= 0.05, Significance := "Not-Significant"]
    }
  }
  compiled_dt$dependent <- gsub(" ~ .+", "", compiled_dt$formula)
  compiled_dt$independent <- gsub(".+ ~ ", "", compiled_dt$formula)
  return(compiled_dt) }

#' Univariate GLM (linear or logistic) screen across many predictors
#'
#' Fits one univariate model per column in \code{key} (\code{DependentColumn
#' ~ predictor}), optionally computing AUC/odds ratios, and flags each
#' result "Significant"/"Not-Significant" at p < 0.05. This is the
#' single-predictor screening step used inside \code{\link{RegressionCompile}}.
#'
#' @param DT A data.table containing the numeric predictor columns and
#'   \code{DependentColumn}.
#' @param key Character or numeric vector of column indices into \code{DT}
#'   to screen, one model per entry.
#' @param significant Character \code{"T"}/\code{"F"}; if \code{"T"}, keep
#'   only significant results.
#' @param type One of \code{"Logistic"} or \code{"Linear"}.
#' @param AUC Logical; compute model AUC (logistic models only).
#' @param Odds Logical; include odds ratios (logistic models only).
#' @param DependentColumn Name of the outcome column.
#' @param Parametric Logical; if \code{FALSE}, drop incomplete cases before
#'   fitting.
#' @param verbose Logical; show a text progress bar.
#' @param ReferenceName Optional reference-level label recorded alongside
#'   each result.
#' @return A data.table of univariate model results, one row per predictor
#'   in \code{key}.
#' @export
univariantglmR <- function(DT, key, significant, type = "Logistic", AUC = TRUE, Odds = TRUE, DependentColumn = "outcome", Parametric = TRUE, verbose = FALSE, ReferenceName = NA){
  if(!type %in% c("Linear", "Logistic")){ stop("type must be 'Linear' or 'Logistic'.") }
  if(type == "Linear" && !Parametric && !requireNamespace("npreg", quietly = TRUE)){
    stop("Package 'npreg' is required for non-parametric linear models.") }
  if(!DependentColumn %in% names(DT)){ stop("DependentColumn not found in DT: ", DependentColumn) }
  compiled_dt <- data.table()
  if(verbose){pb <- txtProgressBar(min = 0, max = length(key), style = 3)}
  for(m in seq_along(key)){
    #### subset data table to contain only required columns
    newdt <- DT[,c(as.numeric(key[m])), with = FALSE]
    if(length(unique(newdt[[1]])[!is.na(unique(newdt[[1]]))]) > 1){
      #### add outcome variable to last column
      newdt[[DependentColumn]] <- DT[[DependentColumn]]
      if(sum(complete.cases(newdt)) > (round(nrow(newdt)/4))){
        if(!Parametric){
          newdt <- newdt[complete.cases(newdt),]
        }
        if(type == "Logistic"){
          tab <- table(newdt[[DependentColumn]])
          Continue <- length(tab) == 2L && all(tab > 2L)
        } else {
          if(length(table(newdt[[1]])) > 2 || sum(table(newdt[[1]]) > 4) > 1){
            Continue <- TRUE
          } else {
          Continue <- FALSE
          }
        }
        if(Continue){
          #### generate formula
          nam <- colnames(newdt[, -(ncol(newdt)), with = FALSE])
          if(length(nam) > 1){
            plus <- "+"
            p1 <- paste(nam[-length(nam)], plus, sep = " ")
            p2 <- paste(p1, nam[length(nam)])
            fmla <- as.formula(paste0(colnames(newdt)[ncol(newdt)], " ~ ", p2))
            fmchr <- paste0(colnames(newdt)[ncol(newdt)], " ~ ", p2)
          }else{
            fmla <- as.formula(paste0(colnames(newdt)[ncol(newdt)], " ~ ", nam))
            fmchr <- paste0(colnames(newdt)[ncol(newdt)], " ~ ", nam) }
          #### fit model ####
          glm.fit <- tryCatch({
            if(type == "Logistic"){
              if(!is.na(ReferenceName)){
                newdt[[DependentColumn]] <- relevel(factor(newdt[[DependentColumn]]), ref = ReferenceName)
              }
              suppressMessages(glm(fmla, data = newdt, family = binomial))
            } else if(Parametric){
              suppressMessages(lm(fmla, data = newdt, na.action = na.omit))
            } else {
              npreg::gsm(formula = fmla, data = newdt)
            }
          }, error = function(e) {
            message(sprintf("Predictor '%s' skipped: %s", nam, conditionMessage(e)))
            NULL
          })
          if(is.null(glm.fit)){ next }
          #### create data table
          if(type == "Logistic"){
            prob <- suppressMessages(predict(glm.fit, type = c("response")))
            if(AUC == TRUE){
              ROC <- suppressMessages(roc(glm.fit$y, prob, quiet = TRUE)) }
            dtf <- data.table(coef(summary(glm.fit)))
            if(AUC == TRUE){
              dtf$AUC <- as.numeric(ROC$auc) }
            if(Odds){
              dtf$OddsRatio <- exp(dtf$Estimate)
              dtf$logOddsRatio <- dtf$Estimate
              dtf$OddsRatioP <- dtf$`Pr(>|z|)`
            }
            dtf$names <- gsub("\\(", "", gsub("\\)", "", gsub("s\\(", "", rownames(summary(glm.fit)$coefficients) )) )
            chr <- as.character(glm.fit$formula)
            dtf$formula <- paste(chr[2], chr[1], chr[3:length(chr)])
          }
          if(type == "Linear"){
            if(Parametric){
              dtf <- data.table(coef(summary(glm.fit)))
              dtf$names <- gsub("\\(", "", gsub("\\)", "", gsub("s\\(", "", rownames(summary(glm.fit)$coefficients) )) )
            } else {
              dtf2 <- summary(glm.fit)$s.table
              dtf <- dtf2
              dtf$names <- gsub("\\(", "", gsub("\\)", "", gsub("s\\(", "", rownames(dtf))) )
            }
            chr <- as.character(fmchr)
            dtf$formula <- chr
            ##### check residual normal distribution ####
            residuals <- glm.fit$residuals[!is.na(glm.fit$residuals)]
            if(length(residuals) >= 3L && length(residuals) <= 5000L){
              test <- shapiro.test(residuals)
              dtf$Residualdistribution <- if (test$p.value < 0.05)
                "Not normally distributed" else "normally distributed"
            } else {
              dtf$Residualdistribution <- "Not tested"
            }
          }
          PValColName <- colnames(dtf)[grepl("Pr(>", colnames(dtf), fixed = TRUE)][1L]
          if(is.na(PValColName)){
            message(sprintf("Predictor '%s' skipped: model summary has no p-value column.", nam))
            next
          }
          if(significant == "T"){
            #### subset data table to retain only significant indepenent variables
            dtf <- dtf[dtf[[PValColName]] < 0.05 & (!dtf$names == "(Intercept)"),]
            dtf <- dtf[dtf[[PValColName]] < 0.05 & (!dtf$names == "Intercept"),]
          }
          if(nrow(dtf) > 0){
            dtf$NcountsNoNA <- length(newdt[[nam]][!is.na(newdt[[nam]])])
            dtf <- dtf[(!dtf$names == "(Intercept)"),]
            dtf <- dtf[(!dtf$names == "Intercept"),]
            #### Add confidence intervals ####
            Conf <- tryCatch(
              suppressMessages(stats::confint.default(glm.fit)),
              error = function(e) NULL
            )
            if(!is.null(Conf)){
              Conf <- Conf[!grepl("knot|Intercept", rownames(Conf)), , drop = FALSE]
              if(nrow(Conf) == nrow(dtf)){
                dtf[[colnames(Conf)[[1]]]] <- Conf[,1]
                dtf[[colnames(Conf)[[2]]]] <- Conf[,2]
              }
            }
            compiled_dt <- rbind(compiled_dt, dtf) }
        } else {
          print(paste("There are not enough dependent variable data for", colnames(newdt)[1], "This data was skipped"))
        }
      } else {
        print(paste("there is a limited sample size for", colnames(newdt)[1], "This data was skipped"))
      }
    }
    if(verbose){setTxtProgressBar(pb, m)}
  }
  if(verbose){close(pb)}
  if(!nrow(compiled_dt) == 0){
    PValColName <- colnames(compiled_dt)[grepl("Pr(>", colnames(compiled_dt), fixed = TRUE)][1L]
    compiled_dt <- compiled_dt[(!compiled_dt$names == "(Intercept)"),]
    compiled_dt <- compiled_dt[(!compiled_dt$names == "Intercept"),]
    if(Parametric){
      compiled_dt <- compiled_dt[!duplicated(compiled_dt[,c(1:6)]),]
      compiled_dt$Significance <- "NA"
      compiled_dt[compiled_dt[[PValColName]] < 0.05, Significance := "Significant"]
      compiled_dt[compiled_dt[[PValColName]] >= 0.05, Significance := "Not-Significant"]
      compiled_dt$p.adjust <- p.adjust(compiled_dt[[PValColName]], n = sum(!is.na(compiled_dt[[PValColName]])), method="bonferroni")
      compiled_dt$BonfSignificance <- "Not-Significant"
      compiled_dt[p.adjust < 0.05, BonfSignificance := "Significant"]
      compiled_dt$FDR <- p.adjust(compiled_dt[[PValColName]], n = sum(!is.na(compiled_dt[[PValColName]])), method="BH")
      compiled_dt$BHSignificance <- "Not-Significant"
      compiled_dt[FDR < 0.05, BHSignificance := "Significant"]
    } else {
      if(type == "Linear"){
        compiled_dt$Significance <- "NA"
        compiled_dt[compiled_dt[[PValColName]] < 0.05, Significance := "Significant"]
        compiled_dt[compiled_dt[[PValColName]] >= 0.05, Significance := "Not-Significant"]
        compiled_dt$p.adjust <- p.adjust(compiled_dt[[PValColName]], n = sum(!is.na(compiled_dt[[PValColName]])), method="bonferroni")
        compiled_dt$BonfSignificance <- "Not-Significant"
        compiled_dt[p.adjust < 0.05, BonfSignificance := "Significant"]
        compiled_dt$FDR <- p.adjust(compiled_dt[[PValColName]], n = sum(!is.na(compiled_dt[[PValColName]])), method="BH")
        compiled_dt$BHSignificance <- "Not-Significant"
        compiled_dt[FDR < 0.05, BHSignificance := "Significant"]
      }
      if(type == "Logistic"){
        compiled_dt$Significance <- "NA"
        compiled_dt[compiled_dt[[PValColName]] < 0.05, Significance := "Significant"]
        compiled_dt[compiled_dt[[PValColName]] >= 0.05, Significance := "Not-Significant"]
        compiled_dt$p.adjust <- p.adjust(compiled_dt[[PValColName]], n = sum(!is.na(compiled_dt[[PValColName]])), method="bonferroni")
        compiled_dt$BonfSignificance <- "Not-Significant"
        compiled_dt[p.adjust < 0.05, BonfSignificance := "Significant"]
        compiled_dt$FDR <- p.adjust(compiled_dt[[PValColName]], n = sum(!is.na(compiled_dt[[PValColName]])), method="BH")
        compiled_dt$BHSignificance <- "Not-Significant"
        compiled_dt[FDR < 0.05, BHSignificance := "Significant"]
      }
    }
  }
  return(compiled_dt)
}

#' Compile univariate regressions (linear or logistic) across dependents and independents
#'
#' For every combination of a dependent and an independent variable (with
#' the dependent itself excluded from its own predictor set), fits a
#' univariate linear or logistic model via \code{\link{univariantglmR}} and
#' compiles the results, adding a per-dependent Benjamini-Hochberg FDR
#' column (\code{p_fdr}) and writing one Excel file per dependent to
#' \code{Fpath}. This is the core "regression screen" used throughout the
#' dataMineR workflow.
#'
#' @param DT A data.table containing every column in \code{dependents} and
#'   \code{independents}.
#' @param dependents Character vector of outcome column names to model.
#' @param independents Character vector of candidate predictor column
#'   names.
#' @param Fpath Directory to write one \code{Sub_DataPredictions_*.xlsx}
#'   file per dependent (must already exist).
#' @param ReturnSignificant Logical; if \code{TRUE}, keep only significant
#'   (p < 0.05) results.
#' @param removeAdditionalDependents Logical; exclude every other column in
#'   \code{dependents} from a given dependent's own predictor set.
#' @param Parametric Logical; if \code{FALSE}, drop incomplete cases before
#'   fitting.
#' @param CatTable Optional data.table of binary grouping variables (as
#'   produced by \code{\link{CategoricalVariableCreate}}), used with
#'   \code{CatReferenceCol} to pick each dependent's reference level.
#' @param CatReferenceCol Column name within \code{CatTable} giving the
#'   reference level to record for each dependent (logistic models only).
#' @param type One of \code{"Linear"} or \code{"Logistic"}.
#' @param WriteFiles Logical; write one workbook per dependent to \code{Fpath}.
#' @return A data.table combining the univariate regression results across
#'   every dependent, including a \code{p_fdr} column.
#' @examples
#' \dontrun{
#' data(dfMigrant)
#' RegressionCompile(DT = dfMigrant, dependents = "Age",
#'                    independents = c("COMORBIDITYINDEX", "hpi2score"),
#'                    Fpath = tempdir(), type = "Linear")
#' }
#' @export
RegressionCompile <- function(DT, dependents, independents, Fpath, 
                               ReturnSignificant = TRUE, 
                               removeAdditionalDependents = TRUE, 
                               Parametric = TRUE, 
                               CatTable = NULL,
                               CatReferenceCol = NA,
                               type = "Linear",
                               WriteFiles = TRUE) {
  #### Input validation ####
  if(!type %in% c("Linear", "Logistic")){ stop('type must be "Linear" or "Logistic"') }
  if(WriteFiles && !dir.exists(Fpath)){ stop(paste("Fpath does not exist:", Fpath)) }
  DT <- copy(as.data.table(DT))
  if(!is.null(CatTable)){ CatTable <- copy(as.data.table(CatTable)) }
  missing_dep <- setdiff(dependents, names(DT))
  if(length(missing_dep) > 0){ stop(paste("Dependents not found in DT:", paste(missing_dep, collapse = ", "))) }
  missing_ind <- setdiff(independents, names(DT))
  if(length(missing_ind) > 0){ stop(paste("Independents not found in DT:", paste(missing_ind, collapse = ", "))) }
  #### Sanitize column names on input DT once (track renames for reporting) ####
  original_names <- names(DT)
  clean_names <- make.names(original_names, unique = FALSE)
  if(anyDuplicated(clean_names)){ stop("Column-name sanitization would create duplicate names; rename the input columns first.") }
  if(!identical(original_names, clean_names)){
    renamed <- original_names[original_names != clean_names]
    message(sprintf( "Non-syntactic column names sanitized for model formulas: %s", paste(renamed, collapse = ", ") ))
    setnames(DT, original_names, clean_names)
    name_map <- setNames(clean_names, original_names)
    dependents <- unname(name_map[dependents])
    independents <- unname(name_map[independents])
    #### Adjust names in CatTable ###
    if(!is.null(CatTable)){
      CatTable[["CategoricalVariable"]] <- ifelse(
        CatTable[["CategoricalVariable"]] %in% names(name_map),
        unname(name_map[CatTable[["CategoricalVariable"]]]),
        make.names(CatTable[["CategoricalVariable"]])
      )
    }
  }
  CompiledResults <- data.table()
  for(i in seq_along(dependents)){
    message(sprintf("[%d / %d] Processing dependent: %s", i, length(dependents), dependents[i]))
    #### Subset columns ####
    cols <- if(removeAdditionalDependents){ unique(c(setdiff(independents, dependents), dependents[i]))
    } else { unique(c(independents, dependents)) }
    glmDT <- tryCatch(
      copy(DT[, cols, with = FALSE]),
      error = function(e) {
        message(sprintf("  Skipped '%s': failed to subset columns - %s", dependents[i], conditionMessage(e)))
        NULL
      }
    )
    if(is.null(glmDT)){ next }
    #### Sanitize glmDT column names (catches any remaining edge cases) ####
    setnames(glmDT, names(glmDT), gsub(" ", "_", names(glmDT)))
    setnames(glmDT, names(glmDT), gsub(";", "_", names(glmDT)))
    #### Reorder so dependent is last ####
    glmDT <- glmDT[, c(setdiff(colnames(glmDT), dependents[i]), dependents[i]), with = FALSE]
    k <- seq_len(ncol(glmDT) - 1)
    #### Obtain reference value ####
    CatReference <- NA_character_
    if(!is.null(CatTable)){
      if(!"CategoricalVariable" %in% names(CatTable)){ stop("CatTable must contain a 'CategoricalVariable' column.") }
      if(!is.character(CatReferenceCol) || length(CatReferenceCol) != 1L || !CatReferenceCol %in% names(CatTable)) {
        stop("CatReferenceCol must name a column in CatTable when CatTable is supplied.")
      }
      reference <- CatTable[CategoricalVariable == dependents[i]][[CatReferenceCol]]
      if(length(reference) > 0L){ CatReference <- as.character(reference[1L]) }
    }
    #### Run univariate regression ####
    uni <- tryCatch({
      suppressWarnings(
        univariantglmR(DT = glmDT,
                       key = k,
                       significant = "F",
                       DependentColumn = dependents[i],
                       ReferenceName = CatReference, 
                       type = type,
                       AUC = TRUE,
                       Odds = TRUE,
                       verbose = FALSE,
                       Parametric = Parametric ) )
    }, error = function(e) {
      message(sprintf("  Skipped '%s': univariantglmR error - %s", dependents[i], conditionMessage(e)))
      NULL
    })
    #### Guard: no results returned ####
    if(is.null(uni) || nrow(uni) == 0) {
      message(sprintf("  Skipped '%s': no results produced.", dependents[i]))
      next
    }
    p_col <- colnames(uni)[grepl("Pr(>", colnames(uni), fixed = TRUE)][1L]
    #### Sort by p-value ####
    if (p_col %in% names(uni)) {
      uni <- uni[order(uni[[p_col]], decreasing = FALSE)]
    } else {
      message(sprintf("  Warning: expected p-value column '%s' not found in results for '%s'.", p_col, dependents[i]))
    }
    #### Multiple testing corrections per dependent ####
    if(p_col %in% names(uni)){
      all_p <- uni[[p_col]]
      uni[, p_bonferroni := p.adjust(all_p, method = "bonferroni")]
      uni[, p_fdr := p.adjust(all_p, method = "BH")]
      uni[, sig_unadjusted := all_p < 0.05]
      uni[, sig_bonferroni := p_bonferroni < 0.05]
      uni[, sig_fdr := p_fdr < 0.05]
    }
    #### Parse formula into dependent / independent columns ####
    if("formula" %in% names(uni)){
      spl <- strsplit(uni$formula, " ~ ", fixed = TRUE)
      formula_dt <- rbindlist(lapply(spl, function(x){
        if(length(x) == 2){
          data.table(dependent = x[1], independent = x[2])
        } else {
          data.table(dependent = NA_character_, independent = NA_character_)
        }
      }))
      uni <- cbind(uni, formula_dt)
    } else {
      message(sprintf("  Warning: 'formula' column not found in results for '%s'.", dependents[i]))
    }
    #### Write per-dependent file ####
    if(WriteFiles) {
      out_path <- file.path(Fpath, paste0("Sub_DataPredictions_", dependents[i], ".xlsx"))
      tryCatch(
        openxlsx::write.xlsx(list(RegressionAnalysis = uni), file = out_path, colWidths = "auto"),
        error = function(e) message(sprintf("  Warning: could not write file for '%s' - %s", dependents[i], conditionMessage(e)))
      ) }
    CompiledResults <- rbind(CompiledResults, uni, fill = TRUE)
  }
  #### Guard: nothing compiled ####
  if(nrow(CompiledResults) == 0){
    message("No results compiled across all dependents.")
    return(CompiledResults)
  }
  p_col <- names(CompiledResults)[grepl("Pr(>", names(CompiledResults), fixed = TRUE)][1L]
  if(!is.na(p_col) && p_col %in% names(CompiledResults)){
    CompiledResults[, p_fdr_overall := p.adjust(get(p_col), method = "BH")]
    CompiledResults[, p_bonferroni_overall := p.adjust(get(p_col), method = "bonferroni")]
    CompiledResults[, sig_fdr_overall := p_fdr_overall < 0.05]
    message(sprintf("Family-wide correction across %d tests: %d significant unadjusted, %d significant after overall FDR (BH).",
                    sum(!is.na(CompiledResults[[p_col]])),
                    sum(CompiledResults[[p_col]] < 0.05, na.rm = TRUE),
                    sum(CompiledResults$sig_fdr_overall, na.rm = TRUE)))
  }
  if(ReturnSignificant){
    if(!is.na(p_col) && p_col %in% names(CompiledResults)){
      CompiledResults <- CompiledResults[CompiledResults[[p_col]] < 0.05]
    } else {
      message(sprintf("Warning: could not filter by significance - column '%s' not found.", p_col))
    }
  }
  return(CompiledResults)
}
