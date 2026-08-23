#' DeLong's test comparing two ROC curves across one or more datasets
#'
#' For each data.frame in \code{DataList}, fits ROC curves for
#' \code{IndVariable1} and \code{IndVariable2} against \code{DepVariable}
#' (via \pkg{pROC}), computes each curve's AUC and Youden index, and runs
#' DeLong's test comparing the two AUCs.
#'
#' @param DataList A named list of data.frames/data.tables, each containing
#'   \code{DepVariable}, \code{IndVariable1}, and \code{IndVariable2}.
#' @param DepVariable Column name of the (typically binary) outcome.
#' @param IndVariable1 Column name of the first predictor.
#' @param IndVariable2 Column name of the second predictor.
#' @param paired Logical; passed to \code{\link[pROC]{roc.test}} (\code{TRUE}
#'   for two ROC curves from the same subjects/dataset).
#' @param Direction Direction passed to both \code{\link[pROC]{roc}} calls:
#'   \code{"auto"} (default) lets each curve pick its own empirical
#'   direction independently, which can disagree between
#'   \code{IndVariable1}/\code{IndVariable2} when both have weak signal and
#'   makes the DeLong comparison invalid (pROC will warn). Pass \code{"<"}
#'   (controls lower than cases) or \code{">"} explicitly once you know the
#'   expected relationship to force both curves onto a consistent, valid
#'   comparison.
#' @return A named list (matching \code{DataList}), each element containing
#'   \code{ROC1}, \code{ROC2} (the \pkg{pROC} objects), \code{ROCDT} (a
#'   data.table with AUC, Youden index, and DeLong test results for both
#'   predictors), \code{Plots}, and \code{PlotNames}.
#' @examples
#' \dontrun{
#' data(dream_data)
#' ROCDeLongTest(DataList = list(Dream = dream_data), DepVariable = "opioids",
#'                IndVariable1 = "ahi", IndVariable2 = "rdi",
#'                paired = TRUE, Direction = "<")
#' }
#' @export
ROCDeLongTest <- function(DataList, DepVariable, IndVariable1, IndVariable2, paired, Direction = "auto"){
  if(!is.list(DataList) || length(DataList) == 0L){ stop("DataList must be a non-empty list.") }
  if(!Direction %in% c("auto", "<", ">")){ stop("Direction must be 'auto', '<', or '>'.") }
  dataset_names <- names(DataList)
  if(is.null(dataset_names)){ dataset_names <- paste0("Dataset", seq_along(DataList)) }
  dataset_names[is.na(dataset_names) | !nzchar(dataset_names)] <- paste0("Dataset", which(is.na(dataset_names) | !nzchar(dataset_names)))
  build_roc <- function(dataset, predictor, dataset_name) {
    roc_obj <- pROC::roc(dataset[[DepVariable]], dataset[[predictor]], direction = Direction, quiet = TRUE)
    youden <- pROC::coords(roc_obj, x = "best", input = "threshold", best.method = "youden")[1L, , drop = FALSE]
    main <- paste(dataset_name, DepVariable, predictor)
    list(roc = roc_obj,
         row = data.table(dataset = dataset_name, DepVar = DepVariable, IndVar = predictor,
                          Ncontrols = length(roc_obj$controls), Ncases = length(roc_obj$cases),
                          AUC = round(as.numeric(roc_obj$auc), 4),
                          YoudenIndex = youden$specificity + youden$sensitivity - 1),
      plot = local({object <- roc_obj; title <- main; function() plot(object, main = title) }),
      main = main
    )
  }
  ROCList <- lapply(seq_along(DataList), function(i) {
    dataset <- as.data.frame(DataList[[i]])
    required <- c(DepVariable, IndVariable1, IndVariable2)
    missing_cols <- setdiff(required, names(dataset))
    if(length(missing_cols) > 0L){ stop(sprintf("%s is missing: %s", dataset_names[i], paste(missing_cols, collapse = ", "))) }
    if(paired) {
      dataset <- dataset[stats::complete.cases(dataset[, required, drop = FALSE]), , drop = FALSE]
    } else {
      dataset <- dataset[!is.na(dataset[[DepVariable]]), , drop = FALSE]
    }
    if(nrow(dataset) == 0L || length(unique(dataset[[DepVariable]])) != 2L){
      stop(sprintf("%s must contain complete observations from exactly two outcome classes.", dataset_names[i])) }
    first <- build_roc(dataset, IndVariable1, dataset_names[i])
    second <- build_roc(dataset, IndVariable2, dataset_names[i])
    if(!identical(first$roc$direction, second$roc$direction)){
      stop("ROC curves selected different directions; pass Direction = '<' or '>' explicitly.") }
    delong <- pROC::roc.test(first$roc, second$roc, method = "delong", paired = paired)
    delong_row <- data.table(DeLongDepVar = DepVariable, DeLongIndVar1 = IndVariable1,
      DeLongIndVar2 = IndVariable2, DeLongD = unname(delong$statistic),
      DeLongdf = if(is.null(delong$parameter)){ NA_real_ } else { unname(delong$parameter) },
      DeLongmethod = delong$method, DeLongp.value = delong$p.value
    )
    list(ROC1 = first$roc, ROC2 = second$roc,
         ROCDT = rbind(cbind(first$row, delong_row), cbind(second$row, delong_row)),
         Plots = list(Plot1 = first$plot, Plot2 = second$plot),
         PlotNames = list(main1 = first$main, main2 = second$main) )
  })
  names(ROCList) <- dataset_names
  ROCList
}

#' Bland-Altman agreement analysis across one or more datasets
#'
#' For each data.frame in \code{DataList}, computes Bland-Altman agreement
#' statistics (bias, limits of agreement, and their confidence intervals via
#' \pkg{BlandAltmanLeh}/\pkg{blandr}) between \code{Var1} and \code{Var2},
#' with an optional mean-difference plot colored by a third feature.
#'
#' @param DataList A named list of data.frames/data.tables, each containing
#'   \code{Var1} and \code{Var2}.
#' @param Var1 Column name of the first measurement.
#' @param Var2 Column name of the second measurement.
#' @param SigLevel Significance level for the agreement confidence
#'   intervals.
#' @param ColorFeature Logical; if \code{TRUE}, also build a mean-difference
#'   plot colored by \code{ColorBy}.
#' @param ColorBy Column name used to color the mean-difference plot when
#'   \code{ColorFeature = TRUE}.
#' @param ScaleData Logical; scale (but do not center) \code{Var1}/\code{Var2}
#'   before computing agreement statistics.
#' @return A named list (matching \code{DataList} entries that contain both
#'   columns) with \code{BARaw} (raw Bland-Altman statistics), \code{BAPlot}
#'   (the standard agreement plot), and, when \code{ColorFeature = TRUE},
#'   \code{BAPlotFeature} (the colored mean-difference plot).
#' @examples
#' \dontrun{
#' data(peru_data)
#' BlandAltman(DataList = list(Peru = peru_data), Var1 = "hcto_1",
#'             Var2 = "hcto_2", SigLevel = 0.95, ColorFeature = FALSE,
#'             ColorBy = NA, ScaleData = FALSE)
#' }
#' @export
BlandAltman <- function(DataList, Var1, Var2, SigLevel, ColorFeature, ColorBy, ScaleData){
  if(!requireNamespace("blandr", quietly = TRUE)){ stop("Package 'blandr' is required for BlandAltman().") }
  if(ColorFeature && !requireNamespace("BlandAltmanLeh", quietly = TRUE)){ stop("Package 'BlandAltmanLeh' is required when ColorFeature = TRUE.") }
  required <- c(Var1, Var2, if (ColorFeature) ColorBy)
  missing_by_dataset <- lapply(DataList, function(x) setdiff(required, names(x)))
  if(any(lengths(missing_by_dataset) > 0L)){
    bad <- which(lengths(missing_by_dataset) > 0L)[1L]
    stop("A dataset is missing required columns: ", paste(missing_by_dataset[[bad]], collapse = ", "))
  }
  if(!is.numeric(SigLevel) || length(SigLevel) != 1L || SigLevel <= 0 || SigLevel >= 1){ stop("SigLevel must be between 0 and 1.") }
  analysis_data <- lapply(DataList, function(x){
    x <- as.data.frame(x)
    keep <- stats::complete.cases(x[, c(Var1, Var2), drop = FALSE])
    x <- x[keep, , drop = FALSE]
    if(nrow(x) < 2L){ stop("Each dataset requires at least two complete measurement pairs.") }
    if(ScaleData){
      scaled <- scale(as.matrix(x[, c(Var1, Var2), drop = FALSE]),
                      center = FALSE, scale = TRUE)
      if(any(!is.finite(scaled))){ stop("ScaleData cannot be used with a zero-variance measurement.") }
      x[[Var1]] <- scaled[, 1L]
      x[[Var2]] <- scaled[, 2L]
    }
    x
  })
  BARaw <- lapply(analysis_data, function(x) { blandr::blandr.statistics(x[[Var1]], x[[Var2]], sig.level = SigLevel) })
  BAPlot <- lapply(analysis_data, function(x) blandr::blandr.draw(x[[Var1]], x[[Var2]]))
  result <- list(BARaw = BARaw, BAPlot = BAPlot)
  if(ColorFeature){
    result$BAPlotFeature <- lapply(analysis_data, function(x) {
      ba <- BlandAltmanLeh::bland.altman.stats(x[[Var1]], x[[Var2]])
      plot_dt <- data.table(mean = ba$means, diffs = ba$diffs, Feature = x[[ColorBy]])
      ggplot(plot_dt, aes(x = .data[["mean"]], y = .data[["diffs"]], color = Feature)) +
        geom_point() +
        geom_hline(yintercept = ba$lines[c(1L, 3L)], linetype = "dashed", color = "red") +
        geom_hline(yintercept = ba$lines[2L], linetype = "dashed", color = "blue") +
        xlab("Means") + ylab("Differences")
    })
  }
  result
}

#' ROC analysis of a multi-level variable's pairwise level comparisons
#'
#' For each data.frame in \code{DataList} and each row of \code{CompDT} (a
#' pair of levels of \code{SubsetVar}), subsets to just those two levels,
#' recodes \code{SubsetVar} to 0/1, and fits an ROC curve of
#' \code{DepVariable} against that recoded predictor. For example, this permits
#' "Severe vs Mild" as a binary predictor of a clinical outcome.
#'
#' @param DataList A named list of data.frames/data.tables, each containing
#'   \code{DepVariable} and \code{SubsetVar}.
#' @param DepVariable Column name of the (typically binary) outcome.
#' @param SubsetVar Column name of the multi-level variable whose pairwise
#'   levels are compared.
#' @param CompDT A data.frame/data.table with (at least) two columns giving
#'   the pairs of \code{SubsetVar} levels to compare, one pair per row.
#' @return A named list (matching \code{DataList}), each element a named
#'   list (one entry per row of \code{CompDT}) with \code{RawROC},
#'   \code{ROCDT}, \code{ROCPlot}, and \code{Main}.
#' @export
FeatureROCAnalysis <- function(DataList, DepVariable, SubsetVar, CompDT){
  if(!is.list(DataList) || length(DataList) == 0L){ stop("DataList must be a non-empty list.") }
  if(ncol(CompDT) < 2L){ stop("CompDT must contain at least two columns.") }
  dataset_names <- names(DataList)
  if(is.null(dataset_names)){ dataset_names <- paste0("Dataset", seq_along(DataList)) }
  ROCList <- vector("list", length(DataList))
  for(i in seq_along(DataList)){
    dataset <- as.data.frame(DataList[[i]])
    missing_cols <- setdiff(c(DepVariable, SubsetVar), names(dataset))
    if(length(missing_cols) > 0L){ stop(dataset_names[i], " is missing: ", paste(missing_cols, collapse = ", ")) }
    setList <- list()
    comparison_names <- character()
    for(b in seq_len(nrow(CompDT))){
        BWDT <- dataset[dataset[[SubsetVar]] == CompDT[[1]][b] | dataset[[SubsetVar]] == CompDT[[2]][b], , drop = FALSE]
        BWDT[[SubsetVar]] <- as.integer(BWDT[[SubsetVar]] == CompDT[[1]][b])
        BWDT <- BWDT[stats::complete.cases(BWDT[, c(DepVariable, SubsetVar), drop = FALSE]), ]
        if(nrow(BWDT) == 0L || length(unique(BWDT[[DepVariable]])) != 2L){ next }
        ROC1 <- pROC::roc(BWDT[[DepVariable]], BWDT[[SubsetVar]], quiet = TRUE)
        #### keep just the first tied "best" threshold so Youden/ROCDT1 stay single-row ####
        youden <- pROC::coords(ROC1, x="best", input="threshold", best.method="youden")[1, , drop = FALSE]
        Youden <- youden$specificity + youden$sensitivity -1
        ROCDT1 <- data.table(dataset = dataset_names[i], DepVar = DepVariable,
                             IndVar = paste(CompDT[[1]][b], "vs", CompDT[[2]][b], sep = "_"),
                             Ncontrols = length(ROC1$controls), Ncases = length(ROC1$cases),
                             AUC = round(as.numeric(ROC1$auc), 4), YoudenIndex = Youden)
        main1 <- paste(dataset_names[i], DepVariable,
                       paste(CompDT[[1]][b], "vs", CompDT[[2]][b], sep = "_"))
        ROCPlot1 <- local({ ROC1 <- ROC1; m <- main1; function() plot(ROC1, main = m) })
        setList[[length(setList) + 1L]] <- list(RawROC=ROC1, ROCDT = ROCDT1, ROCPlot = ROCPlot1, Main = main1)
        comparison_names <- c(comparison_names, paste(CompDT[[1]][b], "vs", CompDT[[2]][b], sep = "_"))
    }
    names(setList) <- comparison_names
    ROCList[[i]] <- setList
  }
  names(ROCList) <- dataset_names
  return(ROCList)
}

#' AUC as a function of an inclusion-threshold sweep
#'
#' For every unique observed value of \code{Predictor} (used as a candidate
#' inclusion cutoff), subsets \code{dt} to records at or beyond that cutoff,
#' refits an ROC curve of \code{Outcome} on the subset, and records the
#' resulting AUC, producing AUC as a function of how aggressively the
#' cohort is filtered by \code{Predictor}.
#'
#' @param dt A data.table containing \code{Outcome} and \code{Predictor}.
#' @param Outcome Column name of the (typically binary) outcome.
#' @param Predictor Column name of the continuous predictor to sweep.
#' @param Direction One of \code{"ge"} (keep records with \code{Predictor >=}
#'   cutoff, sweeping cutoffs from high to low) or \code{"le"} (keep records
#'   \code{<=} cutoff, sweeping low to high).
#' @param MinN Minimum sub-cohort size required to fit an ROC curve at a
#'   given cutoff; cutoffs producing a smaller or single-class sub-cohort are
#'   skipped.
#' @return A list with \code{Table} (a data.table of \code{cutoff},
#'   \code{AUC}, \code{N}) and \code{Plot} (a ggplot of AUC vs. cutoff).
#' @examples
#' \dontrun{
#' data(dream_data)
#' ThresholdSweepAUCR(dream_data, Outcome = "opioids", Predictor = "ahi",
#'                     Direction = "ge")
#' }
#' @export
ThresholdSweepAUCR <- function(dt, Outcome, Predictor, Direction = c("ge", "le"), MinN = 10){
  Direction <- match.arg(Direction)
  dt <- as.data.table(dt)
  missing_cols <- setdiff(c(Outcome, Predictor), names(dt))
  if(length(missing_cols) > 0L){ stop("Columns not found in dt: ", paste(missing_cols, collapse = ", ")) }
  if(!is.numeric(MinN) || length(MinN) != 1L || MinN < 2L){ stop("MinN must be at least 2.") }
  dt <- dt[!is.na(get(Outcome)) & !is.na(get(Predictor))]
  cuts <- sort(unique(dt[[Predictor]]), decreasing = (Direction == "ge"))
  out <- rbindlist(lapply(cuts, function(cut){
    sub <- if(Direction == "ge"){ dt[get(Predictor) >= cut]} else {dt[get(Predictor) <= cut]}
    if(nrow(sub) < MinN || length(unique(sub[[Outcome]])) < 2){ return(NULL) }
    r <- tryCatch(pROC::roc(sub[[Outcome]], sub[[Predictor]], quiet = TRUE), error = function(e) NULL)
    if(is.null(r)){ return(NULL) }
    data.table(cutoff = cut, AUC = round(as.numeric(r$auc), 4), N = nrow(sub))
  }), fill = TRUE)
  if(nrow(out) == 0L){ stop("No threshold produced a sub-cohort with enough rows and two outcome classes.") }
  p <- ggplot2::ggplot(out, ggplot2::aes(cutoff, AUC)) +
    ggplot2::geom_line() + ggplot2::geom_point(ggplot2::aes(size = N)) +
    ggplot2::labs(x = paste0(Predictor, " cutoff (", Direction, ")"), y = "AUC") +
    ggpubr::theme_pubr()
  return(list(Table = out, Plot = p))
}
