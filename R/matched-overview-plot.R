#' Apply a validated numeric comparison to a results table
#'
#' @param DT Data frame or data table.
#' @param column Numeric column to filter.
#' @param rule Comparison such as \code{"<0.05"}, or \code{NA}.
#' @param argument Argument name used in errors.
#' @return The filtered table.
#' @keywords internal
#' @noRd
Matching_ApplySubset <- function(DT, column, rule, argument) {
  if (length(rule) == 1L && is.na(rule)) return(DT)
  if (!is.character(rule) || length(rule) != 1L || is.na(rule))
    stop(argument, " must be NA or one comparison such as '<0.05'.", call. = FALSE)
  hit <- regexec("^\\s*(<=|>=|==|!=|<|>)\\s*(-?(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)(?:[eE][+-]?[0-9]+)?)\\s*$",
                 rule, perl = TRUE)
  pieces <- regmatches(rule, hit)[[1]]
  if (length(pieces) != 3L)
    stop(argument, " must contain one comparison operator and a numeric threshold, e.g. '<0.05'.", call. = FALSE)
  threshold <- as.numeric(pieces[3])
  if (!is.finite(threshold)) stop(argument, " threshold must be finite.", call. = FALSE)
  value <- DT[[column]]
  keep <- switch(pieces[2], "<" = value < threshold, "<=" = value <= threshold,
                 ">" = value > threshold, ">=" = value >= threshold,
                 "==" = value == threshold, "!=" = value != threshold)
  rows <- which(!is.na(keep) & keep)
  if (data.table::is.data.table(DT)) DT[rows] else DT[rows, , drop = FALSE]
}

#' Overview bubble plot of matching and weighting results
#'
#' The layout of \code{\link{ComparisonOverviewBubblePlot}} (outcomes across,
#' exposures down, bubble size for significance, colour for the effect, marginal
#' sample-size bars) for the folders written by \code{PSRegression()},
#' \code{IPTWRegression()} and \code{EBRegression()}. Only the
#' \code{*_Sub_DataPredictions_*.xlsx} files are read, from their
#' RegressionAnalysis sheet. Each fitted model gets its own row: when an
#' exposure x outcome pair was fitted more than once, the \code{LabelColumns}
#' that differ are added to its label. \code{"Pr(>|z|)"} is read as
#' \code{"Pr(>|t|)"}. The marginal bars come from \code{NColumn} in the
#' results.
#'
#' @param Fpath Folder of results.
#' @param SizeVariable Column mapped to bubble size (typically the p-value).
#' @param log2SizeVariable Logical; plot \code{-log2} of \code{SizeVariable}.
#' @param SizeSubsetLevel Optional filter on \code{SizeVariable}, for example
#'   \code{"<0.05"}.
#' @param ColorVariable Column mapped to colour (typically \code{"Estimate"}).
#' @param log2ColorVariable Logical; plot \code{-log2} of \code{ColorVariable}.
#' @param ColorSubsetLevel Optional filter on \code{ColorVariable}.
#' @param RemoveAllEmpty Logical; drop rows and columns with no results.
#' @param XColumn Column for the x axis.
#' @param YColumn Column for the y axis.
#' @param LabelColumns Columns that can tell results in one cell apart; only
#'   those that differ are added to the row label.
#' @param LabelSeparator Text between the parts of a row label.
#' @param SignificantDigits Digits in the overview table.
#' @param NColumn Results column drawn in the marginal bars.
#' @param ShowMarginalBars Logical; draw the marginal bars.
#' @param NAFill Colour of cells without a result.
#' @param ColorLimitQuantile Quantile of the absolute colour values that caps
#'   the colour scale.
#' @param ClusterMethod \code{"correlation"}, \code{"signed"},
#'   \code{"tanimoto"} or \code{"none"}.
#' @param ClusterAxes \code{"both"}, \code{"x"}, \code{"y"} or \code{"none"}.
#' @param MarginalAxes \code{"both"}, \code{"y"}, \code{"x"} or \code{"none"}.
#' @param MarginalPanelRatio Main panel size relative to a marginal panel.
#' @param SigColumn Optional column; bubbles below \code{SigThreshold} are
#'   outlined.
#' @param SigThreshold Threshold for \code{SigColumn}.
#' @param Xlab,Ylab Optional axis titles.
#' @param RowLabelStyle \code{"compact"} (short analysis IDs with the exposure)
#'   or \code{"full"} (the full row labels).
#' @param LabelWidth Maximum characters of a compact label.
#' @param MaxAxisLabels Maximum number of y-axis labels shown.
#' @return A list: \code{Plot}, \code{BubblePlot}, \code{Table},
#'   \code{ClusterOrder}, \code{MarginalN}, \code{ColorLimits} and
#'   \code{Labels} (each row label and the columns it was built from).
#' @examples
#' \dontrun{
#' bp <- MatchedComparisonOverviewBubblePlot(Fpath = file.path(tempdir(), "PS"),
#'   SizeSubsetLevel = "<0.05", SigColumn = "p_fdr")
#' bp$Plot
#' }
#' @export
MatchedComparisonOverviewBubblePlot <- function(Fpath, SizeVariable = "Pr(>|t|)", log2SizeVariable = TRUE, SizeSubsetLevel = NA,
                                                ColorVariable = "Estimate", log2ColorVariable = FALSE, ColorSubsetLevel = NA,
                                                RemoveAllEmpty = TRUE, XColumn = "dependent", YColumn = "independent",
                                                LabelColumns = c("MatchingVariable", "Method", "Type", "estimand", "trim", "moments",
                                                                 "IndependentVariable", "Ratio", "PSCaliperSD", "CovariateCalipers"),
                                                LabelSeparator = " | ", SignificantDigits = 1, NColumn = "NcountsNoNA",
                                                ShowMarginalBars = TRUE, NAFill = "grey85", ColorLimitQuantile = 0.95,
                                                ClusterMethod = c("correlation", "signed", "tanimoto", "none"),
                                                ClusterAxes  = c("both", "x", "y", "none"),
                                                MarginalAxes = c("both", "y", "x", "none"),
                                                MarginalPanelRatio = 8,
                                                SigColumn = NULL, SigThreshold = 0.05,
                                                Xlab = NULL, Ylab = NULL,
                                                RowLabelStyle = c("compact", "full"), LabelWidth = 24, MaxAxisLabels = 30){

  RowLabelStyle <- match.arg(RowLabelStyle)
  if(length(MaxAxisLabels)!=1L || !is.numeric(MaxAxisLabels) || !is.finite(MaxAxisLabels) || MaxAxisLabels<2)
    stop("MaxAxisLabels must be at least 2.")
  if(length(LabelWidth)!=1L || !is.numeric(LabelWidth) || !is.finite(LabelWidth) || LabelWidth<8)
    stop("LabelWidth must be a finite number of at least 8 characters.")
  ClusterMethod <- match.arg(ClusterMethod)
  ClusterAxes   <- match.arg(ClusterAxes)
  MarginalAxes  <- match.arg(MarginalAxes)
  if(!is.numeric(MarginalPanelRatio) || length(MarginalPanelRatio) != 1 ||
     !is.finite(MarginalPanelRatio) || MarginalPanelRatio <= 0)
    stop("MarginalPanelRatio must be a single positive number.")
  p_alias       <- function(v) if(!is.null(v) && identical(v, "Pr(>|z|)")) "Pr(>|t|)" else v
  SizeVariable  <- p_alias(SizeVariable); ColorVariable <- p_alias(ColorVariable); SigColumn <- p_alias(SigColumn)

  files <- list.files(Fpath, pattern = "_Sub_DataPredictions_.*\\.xlsx$", ignore.case = TRUE, full.names = TRUE)
  if(length(files) == 0) stop(sprintf("No <prefix>_Sub_DataPredictions_<outcome>.xlsx files found in Fpath: %s", Fpath))
  PlotDT <- rbindlist(lapply(files, function(f)
    as.data.table(openxlsx::read.xlsx(f, sheet = "RegressionAnalysis", check.names = FALSE, sep.names = " "))), fill = TRUE)
  if(nrow(PlotDT) == 0) stop("No rows read from the result files in Fpath.")
  missing_cols <- setdiff(c(XColumn, YColumn, SizeVariable, ColorVariable, SigColumn), names(PlotDT))
  if(length(missing_cols))
    stop(sprintf("Column(s) not in the RegressionAnalysis sheets: %s. Numeric columns available: %s",
                 paste(missing_cols, collapse = ", "),
                 paste(names(PlotDT)[vapply(PlotDT, is.numeric, logical(1))], collapse = ", ")))

  LabelColumns <- intersect(LabelColumns, names(PlotDT))
  cell   <- paste(PlotDT[[XColumn]], PlotDT[[YColumn]], sep = "\r")
  varies <- LabelColumns[vapply(LabelColumns, function(cn)
    any(tapply(as.character(PlotDT[[cn]]), cell, function(v) length(unique(v)) > 1)), logical(1))]
  PlotDT[, .row_label := do.call(paste, c(unname(lapply(.SD, as.character)), sep = LabelSeparator)), .SDcols = c(YColumn, varies)]
  dup <- duplicated(PlotDT[, c(XColumn, ".row_label"), with = FALSE])
  if(any(dup)) {
    # Example rows can intentionally repeat a model. Preserve each analysis using
    # its stable setup ID when the requested descriptive fields are insufficient.
    if (!"ModelID" %in% names(PlotDT))
      stop("Results have duplicate labels and no ModelID; include a distinguishing field in LabelColumns.")
    cells <- PlotDT[, c(XColumn, ".row_label"), with=FALSE]
    collided <- duplicated(cells) | duplicated(cells, fromLast=TRUE)
    ambiguous_labels <- unique(PlotDT$.row_label[collided])
    needs_id <- PlotDT$.row_label %in% ambiguous_labels
    ids <- as.character(PlotDT$ModelID)
    if(any(is.na(ids[needs_id]) | !nzchar(trimws(ids[needs_id]))))
      stop("Duplicate plot labels require nonblank ModelID values.")
    PlotDT[needs_id, .row_label := paste0(.row_label, LabelSeparator, "ModelID=", ModelID)]
    if(anyDuplicated(PlotDT[, c(XColumn, ".row_label"), with=FALSE]))
      stop("Results repeat the same outcome and ModelID. Check for duplicated result files or setup IDs.")
    varies <- unique(c(varies, "ModelID"))
    message("Added ModelID to ambiguous row labels; all analysis rows retained.")
  }
  Labels  <- unique(PlotDT[, c(".row_label", YColumn, varies), with = FALSE])
  setnames(Labels, ".row_label", "label")
  Labels[, analysis_id := sprintf("R%03d", seq_len(.N))]
  short_name <- as.character(Labels[[YColumn]])
  short_name <- ifelse(nchar(short_name)>LabelWidth, paste0(substr(short_name,1,LabelWidth-3),"..."),short_name)
  Labels[, display_label := if(RowLabelStyle=="compact") paste(analysis_id,short_name) else label]
  display_labels <- setNames(Labels$display_label, Labels$label)
  YLabCol <- ".row_label"
  PlotDT_full <- copy(PlotDT)

  fmla    <- as.formula(paste0("`", YLabCol, "` ~ `", XColumn, "`"))
  to_mat  <- function(vv){
    w  <- dcast(PlotDT, fmla, value.var = vv)
    rn <- as.character(w[[YLabCol]])
    m  <- as.matrix(w[, !YLabCol, with = FALSE])
    rownames(m) <- rn
    m
  }
  SizeMat  <- to_mat(SizeVariable)
  ColorMat <- to_mat(ColorVariable)[rownames(SizeMat), colnames(SizeMat), drop = FALSE]

  CombinedDF <- matrix(paste(round(SizeMat, SignificantDigits),
                             paste0("(", round(ColorMat, SignificantDigits), ")")),
                       nrow = nrow(SizeMat), dimnames = dimnames(SizeMat))
  CombinedDF <- as.data.frame(CombinedDF, stringsAsFactors = FALSE)
  empty_mat  <- is.na(SizeMat) | is.na(ColorMat)
  if(RemoveAllEmpty){
    keep_rows <- rowSums(empty_mat) != ncol(empty_mat)
    keep_cols <- colSums(empty_mat) != nrow(empty_mat)
    CombinedDF <- CombinedDF[keep_rows, keep_cols, drop = FALSE]
  }
  CombinedDF <- cbind(data.frame(label = rownames(CombinedDF), stringsAsFactors = FALSE), CombinedDF)
  rownames(CombinedDF) <- NULL

  cluster_order <- function(mat, method){
    if(method == "none" || nrow(mat) < 3) return(rownames(mat))
    m <- mat; m[!is.finite(m)] <- 0                       # neutral fill for distance only
    d <- switch(method,
      "correlation" = { cc <- suppressWarnings(cor(t(m))); cc[!is.finite(cc)] <- 0; as.dist(1 - cc) },
      "signed"      = dist(sign(m), method = "manhattan"),
      "tanimoto"    = proxy::dist(m != 0, method = "Jaccard"))
    if(anyNA(d)) d[is.na(d)] <- max(d, na.rm = TRUE)
    rownames(m)[stats::hclust(d, method = "average")$order]
  }
  y_order <- if(ClusterAxes %in% c("both","y")) cluster_order(ColorMat, ClusterMethod) else rownames(ColorMat)
  x_order <- if(ClusterAxes %in% c("both","x")) cluster_order(t(ColorMat), ClusterMethod) else colnames(ColorMat)

  MarginalN <- NULL
  if(ShowMarginalBars){
    ndt <- if(NColumn %in% names(PlotDT_full)) PlotDT_full[!is.na(get(NColumn))] else PlotDT_full[0]
    if(nrow(ndt) == 0){
      message(sprintf("No '%s' values in the results -- marginal bars skipped.", NColumn))
      ShowMarginalBars <- FALSE
    } else {
      agg <- function(grp) ndt[, .(medN = as.numeric(median(get(NColumn))),
                                   minN = as.numeric(min(get(NColumn))),
                                   maxN = as.numeric(max(get(NColumn)))), by = grp]
      MarginalN <- list(x = agg(XColumn), y = agg(YLabCol))
    }
  }

  PlotDT <- Matching_ApplySubset(PlotDT, ColorVariable, ColorSubsetLevel, "ColorSubsetLevel")
  PlotDT <- Matching_ApplySubset(PlotDT, SizeVariable, SizeSubsetLevel, "SizeSubsetLevel")
  if(nrow(PlotDT) == 0) stop("No rows remain after subsetting -- check SizeSubsetLevel / ColorSubsetLevel.")
  if(log2ColorVariable) PlotDT[[ColorVariable]] <- -log2(PlotDT[[ColorVariable]])
  if(log2SizeVariable)  PlotDT[[SizeVariable]]  <- -log2(PlotDT[[SizeVariable]])
  if(!is.null(SigColumn))
    set(PlotDT, j = ".SigFlag", value = PlotDT[[SigColumn]] < SigThreshold)

  PlotDT <- PlotDT[!is.na(get(SizeVariable)) & !is.na(get(ColorVariable))]
  if(nrow(PlotDT) == 0) stop("No plottable rows remain after removing NA size/color values.")
  x_lv <- x_order[x_order %in% unique(PlotDT[[XColumn]])]
  y_lv <- y_order[y_order %in% unique(PlotDT[[YLabCol]])]
  PlotDT[[XColumn]] <- factor(PlotDT[[XColumn]], levels = x_lv)
  PlotDT[[YLabCol]] <- factor(PlotDT[[YLabCol]], levels = y_lv)

  ColorLabel <- if(log2ColorVariable) paste0("-log2(", ColorVariable, ")") else ColorVariable
  SizeLabel  <- if(log2SizeVariable)  paste0("-log2(", SizeVariable,  ")") else SizeVariable
  if(is.null(Xlab)) Xlab <- gsub("_", " ", XColumn)
  if(is.null(Ylab)) Ylab <- if(RowLabelStyle=="compact") "Analysis" else gsub("_", " ", YColumn)

  cvals <- PlotDT[[ColorVariable]]; cvals <- cvals[is.finite(cvals)]
  lim   <- stats::quantile(abs(cvals), ColorLimitQuantile, na.rm = TRUE)
  if(!is.finite(lim) || lim == 0) lim <- max(abs(cvals), 1, na.rm = TRUE)
  color_limits <- if(length(cvals) && min(cvals) < 0) c(-lim, lim) else c(0, lim)

  shown_empty <- empty_mat[y_lv, x_lv, drop = FALSE]
  if(length(shown_empty) && any(shown_empty)){
    idx <- which(shown_empty, arr.ind = TRUE)
    grid_bg <- data.table(
      X = factor(colnames(shown_empty)[idx[, "col"]], levels = x_lv),
      Y = factor(rownames(shown_empty)[idx[, "row"]], levels = y_lv)
    )
  } else {
    grid_bg <- data.table(X = factor(character(0), levels = x_lv),
                          Y = factor(character(0), levels = y_lv))
  }

  P2 <- ggplot(PlotDT, aes(x = .data[[XColumn]], y = .data[[YLabCol]])) +
    geom_point(data = grid_bg, aes(x = X, y = Y), shape = 21,
               size = 4, fill = NAFill, colour = NAFill, alpha = 0.65,
               inherit.aes = FALSE, show.legend = FALSE) +
    geom_point(aes(size = .data[[SizeVariable]], color = .data[[ColorVariable]])) +
    scale_color_distiller(palette = "RdBu", direction = 1,
                          limits = color_limits, oob = scales::squish) +
    scale_size_continuous(range = c(1, 10)) +
    scale_x_discrete(limits = x_lv, drop = FALSE) +
    scale_y_discrete(limits = y_lv, labels = function(x) unname(display_labels[x]),
                     breaks = y_lv[unique(round(seq(1,length(y_lv),length.out=min(length(y_lv),MaxAxisLabels))))], drop = FALSE) +
    theme_minimal() +
    labs(x = Xlab, y = Ylab, color = ColorLabel, size = SizeLabel) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size=8))
  if(!is.null(SigColumn))
    P2 <- P2 + geom_point(data = PlotDT[.SigFlag == TRUE],
                          aes(size = .data[[SizeVariable]]),
                          shape = 21, fill = NA, colour = "black", stroke = 1.1,
                          show.legend = FALSE)

  FinalPlot <- P2
  if(ShowMarginalBars && requireNamespace("patchwork", quietly = TRUE)){
    xN <- MarginalN$x[get(XColumn) %in% x_lv]; xN[[XColumn]] <- factor(xN[[XColumn]], levels = x_lv)
    yN <- MarginalN$y[get(YLabCol) %in% y_lv]; yN[[YLabCol]] <- factor(yN[[YLabCol]], levels = y_lv)
    show_top   <- MarginalAxes %in% c("both", "x") && nrow(xN) > 0
    show_right <- MarginalAxes %in% c("both", "y") && nrow(yN) > 0
    if(show_top){
      top_bar <- ggplot(xN, aes(x = .data[[XColumn]], y = medN)) +
        geom_col(fill = "grey55", width = 0.45) +
        geom_errorbar(aes(ymin = minN, ymax = maxN), width = 0.3, colour = "grey30") +
        scale_x_discrete(limits = x_lv, drop = FALSE) +
        labs(y = "Analysis N") + theme_minimal() +
        theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
              panel.grid.major.x = element_blank())
    }
    if(show_right){
      right_bar <- ggplot(yN, aes(y = .data[[YLabCol]], x = medN)) +
        geom_col(fill = "grey55", width = 0.45) +
        geom_errorbar(aes(xmin = minN, xmax = maxN), width = 0.3, colour = "grey30") +
        scale_y_discrete(limits = y_lv, drop = FALSE) +
        labs(x = "Analysis N") + theme_minimal() +
        theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
              panel.grid.major.y = element_blank())
    }
    if(show_top && show_right){
      FinalPlot <- (top_bar + patchwork::plot_spacer() + P2 + right_bar) +
        patchwork::plot_layout(ncol = 2,
                               widths  = c(MarginalPanelRatio, 1),
                               heights = c(1, MarginalPanelRatio),
                               guides = "collect") &
        ggplot2::theme(legend.position = "right")
    } else if(show_top){
      FinalPlot <- (top_bar / P2) +
        patchwork::plot_layout(heights = c(1, MarginalPanelRatio), guides = "collect") &
        ggplot2::theme(legend.position = "right")
    } else if(show_right){
      FinalPlot <- (P2 + right_bar) +
        patchwork::plot_layout(widths = c(MarginalPanelRatio, 1), guides = "collect") &
        ggplot2::theme(legend.position = "right")
    }
  }

  return(list(Plot = FinalPlot, BubblePlot = P2, Table = CombinedDF,
              ClusterOrder = list(x = x_lv, y = y_lv), MarginalN = MarginalN,
              ColorLimits = color_limits, Labels = as.data.frame(Labels)))
}
