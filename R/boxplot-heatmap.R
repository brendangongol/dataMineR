#' Grouped boxplots with pairwise significance brackets
#'
#' For each column in \code{ValueCols}, draws a \pkg{ggpubr} boxplot (with
#' jittered points) grouped by \code{GroupCol}, annotated with pairwise
#' significance-test brackets for every pair of group levels.
#'
#' @param dt A data.table containing \code{GroupCol} and every column in
#'   \code{ValueCols}.
#' @param ValueCols Character vector of numeric column names to plot, one
#'   boxplot per column.
#' @param GroupCol Column name of the grouping variable.
#' @param OutPath Directory to save PNGs into when \code{SavePlots = TRUE}.
#' @param method Test passed to \code{ggpubr::stat_compare_means()}, e.g.
#'   \code{"wilcox.test"} or \code{"t.test"}.
#' @param SavePlots Logical; if \code{TRUE}, write one PNG per column to
#'   \code{OutPath}.
#' @return A named list of ggplot objects, one per element of
#'   \code{ValueCols} with at least two non-missing group levels.
#' @examples
#' \dontrun{
#' load(system.file("extdata", "copd_data.rda", package = "dataMineR"))
#' GroupedBoxplotStatR(copd_data, ValueCols = c("FEV1", "FVC"),
#'                      GroupCol = "group", SavePlots = FALSE)
#' }
#' @export
GroupedBoxplotStatR <- function(dt, ValueCols, GroupCol, OutPath = NULL, method = "wilcox.test", SavePlots = TRUE){
  dt <- as.data.table(dt)
  missing_cols <- setdiff(c(ValueCols, GroupCol), names(dt))
  if(length(missing_cols) > 0L){ stop("Columns not found in dt: ", paste(missing_cols, collapse = ", ")) }
  if(SavePlots && is.null(OutPath)){ stop("OutPath is required when SavePlots = TRUE.") }
  if(SavePlots && !dir.exists(OutPath)){ dir.create(OutPath, recursive = TRUE, showWarnings = FALSE) }
  PlotList <- list()
  for(col in ValueCols){
    PlotDT <- data.table(value = dt[[col]], group = as.character(dt[[GroupCol]]))
    PlotDT <- PlotDT[!is.na(value) & !is.na(group)]
    groups <- sort(unique(PlotDT$group))
    if(nrow(PlotDT) == 0 || length(groups) < 2) next
    comps <- utils::combn(groups, 2, simplify = FALSE)
    p <- ggpubr::ggboxplot(PlotDT, x = "group", y = "value", color = "group", add = "jitter") +
      ggpubr::stat_compare_means(comparisons = comps, method = method) +
      ggplot2::labs(x = GroupCol, y = col, title = col) +
      ggplot2::theme(legend.position = "none")
    PlotList[[col]] <- p
    if(SavePlots && !is.null(OutPath)){
      grDevices::png(file.path(OutPath, paste0(gsub("[^A-Za-z0-9_]", "_", col), "_boxplot.png")), width = 900, height = 700, res = 120)
      print(p); grDevices::dev.off()
    }
  }
  return(PlotList)
}

#' Clustered heatmap of a long-format results table
#'
#' Casts a long-format results table (e.g. the output of
#' \code{\link{RegressionCompile}} or \code{\link{continuous_compare_Compile}})
#' into a row-by-column matrix, optionally -log10-transforming the value
#' column, and draws a clustered \pkg{ggplot2} tile heatmap as an alternative
#' view to \code{\link{ComparisonOverviewBubblePlot}} for dense
#' dependent-by-independent grids.
#'
#' @param ResultsDT A long-format data.table with (at least) the columns
#'   named in \code{RowVar}, \code{ColVar}, and \code{ValueVar}.
#' @param RowVar Column name to use for heatmap rows.
#' @param ColVar Column name to use for heatmap columns.
#' @param ValueVar Column name of the value to plot (typically a p-value or
#'   effect size).
#' @param Transform One of \code{"neglog10"} (plot \code{-log10(ValueVar)},
#'   suited to p-values) or \code{"none"} (plot \code{ValueVar} as-is).
#' @param NAFill Value substituted for missing row/column combinations
#'   before clustering and plotting.
#' @return A list with \code{Matrix} (the row-by-column numeric matrix) and
#'   \code{Plot} (a ggplot object).
#' @examples
#' \dontrun{
#' load(system.file("extdata", "dfMigrant.rda", package = "dataMineR"))
#' CatTable <- CategoricalVariableCreate(dfMigrant, c("Migrant", "Homeless"))
#' res <- RegressionCompile(dfMigrant, dependents = CatTable$CategoricalVariable,
#'                           independents = "Age", Fpath = tempdir(),
#'                           CatTable = CatTable, CatReferenceCol = "Table2Subset",
#'                           type = "Logistic")
#' SignificanceHeatmapR(res, RowVar = "dependent", ColVar = "independent",
#'                       ValueVar = "Pr(>|z|)")
#' }
#' @export
SignificanceHeatmapR <- function(ResultsDT, RowVar, ColVar, ValueVar, Transform = c("neglog10", "none"), NAFill = 0){
  Transform <- match.arg(Transform)
  ResultsDT <- as.data.table(ResultsDT)
  missing_cols <- setdiff(c(RowVar, ColVar, ValueVar), names(ResultsDT))
  if(length(missing_cols) > 0L){ stop("Columns not found in ResultsDT: ", paste(missing_cols, collapse = ", ")) }
  if(!is.numeric(ResultsDT[[ValueVar]])){ stop("ValueVar must name a numeric column.") }
  ResultsDT <- ResultsDT[!is.na(get(RowVar)) & !is.na(get(ColVar)) & !is.na(get(ValueVar))]
  if(nrow(ResultsDT) == 0L){ stop("ResultsDT has no complete row/column/value observations.") }
  if(Transform == "neglog10" && any(ResultsDT[[ValueVar]] < 0, na.rm = TRUE)){
    stop("ValueVar contains negative values and cannot be -log10 transformed.") }
  ResultsDT[, PlotValue := if(Transform == "neglog10"){
    -log10(pmax(get(ValueVar), 1e-300))
  } else {
    get(ValueVar)
  }]
  mat_dt <- data.table::dcast(ResultsDT, stats::as.formula(paste0("`", RowVar, "` ~ `", ColVar, "`")),
                              value.var = "PlotValue", fun.aggregate = mean)
  rn <- mat_dt[[1]]; mat_dt[, (1) := NULL]
  m <- as.matrix(mat_dt); rownames(m) <- rn
  m[is.na(m)] <- NAFill
  row_ord <- tryCatch(rownames(m)[stats::hclust(stats::dist(m))$order], error = function(e) rownames(m))
  col_ord <- tryCatch(colnames(m)[stats::hclust(stats::dist(t(m)))$order], error = function(e) colnames(m))
  long <- data.table::as.data.table(as.table(m))
  colnames(long) <- c("Row", "Col", "Value")
  long[, Row := factor(Row, levels = row_ord)]
  long[, Col := factor(Col, levels = col_ord)]
  p <- ggplot2::ggplot(long, ggplot2::aes(Col, Row, fill = Value)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient(low = "white", high = "firebrick",
                                 name = if(Transform == "neglog10") "-log10(p)" else ValueVar) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(x = ColVar, y = RowVar)
  return(list(Matrix = m, Plot = p))
}
