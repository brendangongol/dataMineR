#' Clustered bubble-plot overview of many compiled comparisons/regressions
#'
#' Reads every Excel file in \code{Fpath} (as written by
#' \code{\link{continuous_compare_Compile}}, \code{\link{Categorical_Compare_compile}},
#' or \code{\link{RegressionCompile}}), compiles them into one long table, and
#' draws a clustered bubble-plot grid (dependent x independent, or comparison
#' x parameter) where bubble size encodes a p-value-like statistic and color
#' encodes an effect size. Optional marginal bar panels show per-row/column
#' sample sizes, and results can be outlined by a significance column. Supply
#' \code{ResultsDT} to plot an already-compiled result without rereading files.
#'
#' @param Fpath Directory of Excel files to read and compile. May be
#'   \code{NULL} when \code{ResultsDT} is supplied.
#' @param SizeVariable Column name whose (optionally log2-transformed) value
#'   maps to bubble size (typically a p-value).
#' @param log2SizeVariable Logical; log2-transform \code{SizeVariable} before
#'   mapping to size.
#' @param SizeSubsetLevel Optional filter string (e.g. \code{"<0.05"}) applied
#'   to \code{SizeVariable} before plotting.
#' @param ColorVariable Column name whose (optionally log2-transformed) value
#'   maps to bubble color (typically an effect size or estimate).
#' @param log2ColorVariable Logical; log2-transform \code{ColorVariable}
#'   before mapping to color.
#' @param ColorSubsetLevel Optional filter string applied to
#'   \code{ColorVariable} before plotting.
#' @param RemoveAllEmpty Logical; drop rows/columns that are entirely empty
#'   after filtering.
#' @param XColumn Column name to use for the plot's x axis (e.g.
#'   \code{"comparison"} or \code{"dependent"}).
#' @param YColumn Column name to use for the plot's y axis (e.g.
#'   \code{"parameter"} or \code{"independent"}).
#' @param SignificantDigits Number of digits used when formatting values in
#'   the compiled overview table.
#' @param SourceData Optional data.frame/data.table used to compute marginal
#'   sample-size bars (complete-case N per row/column).
#' @param ShowMarginalBars Logical; draw marginal N-bar panels alongside the
#'   main bubble grid (requires \code{SourceData}).
#' @param NAFill Fill color used for missing/empty grid cells.
#' @param ColorLimitQuantile Quantile used to cap the color scale, reducing
#'   the influence of extreme outliers.
#' @param ClusterMethod One of \code{"correlation"}, \code{"signed"},
#'   \code{"tanimoto"}, or \code{"none"} controlling how rows/columns are
#'   reordered by similarity.
#' @param ClusterAxes One of \code{"both"}, \code{"x"}, \code{"y"}, or
#'   \code{"none"} controlling which axes are reordered by clustering.
#' @param MarginalAxes One of \code{"both"}, \code{"y"}, \code{"x"}, or
#'   \code{"none"} controlling which marginal N-bar panels are drawn.
#' @param MarginalPanelRatio Positive number giving the main-panel size
#'   relative to a marginal bar panel (higher = thinner bars).
#' @param SigColumn Optional column name of a logical/threshold column used
#'   to outline bubbles that meet \code{SigThreshold}.
#' @param SigThreshold Threshold applied to \code{SigColumn} for outlining.
#' @param Xlab Optional x-axis label override.
#' @param Ylab Optional y-axis label override.
#' @param ResultsDT Optional precompiled data.frame/data.table. When supplied,
#'   it is used instead of reading Excel files from \code{Fpath}.
#' @return A list with (at least) \code{Plot} (the ggplot object) and
#'   \code{Table} (the compiled overview data.table).
#' @examples
#' \dontrun{
#' ComparisonOverviewBubblePlot(Fpath = "ContinuousComparisonResults",
#'   XColumn = "comparison", YColumn = "parameter",
#'   SizeVariable = "mann_whitney_u_pvalue", ColorVariable = "CliffsDelta")
#' }
#' @export
ComparisonOverviewBubblePlot <- function(Fpath = NULL, SizeVariable, log2SizeVariable = FALSE, SizeSubsetLevel = NA,
                                         ColorVariable, log2ColorVariable = FALSE, ColorSubsetLevel = NA, RemoveAllEmpty = TRUE,
                                         XColumn = "comparison", YColumn = "parameter", SignificantDigits = 1,
                                         SourceData = NULL, ShowMarginalBars = NULL, NAFill = "grey85",
                                         ColorLimitQuantile = 0.95,
                                         ClusterMethod = c("correlation", "signed", "tanimoto", "none"),
                                         ClusterAxes  = c("both", "x", "y", "none"),
                                         MarginalAxes = c("both", "y", "x", "none"),
                                         MarginalPanelRatio = 8,
                                         SigColumn = NULL, SigThreshold = 0.05,
                                         Xlab = NULL, Ylab = NULL, ResultsDT = NULL){
  ClusterMethod <- match.arg(ClusterMethod)
  ClusterAxes <- match.arg(ClusterAxes)
  MarginalAxes <- match.arg(MarginalAxes)
  if(!is.numeric(MarginalPanelRatio) || length(MarginalPanelRatio) != 1 || !is.finite(MarginalPanelRatio) || MarginalPanelRatio <= 0){
    stop("MarginalPanelRatio must be a single positive number.") }
  if (!is.numeric(ColorLimitQuantile) || length(ColorLimitQuantile) != 1L || !is.finite(ColorLimitQuantile) || ColorLimitQuantile <= 0 || ColorLimitQuantile > 1){
    stop("ColorLimitQuantile must be greater than 0 and at most 1.") }
  if (!is.null(ShowMarginalBars) && (!is.logical(ShowMarginalBars) || length(ShowMarginalBars) != 1L || is.na(ShowMarginalBars))){
    stop("ShowMarginalBars must be NULL, TRUE, or FALSE.") }
  sanitize <- function(v) gsub(";", "_", gsub(" ", "_", as.character(v)))
  name_key <- function(v){ v <- gsub("`", "", sanitize(v))
                           v <- make.names(v, unique = FALSE)
                           tolower(gsub("[^[:alnum:]]+", "_", v))
                          }
  read_cols <- c(XColumn, YColumn, SizeVariable, ColorVariable)
  if(!is.null(ResultsDT)){
    PlotDT <- as.data.table(copy(ResultsDT))
    missing_cols <- setdiff(read_cols, names(PlotDT))
    if(length(missing_cols) > 0L){
      stop("ResultsDT is missing required columns: ", paste(missing_cols, collapse = ", "))
    }
    keep <- intersect(c(read_cols, SigColumn), names(PlotDT))
    PlotDT <- PlotDT[, keep, with = FALSE]
  } else {
    if(is.null(Fpath) || length(Fpath) != 1L || !dir.exists(Fpath)){
      stop("Fpath must name an existing directory when ResultsDT is NULL.")
    }
    files <- list.files(Fpath, pattern = "\\.xlsx$|\\.xlsm$", ignore.case = TRUE, full.names = TRUE)
    if(length(files) == 0){ stop(sprintf("No .xlsx/.xlsm result files found in Fpath: %s", Fpath)) }
    PlotDT <- data.table()
    for(i in seq_along(files)){
      temp <- as.data.table(openxlsx::read.xlsx(files[i], sheet = 1))
      keep <- intersect(c(read_cols, SigColumn), names(temp))
      if(!all(read_cols %in% keep)) next
      PlotDT <- rbind(PlotDT, temp[, keep, with = FALSE], fill = TRUE)
    }
  }
  if(nrow(PlotDT) == 0){ stop("No usable rows read from result files in Fpath.") }
  if (!is.numeric(PlotDT[[SizeVariable]]) || !is.numeric(PlotDT[[ColorVariable]])){
    stop("SizeVariable and ColorVariable must be numeric columns in the result files.") }
  PlotDT_full <- copy(PlotDT)
  fmla <- as.formula(paste0("`", YColumn, "` ~ `", XColumn, "`"))
  to_mat <- function(vv){
    w <- dcast(PlotDT, fmla, value.var = vv, fun.aggregate = mean, fill = NA_real_, na.rm = TRUE)
    rn <- as.character(w[[YColumn]])
    m <- as.matrix(w[, !YColumn, with = FALSE])
    rownames(m) <- rn
    m
  }
  SizeMat <- to_mat(SizeVariable)
  ColorMat <- to_mat(ColorVariable)[rownames(SizeMat), colnames(SizeMat), drop = FALSE]
  CombinedDF <- matrix(paste(round(SizeMat, SignificantDigits), paste0("(", round(ColorMat, SignificantDigits), ")")),
                       nrow = nrow(SizeMat), dimnames = dimnames(SizeMat))
  CombinedDF <- as.data.frame(CombinedDF, stringsAsFactors = FALSE)
  empty_mat <- is.na(SizeMat) | is.na(ColorMat)
  if(RemoveAllEmpty){
    keep_rows <- rowSums(empty_mat) != ncol(empty_mat)
    keep_cols <- colSums(empty_mat) != nrow(empty_mat)
    CombinedDF <- CombinedDF[keep_rows, keep_cols, drop = FALSE]
  }
  CombinedDF <- cbind(setNames(data.frame(rownames(CombinedDF), stringsAsFactors = FALSE), YColumn), CombinedDF)
  rownames(CombinedDF) <- NULL
  cluster_order <- function(mat, method){
    if(method == "none" || nrow(mat) < 3){ return(rownames(mat)) }
    if (method == "tanimoto" && !requireNamespace("proxy", quietly = TRUE)){
      stop("Package 'proxy' is required for ClusterMethod = 'tanimoto'.") }
    m <- mat; m[!is.finite(m)] <- 0 
    d <- switch(method,
      "correlation" = { cc <- suppressWarnings(cor(t(m))); cc[!is.finite(cc)] <- 0; as.dist(1 - cc) },
      "signed" = dist(sign(m), method = "manhattan"),
      "tanimoto" = proxy::dist(m != 0, method = "Jaccard"))
    if(anyNA(d)) d[is.na(d)] <- max(d, na.rm = TRUE)
    rownames(m)[stats::hclust(d, method = "average")$order]
  }
  y_order <- if(ClusterAxes %in% c("both","y")){ cluster_order(ColorMat, ClusterMethod) } else { rownames(ColorMat) }
  x_order <- if(ClusterAxes %in% c("both","x")){ cluster_order(t(ColorMat), ClusterMethod) } else { colnames(ColorMat) }
  MarginalN <- NULL
  if(is.null(ShowMarginalBars)){ ShowMarginalBars <- !is.null(SourceData) }
  if(ShowMarginalBars && !is.null(SourceData)){
    sd <- as.data.table(copy(SourceData))
    source_lookup <- setNames(names(sd), name_key(names(sd)))
    source_lookup <- source_lookup[!duplicated(names(source_lookup))]
    resolve_source_col <- function(v) {
      key <- name_key(v)
      if(key %in% names(source_lookup)) { unname(source_lookup[[key]]) } else { NA_character_ }
    }
    pairs <- unique(PlotDT_full[, c(XColumn, YColumn), with = FALSE])
    pairs[, `:=`(.dep_source = vapply(get(XColumn), resolve_source_col, character(1)), .ind_source = vapply(get(YColumn), resolve_source_col, character(1)) )]
    pairs[, N_pair := mapply(function(dep, ind){
      if(anyNA(c(dep, ind))){ return(NA_integer_) }
      sum(stats::complete.cases(sd[, c(dep, ind), with = FALSE]))
    }, .dep_source, .ind_source)]
    pairs[, N_x := vapply(.dep_source, function(col){
      if(is.na(col)){ return(NA_integer_) }
      sum(!is.na(sd[[col]]))
    }, integer(1))]
    pairs[, N_y := vapply(.ind_source, function(col){
      if(is.na(col)) return(NA_integer_)
      sum(!is.na(sd[[col]]))
    }, integer(1))]
    pairs[, N := fifelse(!is.na(N_pair), N_pair, fifelse(!is.na(N_x), N_x, N_y))]
    agg <- function(col, n_col) {
      sub <- pairs[!is.na(get(n_col))]
      if(nrow(sub) == 0){
        empty <- data.table(character(0), numeric(0), numeric(0), numeric(0))
        setnames(empty, c(col, "medN", "minN", "maxN"))
        return(empty)
      }
      sub[, .(medN = as.numeric(median(get(n_col))),
              minN = as.numeric(min(get(n_col))),
              maxN = as.numeric(max(get(n_col)))), by = col]
    }
    MarginalN <- list(x = agg(XColumn, "N_x"), y = agg(YColumn, "N_y"))
    if(all(is.na(pairs$N_x)) && all(is.na(pairs$N_y))){
      unmatched <- unique(c(pairs[is.na(.dep_source), get(XColumn)], pairs[is.na(.ind_source), get(YColumn)]))
      message(sprintf("SourceData supplied but no feature names matched its columns; marginal bars skipped. First unmatched names: %s",
                      paste(utils::head(unmatched, 8), collapse = ", ")))
      ShowMarginalBars <- FALSE }
  } else if(isTRUE(ShowMarginalBars) && is.null(SourceData)){
    message("ShowMarginalBars = TRUE but SourceData is NULL; marginal bars skipped.")
    ShowMarginalBars <- FALSE
  }
  apply_subset <- function(data, column, expression, argument) {
    if(length(expression) == 1L && is.na(expression)){ return(data) }
    if(!is.character(expression) || length(expression) != 1L){
      stop(argument, " must be NA or one comparison such as '<0.05'.") }
    match <- regexec( "^(<=|>=|==|!=|<|>)\\s*(-?(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)(?:[eE][+-]?[0-9]+)?)$", trimws(expression), perl = TRUE )
    parts <- regmatches(trimws(expression), match)[[1L]]
    if(length(parts) != 3L){ stop(argument, " must be a comparison such as '<0.05' or '>=1'.") }
    value <- as.numeric(parts[3L])
    x <- data[[column]]
    keep <- switch(parts[2L], "<" = x < value, "<=" = x <= value,
                   ">" = x > value, ">=" = x >= value,
                   "==" = x == value, "!=" = x != value)
    data[!is.na(keep) & keep]
  }
  PlotDT <- apply_subset(PlotDT, ColorVariable, ColorSubsetLevel, "ColorSubsetLevel")
  PlotDT <- apply_subset(PlotDT, SizeVariable, SizeSubsetLevel, "SizeSubsetLevel")
  if(nrow(PlotDT) == 0){ stop("No rows remain after subsetting. Check SizeSubsetLevel and ColorSubsetLevel.") }
  if(log2ColorVariable){
    if (any(PlotDT[[ColorVariable]] < 0, na.rm = TRUE)){
      stop("ColorVariable contains negative values and cannot be -log2 transformed.") }
    PlotDT[[ColorVariable]] <- -log2(pmax(PlotDT[[ColorVariable]], .Machine$double.xmin))
  }
  if(log2SizeVariable){
    if (any(PlotDT[[SizeVariable]] < 0, na.rm = TRUE)){
      stop("SizeVariable contains negative values and cannot be -log2 transformed.") }
    PlotDT[[SizeVariable]] <- -log2(pmax(PlotDT[[SizeVariable]], .Machine$double.xmin))
  }
  if(!is.null(SigColumn) && SigColumn %in% names(PlotDT)) {
    sig_values <- PlotDT[[SigColumn]]
    flag <- if(is.logical(sig_values)){ sig_values } else { sig_values < SigThreshold }
    set(PlotDT, j = ".SigFlag", value = flag)
  }
  PlotDT <- PlotDT[!is.na(get(SizeVariable)) & !is.na(get(ColorVariable))]
  if(nrow(PlotDT) == 0){ stop("No plottable rows remain after removing NA size/color values.") }
  x_lv <- x_order[x_order %in% unique(PlotDT[[XColumn]])]
  y_lv <- y_order[y_order %in% unique(PlotDT[[YColumn]])]
  PlotDT[[XColumn]] <- factor(PlotDT[[XColumn]], levels = x_lv)
  PlotDT[[YColumn]] <- factor(PlotDT[[YColumn]], levels = y_lv)
  ColorLabel <- if(log2ColorVariable){ paste0("-log2(", ColorVariable, ")") } else { ColorVariable }
  SizeLabel <- if(log2SizeVariable){ paste0("-log2(", SizeVariable,  ")") } else { SizeVariable }
  if(is.null(Xlab)){ Xlab <- gsub("_", " ", XColumn) }
  if(is.null(Ylab)){ Ylab <- gsub("_", " ", YColumn) }
  cvals <- PlotDT[[ColorVariable]]; cvals <- cvals[is.finite(cvals)]
  if (length(cvals) == 0L){ stop("ColorVariable has no finite values to plot.") }
  lim   <- stats::quantile(abs(cvals), ColorLimitQuantile, na.rm = TRUE)
  if(!is.finite(lim) || lim == 0){ lim <- max(abs(cvals), 1, na.rm = TRUE) }
  color_limits <- if(length(cvals) && min(cvals) < 0){ c(-lim, lim) } else { c(0, lim) }
  shown_empty <- empty_mat[y_lv, x_lv, drop = FALSE]
  if(length(shown_empty) && any(shown_empty)){
    idx <- which(shown_empty, arr.ind = TRUE)
    grid_bg <- data.table(X = factor(colnames(shown_empty)[idx[, "col"]], levels = x_lv),
                          Y = factor(rownames(shown_empty)[idx[, "row"]], levels = y_lv) )
  } else {
    grid_bg <- data.table(X = factor(character(0), levels = x_lv),
                          Y = factor(character(0), levels = y_lv))
  }
  P2 <- ggplot(PlotDT, aes(x = .data[[XColumn]], y = .data[[YColumn]])) +
    geom_point(data = grid_bg, aes(x = X, y = Y), shape = 21,
               size = 4, fill = NAFill, colour = NAFill, alpha = 0.65,
               inherit.aes = FALSE, show.legend = FALSE) +
    geom_point(aes(size = .data[[SizeVariable]], color = .data[[ColorVariable]])) +
    scale_color_distiller(palette = "RdBu", direction = 1,
                          limits = color_limits, oob = scales::squish) +
    scale_size_continuous(range = c(1, 10)) +
    scale_x_discrete(limits = x_lv, drop = FALSE) +
    scale_y_discrete(limits = y_lv, drop = FALSE) +
    theme_minimal() +
    labs(x = Xlab, y = Ylab, color = ColorLabel, size = SizeLabel) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  if(!is.null(SigColumn) && ".SigFlag" %in% names(PlotDT)){
    P2 <- P2 + geom_point(data = PlotDT[.SigFlag == TRUE],
                          aes(size = .data[[SizeVariable]]),
                          shape = 21, fill = NA, colour = "black", stroke = 1.1,
                          show.legend = FALSE)
    } 
  FinalPlot <- P2
  if(ShowMarginalBars && requireNamespace("patchwork", quietly = TRUE)){
    xN <- MarginalN$x[get(XColumn) %in% x_lv]; xN[[XColumn]] <- factor(xN[[XColumn]], levels = x_lv)
    yN <- MarginalN$y[get(YColumn) %in% y_lv]; yN[[YColumn]] <- factor(yN[[YColumn]], levels = y_lv)
    show_top   <- MarginalAxes %in% c("both", "x") && nrow(xN) > 0
    show_right <- MarginalAxes %in% c("both", "y") && nrow(yN) > 0
    if(show_top){
      top_bar <- ggplot(xN, aes(x = .data[[XColumn]], y = medN)) +
        geom_col(fill = "grey55", width = 0.45) +
        geom_errorbar(aes(ymin = minN, ymax = maxN), width = 0.3, colour = "grey30") +
        scale_x_discrete(limits = x_lv, drop = FALSE) +
        labs(y = "Feature N") + theme_minimal() +
        theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
              panel.grid.major.x = element_blank())
    }
    if(show_right){
      right_bar <- ggplot(yN, aes(y = .data[[YColumn]], x = medN)) +
        geom_col(fill = "grey55", width = 0.45) +
        geom_errorbar(aes(xmin = minN, xmax = maxN), width = 0.3, colour = "grey30") +
        scale_y_discrete(limits = y_lv, drop = FALSE) +
        labs(x = "Feature N") + theme_minimal() +
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
  } else if (ShowMarginalBars) {
    message("Package 'patchwork' is not installed; marginal bars were skipped.")
  }
  return(list(Plot = FinalPlot, BubblePlot = P2, Table = CombinedDF,
              ClusterOrder = list(x = x_lv, y = y_lv), MarginalN = MarginalN,
              ColorLimits = color_limits))
}
