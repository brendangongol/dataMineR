#' PCA plot of a numeric feature matrix
#'
#' Runs \code{\link[stats]{prcomp}} on the complete-case rows of
#' \code{features}, and plots PC1 vs PC2, optionally colored by
#' \code{GroupCol} and labeled by \code{LabelCol} (via \pkg{ggrepel}).
#'
#' @param dt A data.table containing \code{features} and, if supplied,
#'   \code{GroupCol}/\code{LabelCol}.
#' @param features Character vector of numeric column names to include in
#'   the PCA.
#' @param GroupCol Optional column name used to color points by group.
#' @param LabelCol Optional column name used to label points.
#' @param ScaleData Logical; scale features to unit variance before PCA
#'   (passed to \code{prcomp(scale. = )}).
#' @return A list with \code{pca} (the \code{prcomp} object),
#'   \code{VarianceExplained} (percent variance explained by PC1/PC2),
#'   \code{Scores} (the plotted data.frame), and \code{Plot} (a ggplot
#'   object).
#' @examples
#' \dontrun{
#' data(copd_data)
#' PCAPlotR(copd_data, features = c("FEV1", "FVC", "sao2", "AHI"),
#'          GroupCol = "group")
#' }
#' @export
PCAPlotR <- function(dt, features, GroupCol = NULL, LabelCol = NULL, ScaleData = TRUE){
  dt <- as.data.table(dt)
  needed <- c(features, GroupCol, LabelCol)
  missing_cols <- setdiff(needed, names(dt))
  if(length(missing_cols) > 0L){ stop("Columns not found in dt: ", paste(missing_cols, collapse = ", ")) }
  if(length(features) < 2L){ stop("PCA requires at least two feature columns.") }
  keep_idx <- stats::complete.cases(dt[, features, with = FALSE])
  mat <- as.matrix(dt[keep_idx, features, with = FALSE])
  if(!is.numeric(mat) || nrow(mat) < 2L){ stop("PCA requires at least two complete numeric rows.") }
  if(ScaleData && any(apply(mat, 2L, stats::sd) == 0)){ stop("PCA cannot scale a feature with zero variance.") }
  pca <- stats::prcomp(mat, scale. = ScaleData, center = TRUE)
  var_explained <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  scores <- as.data.frame(pca$x[, 1:2, drop = FALSE]); colnames(scores) <- c("PC1", "PC2")
  if(!is.null(GroupCol)){ scores$Group <- as.character(dt[[GroupCol]][keep_idx]) }
  if(!is.null(LabelCol)){ scores$Label <- as.character(dt[[LabelCol]][keep_idx]) }
  p <- if(!is.null(GroupCol)){
    ggplot2::ggplot(scores, ggplot2::aes(PC1, PC2, fill = Group)) + ggplot2::geom_point(shape = 21, size = 3, color = "black")
  } else {
    ggplot2::ggplot(scores, ggplot2::aes(PC1, PC2)) + ggplot2::geom_point(shape = 21, size = 3, fill = "grey40", color = "black")
  }
  p <- p + ggplot2::labs(x = paste0("PC1 (", var_explained[1], "%)"), y = paste0("PC2 (", var_explained[2], "%)"), fill = GroupCol) +
    ggpubr::theme_pubr()
  if(!is.null(LabelCol)){
    if(!requireNamespace("ggrepel", quietly = TRUE)) stop("Package 'ggrepel' is required when LabelCol is supplied.")
    p <- p + ggrepel::geom_label_repel(ggplot2::aes(label = Label), show.legend = FALSE)
  }
  return(list(pca = pca, VarianceExplained = var_explained, Scores = scores, Plot = p))
}

#' UMAP embedding of a numeric feature matrix
#'
#' Runs \pkg{umap} on the complete-case rows of \code{features} (optionally
#' scaled), and plots the 2D embedding, optionally colored by
#' \code{GroupCol} and labeled by \code{LabelCol} (via \pkg{ggrepel}).
#'
#' @param dt A data.table containing \code{features} and, if supplied,
#'   \code{GroupCol}/\code{LabelCol}.
#' @param features Character vector of numeric column names to embed.
#' @param GroupCol Optional column name used to color points by group.
#' @param LabelCol Optional column name used to label points.
#' @param ScaleData Logical; scale features (to mean 0, unit variance)
#'   before embedding.
#' @param n_neighbors UMAP neighborhood size (capped at \code{nrow - 1}).
#' @param seed Random seed for reproducibility.
#' @return A list with \code{umap} (the \pkg{umap} result object),
#'   \code{Embedding} (the plotted data.frame), and \code{Plot} (a ggplot
#'   object).
#' @examples
#' \dontrun{
#' data(copd_data)
#' UMAPPlotR(copd_data, features = c("FEV1", "FVC", "sao2", "AHI"),
#'           GroupCol = "group")
#' }
#' @export
UMAPPlotR <- function(dt, features, GroupCol = NULL, LabelCol = NULL, ScaleData = TRUE, n_neighbors = 15, seed = 42){
  if(!requireNamespace("umap", quietly = TRUE)){ stop("Package 'umap' is required for UMAPPlotR().") }
  dt <- as.data.table(dt)
  needed <- c(features, GroupCol, LabelCol)
  missing_cols <- setdiff(needed, names(dt))
  if(length(missing_cols) > 0L){ stop("Columns not found in dt: ", paste(missing_cols, collapse = ", ")) }
  if(length(features) < 2L){ stop("UMAP requires at least two feature columns.") }
  if(!is.numeric(n_neighbors) || length(n_neighbors) != 1L || n_neighbors < 2L || n_neighbors != as.integer(n_neighbors)){
    stop("n_neighbors must be an integer of at least 2.") }
  keep_idx <- stats::complete.cases(dt[, features, with = FALSE])
  mat <- as.matrix(dt[keep_idx, features, with = FALSE])
  if(!is.numeric(mat) || nrow(mat) < 3L){ stop("UMAP requires at least three complete numeric rows.") }
  if(ScaleData && any(apply(mat, 2L, stats::sd) == 0)){ stop("UMAP cannot scale a feature with zero variance.") }
  if(ScaleData){ mat <- scale(mat) }
  set.seed(seed)
  cfg <- umap::umap.defaults; cfg$n_neighbors <- min(as.integer(n_neighbors), nrow(mat) - 1L)
  ur <- umap::umap(mat, config = cfg)
  embed <- as.data.frame(ur$layout); colnames(embed) <- c("UMAP1", "UMAP2")
  if(!is.null(GroupCol)){ embed$Group <- as.character(dt[[GroupCol]][keep_idx]) }
  if(!is.null(LabelCol)){ embed$Label <- as.character(dt[[LabelCol]][keep_idx]) }
  p <- if(!is.null(GroupCol)){
    ggplot2::ggplot(embed, ggplot2::aes(UMAP1, UMAP2, color = Group)) + ggplot2::geom_point(size = 2.5)
  } else {
    ggplot2::ggplot(embed, ggplot2::aes(UMAP1, UMAP2)) + ggplot2::geom_point(size = 2.5, color = "grey30")
  }
  p <- p + ggpubr::theme_pubr()
  if(!is.null(LabelCol)){
    if(!requireNamespace("ggrepel", quietly = TRUE)) stop("Package 'ggrepel' is required when LabelCol is supplied.")
    p <- p + ggrepel::geom_label_repel(ggplot2::aes(label = Label), show.legend = FALSE)
  }
  return(list(umap = ur, Embedding = embed, Plot = p))
}

#' Scatterplot-matrix (pairs plot) of a set of features
#'
#' Wraps \code{\link[GGally]{ggpairs}} with correlations in the upper
#' triangle, scatter/histogram panels in the lower triangle, and density
#' curves on the diagonal, optionally colored by \code{GroupCol}.
#'
#' @param dt A data.frame/data.table containing \code{features} and, if
#'   supplied, \code{GroupCol}.
#' @param features Character vector of column names to include in the
#'   matrix.
#' @param GroupCol Optional column name used to color panels by group.
#' @return A \code{ggmatrix} object (from \pkg{GGally}).
#' @examples
#' \dontrun{
#' data(copd_data)
#' PairsPlotR(copd_data, features = c("FEV1", "FVC", "sao2", "AHI"),
#'            GroupCol = "group")
#' }
#' @export
PairsPlotR <- function(dt, features, GroupCol = NULL){
  if(!requireNamespace("GGally", quietly = TRUE)) stop("Package 'GGally' is required for PairsPlotR().")
  dt <- as.data.frame(dt)
  missing_cols <- setdiff(c(features, GroupCol), names(dt))
  if(length(missing_cols) > 0L){ stop("Columns not found in dt: ", paste(missing_cols, collapse = ", ")) }
  if(!is.null(GroupCol)){
    plot_dt <- dt[, c(features, GroupCol)]
    p <- GGally::ggpairs(plot_dt, columns = seq_along(features),
                         mapping = ggplot2::aes(color = .data[[GroupCol]]),
                         upper = list(continuous = GGally::wrap("cor", size = 3)),
                         lower = list(combo = GGally::wrap("facethist", bins = 30),
                                     continuous = GGally::wrap("points", alpha = 0.5, size = 0.8)),
                         diag = list(continuous = GGally::wrap("densityDiag", alpha = 0.5)))
  } else {
    plot_dt <- dt[, features]
    p <- GGally::ggpairs(plot_dt,
                         upper = list(continuous = GGally::wrap("cor", size = 3)),
                         lower = list(continuous = GGally::wrap("points", alpha = 0.5, size = 0.8)),
                         diag = list(continuous = GGally::wrap("densityDiag", alpha = 0.5)))
  }
  return(p)
}
