#' Scatterplot of two continuous variables
#'
#' Simple \pkg{ggplot2} scatterplot helper used to visualize individual
#' variable pairs, e.g. the significant results from
#' \code{\link{flattenCorrMatrix}}.
#'
#' @param dataset A data.frame or data.table containing \code{xaxis} and
#'   \code{yaxis}.
#' @param xaxis Character column name to plot on the x axis.
#' @param yaxis Character column name to plot on the y axis.
#' @param title Plot title.
#' @param xlabel X-axis label.
#' @param ylabel Y-axis label.
#' @return A ggplot object.
#' @export
scatterplotgraph <- function(dataset, xaxis, yaxis, title, xlabel, ylabel){
  missing_cols <- setdiff(c(xaxis, yaxis), names(dataset))
  if(length(missing_cols) > 0L){ stop("Columns not found in dataset: ", paste(missing_cols, collapse = ", ")) }
  p <- ggplot(dataset, aes(x = .data[[xaxis]], y = .data[[yaxis]]))+
  geom_point(shape = 21, size = 3, colour = "black", fill = "#08519C")+
  geom_smooth(method = "lm")+
  ggtitle(title)+ xlab(xlabel) + ylab(ylabel) +
  theme_pubr()
return(p)
}
