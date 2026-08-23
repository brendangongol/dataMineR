#' Flatten a correlation matrix into a long-format table
#'
#' Converts the upper triangle of a correlation matrix (and its matching
#' p-value, and optionally sample-size matrices, such as the \code{r}, \code{P},
#' and \code{n} components returned by \code{Hmisc::rcorr()}) into one row
#' per variable pair.
#'
#' @param cormat A square correlation matrix with row/column names.
#' @param pmat A matrix of p-values, same dimensions as \code{cormat}.
#' @param nmat Optional matrix of pairwise sample sizes, same dimensions as
#'   \code{cormat}. If supplied, an \code{n} column is included in the output.
#' @return A data.frame with one row per variable pair: \code{row},
#'   \code{column}, \code{cor}, optionally \code{n}, and \code{p}.
#' @examples
#' \dontrun{
#' res <- Hmisc::rcorr(as.matrix(dfMigrant[, c("Age", "COMORBIDITYINDEX")]))
#' flattenCorrMatrix(res$r, res$P, res$n)
#' }
#' @export
flattenCorrMatrix <- function(cormat, pmat, nmat = NULL){
  if (!is.matrix(cormat) || nrow(cormat) != ncol(cormat))
    stop("cormat must be a square matrix.")
  if (!is.matrix(pmat) || !identical(dim(pmat), dim(cormat)))
    stop("pmat must be a matrix with the same dimensions as cormat.")
  if (!is.null(nmat) && (!is.matrix(nmat) || !identical(dim(nmat), dim(cormat))))
    stop("nmat must be NULL or a matrix with the same dimensions as cormat.")
  if (is.null(rownames(cormat)) || is.null(colnames(cormat)))
    stop("cormat must have row and column names.")
  ut <- upper.tri(cormat)
  if(is.null(nmat)){
    data.frame(row = rownames(cormat)[row(cormat)[ut]],
               column = colnames(cormat)[col(cormat)[ut]],
               cor  =(cormat)[ut],
               p = pmat[ut] )
  } else {
    data.frame(row = rownames(cormat)[row(cormat)[ut]],
               column = colnames(cormat)[col(cormat)[ut]],
               cor =(cormat)[ut],
               n =  (nmat)[ut],
               p = pmat[ut] ) }
}
