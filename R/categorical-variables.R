#' Auto-detect binary categorical variables and build a comparison table
#'
#' Scans each column named in \code{CategoricalDependents} and keeps only
#' those with exactly two distinct non-missing, non-blank values, recording
#' the two levels as \code{Table1Subset}/\code{Table2Subset} and a combined
#' \code{NameOrientation} label (e.g. \code{"Migrant1-Migrant0"}). The result
#' drives \code{\link{continuous_compare_Compile}},
#' \code{\link{Categorical_Compare_compile}}, and the logistic-regression
#' reference level in \code{\link{RegressionCompile}}.
#'
#' @param dt A data.table containing the candidate columns.
#' @param CategoricalDependents Character vector of column names to screen.
#' @return A data.table with columns \code{CategoricalVariable},
#'   \code{Table1Subset}, \code{Table2Subset}, and \code{NameOrientation}, one
#'   row per binary variable found.
#' @examples
#' \dontrun{
#' data(dfMigrant)
#' CategoricalVariableCreate(dfMigrant, c("Migrant", "Homeless", "Sex"))
#' }
#' @export
CategoricalVariableCreate <- function(dt, CategoricalDependents) {
  dt <- as.data.table(dt)
  missing_cols <- setdiff(CategoricalDependents, names(dt))
  if (length(missing_cols) > 0L) {
    stop("CategoricalDependents not found in dt: ",
         paste(missing_cols, collapse = ", "))
  }
  rows <- lapply(CategoricalDependents, function(column){
    levels <- sort(unique(as.character(dt[[column]])))
    levels <- levels[!is.na(levels) & nzchar(trimws(levels))]
    if(length(levels) != 2L){ return(NULL) }
    data.table(CategoricalVariable = column,
               Table1Subset = levels[1L],
               Table2Subset = levels[2L],
               NameOrientation = paste0(column, levels, collapse = "-") )  })
  rbindlist(rows, fill = TRUE)
}
