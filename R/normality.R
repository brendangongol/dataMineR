#' Shapiro-Wilk normality test across many columns
#'
#' Runs \code{\link[stats]{shapiro.test}} on each column in \code{cols} and
#' reports whether the normality null hypothesis is rejected at p < 0.05.
#'
#' @param DT A data.table (or data.frame) containing the columns in
#'   \code{cols}.
#' @param cols Character vector of numeric column names to test.
#' @return A data.table with columns \code{method}, \code{pvalue},
#'   \code{variable}, and \code{normality} (\code{"normally distributed"} or
#'   \code{"Not normally distributed"}), one row per column in \code{cols}.
#' @examples
#' \dontrun{
#' load(system.file("extdata", "dfMigrant.rda", package = "dataMineR"))
#' normalityfindR(dfMigrant, cols = c("Age", "COMORBIDITYINDEX"))
#' }
#' @export
normalityfindR <- function(DT, cols){
  missing_cols <- setdiff(cols, names(DT))
  if(length(missing_cols) > 0L){ stop("Columns not found in DT: ", paste(missing_cols, collapse = ", ")) }
  non_numeric <- cols[!vapply(cols, function(column) is.numeric(DT[[column]]), logical(1))]
  if(length(non_numeric) > 0L){ stop("cols must contain numeric columns: ", paste(non_numeric, collapse = ", ")) }
  finDT <- data.table()
  for(i in seq_along(cols)){
    select <- DT[[cols[i]]]
    select <- select[!is.na(select)]
    variable_name <- colnames(DT[, cols[i], with = FALSE])
    n_valid <- length(select)
    if(length(unique(select)) <= 1){
      DT2 <- data.table(method = NA_character_,
                        pvalue = NA_real_,
                        variable = variable_name,
                        normality = "Not tested: fewer than 2 unique non-missing values")
      finDT <- rbind(finDT, DT2)
      next
    }
    if(n_valid < 3 || n_valid > 5000){
      DT2 <- data.table(method = "Shapiro-Wilk normality test",
                        pvalue = NA_real_,
                        variable = variable_name,
                        normality = sprintf("Not tested: Shapiro-Wilk requires 3 to 5000 non-missing values (n=%d)", n_valid))
      finDT <- rbind(finDT, DT2)
      next
    }
    test <- shapiro.test(select)
    if(test$p.value < 0.05){
      distribution <- "Not normally distributed"
    } else { distribution <- "normally distributed" }
    DT2 <- data.table(method = test$method,
                      pvalue = test$p.value,
                      variable = variable_name,
                      normality = distribution)
    finDT <- rbind(finDT, DT2)
  }
  return(finDT) }
