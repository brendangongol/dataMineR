#' Compare continuous variables between two pre-split data tables
#'
#' For each column named in \code{key}, runs a Wilcoxon rank-sum test (and,
#' where computable, a Cliff's delta effect size via
#' \code{\link{cliffs_delta_report}}) comparing its values in \code{DT1}
#' versus \code{DT2}. \code{\link{continuous_compare_Compile}} calls this
#' once per row of a \code{CatTable} to run it across every binary grouping
#' variable at once.
#'
#' @param DT1 First data.table (e.g. one level of a grouping variable).
#' @param DT2 Second data.table (e.g. the other level).
#' @param compair Character label identifying this comparison (e.g.
#'   \code{"Migrant1-Migrant0"}), carried through to the output.
#' @param key Character vector of continuous column names to test, present
#'   in both \code{DT1} and \code{DT2}.
#' @return A data.table with one row per element of \code{key}, including
#'   the Mann-Whitney U p-value and Cliff's delta effect size.
#' @export
continuous_compare <- function(DT1, DT2, compair, key) {
  missing_keys_DT1 <- setdiff(key, names(DT1))
  missing_keys_DT2 <- setdiff(key, names(DT2))
  if (length(missing_keys_DT1) > 0){
    stop(paste("Keys not found in DT1:", paste(missing_keys_DT1, collapse = ", "))) }
  if (length(missing_keys_DT2) > 0){
    stop(paste("Keys not found in DT2:", paste(missing_keys_DT2, collapse = ", "))) }
  non_numeric <- key[!vapply(key, function(column) {
    is.numeric(DT1[[column]]) && is.numeric(DT2[[column]])
  }, logical(1))]
  if (length(non_numeric) > 0L){
    stop("Continuous comparison columns must be numeric: ", paste(non_numeric, collapse = ", ")) }
  if (length(key) == 0L) {
    return(data.table( comparison = character(), parameter = character(), summary = character(),
                       W = numeric(), CliffsDelta = numeric(), VarghaDelaneyA = numeric(),
                       mann_whitney_u_pvalue = numeric(), p_bonferroni = numeric(), p_fdr = numeric(),
                       sig_unadjusted = logical(), sig_bonferroni = logical(), sig_fdr = logical(),
                       skipped = logical(), skip_reason = character() ))
  }
  test_dt <- NULL
  for (i in seq_along(key)){
    x <- DT1[[key[i]]]
    y <- DT2[[key[i]]]
    x <- x[!is.na(x)]
    y <- y[!is.na(y)]
    summarize_group <- function(values, label){
      if(length(values) == 0L){ return(sprintf("%s(n=0)", label)) }
      qs <- quantile(values, c(0.25, 0.75), names = FALSE)
      sprintf("%s(n=%d, median=%.3f, IQR=%.3f-%.3f)", label, length(values), median(values), qs[1L], qs[2L])
    }
    summary_str <- paste(summarize_group(x, "DT1"), summarize_group(y, "DT2"), sep = " | ")
    if (length(x) < 2 || length(y) < 2) {
      message(sprintf( "Key '%s': skipped - not enough non-NA observations (DT1 n=%d, DT2 n=%d).", key[i], length(x), length(y) ))
      test_dt <- rbind(test_dt, data.table( comparison = compair,
                                            parameter = key[i],
                                            summary = summary_str,
                                            W = NA_real_,
                                            CliffsDelta = NA_real_,
                                            VarghaDelaneyA = NA_real_,
                                            mann_whitney_u_pvalue = NA_real_,
                                            p_bonferroni = NA_real_,
                                            p_fdr = NA_real_,
                                            sig_unadjusted = NA,
                                            sig_bonferroni = NA,
                                            sig_fdr = NA,
                                            skipped = TRUE,
                                            skip_reason = sprintf("insufficient observations (DT1 n=%d, DT2 n=%d)", length(x), length(y))
                                          ), fill = TRUE)
      next
    }
    if(length(unique(x)) == 1 && length(unique(y)) == 1) {
      message(sprintf("Key '%s': skipped - zero variance in both groups.", key[i]))
      test_dt <- rbind(test_dt, data.table(
        comparison = compair,
        parameter = key[i],
        summary = summary_str,
        W = NA_real_,
        CliffsDelta = NA_real_,
        VarghaDelaneyA = NA_real_,
        mann_whitney_u_pvalue = NA_real_,
        p_bonferroni = NA_real_,
        p_fdr = NA_real_,
        sig_unadjusted = NA,
        sig_bonferroni = NA,
        sig_fdr = NA,
        skipped = TRUE,
        skip_reason = "zero variance in both groups"
      ), fill = TRUE)
      next
    }
    test_result <- tryCatch({
      test <- suppressWarnings( wilcox.test(x, y, conf.int = TRUE, paired = FALSE) )
      cliffsOutput <- cliffs_delta_report(x = x, y = y, wilcox_result = test)
      list( p_val = test$p.value, 
            W = test$statistic, 
            CliffsDelta = cliffsOutput$delta, 
            VarghaDelaneyA = cliffsOutput$v_a,
            error = NULL )
    }, error = function(e){
      list(p_val = NA_real_, W = NA_real_, CliffsDelta = NA_real_, 
           VarghaDelaneyA = NA_real_, error = conditionMessage(e))
    })
    skipped <- !is.null(test_result$error)
    skip_reason <- if(skipped){test_result$error } else {NA_character_ }
    if (skipped) message(sprintf("Key '%s': test error - %s", key[i], skip_reason))
    #### Assemble result row ####
    test_dt <- rbind(test_dt, data.table(comparison = compair,
                                         parameter = key[i],
                                         summary = summary_str,
                                         W = test_result$W,
                                         CliffsDelta = test_result$CliffsDelta,
                                         VarghaDelaneyA = test_result$VarghaDelaneyA,
                                         mann_whitney_u_pvalue = test_result$p_val,
                                         p_bonferroni = NA_real_,   # filled after loop
                                         p_fdr = NA_real_,   # filled after loop
                                         sig_unadjusted = NA,
                                         sig_bonferroni = NA,
                                         sig_fdr = NA,
                                         skipped = skipped,
                                         skip_reason = skip_reason ), fill = TRUE)
  }
  #### Multiple testing corrections (non-skipped rows only) ####
  all_p <- test_dt[["mann_whitney_u_pvalue"]]
  test_dt[, p_bonferroni := p.adjust(all_p, method = "bonferroni")]
  test_dt[, p_fdr := p.adjust(all_p, method = "BH")]
  test_dt[, sig_unadjusted := all_p < 0.05]
  test_dt[, sig_bonferroni := p_bonferroni < 0.05]
  test_dt[, sig_fdr := p_fdr < 0.05]
  return(test_dt)
}

#' Cliff's delta effect size from a completed Wilcoxon test
#'
#' Computes Vargha and Delaney's A (probability of superiority) and Cliff's
#' delta directly from the U statistic of an already-run
#' \code{\link[stats]{wilcox.test}}, avoiding an O(n^2) recomputation.
#'
#' @param x Numeric vector, group 1 (same vector passed to the original
#'   \code{wilcox.test} call).
#' @param y Numeric vector, group 2.
#' @param wilcox_result The \code{htest} object returned by
#'   \code{\link[stats]{wilcox.test}(x, y, ...)}.
#' @return A data.frame/list summarizing Vargha and Delaney's A and Cliff's
#'   delta (magnitude and interpretation) for the comparison.
#' @export
cliffs_delta_report <- function(x, y, wilcox_result){
  # Clean missing values
  x <- na.omit(x)
  y <- na.omit(y)
  # Calculate sample sizes
  n_x <- length(x)
  n_y <- length(y)
  # Extract values from the wilcox.test object
  w_stat <- wilcox_result$statistic
  p_val <- wilcox_result$p.value
  # Vargha and Delaney's A (Probability of Superiority). R's wilcox.test statistic W is exactly the Mann-Whitney U statistic
  v_a <- w_stat / (n_x * n_y)
  # Cliff's Delta. Directly calculated from Vargha & Delaney's A to avoid expensive O(N^2) outer loops
  delta <- 2 * v_a - 1
  # Rank-Biserial Correlation (r_rb). For two independent samples, rank-biserial correlation is identical to Cliff's Delta
  r_rb <- delta
  # Format the final APA style snippet
  report_str <- sprintf("W = %.1f, p = %.4f, Cliff's delta = %.3f, r = %.3f, A = %.3f", w_stat, p_val, delta, r_rb, v_a )
  # Return a list with all metrics
  return(list(
    delta = delta,
    r_biserial = r_rb,
    v_a = v_a,
    sample_sizes = c(n_x = n_x, n_y = n_y),
    report = report_str
  ))
}

#' Run \code{\link{continuous_compare}} across every binary grouping variable
#'
#' For each row of \code{CatTable} (as produced by
#' \code{\link{CategoricalVariableCreate}}), splits \code{dt} into its two
#' levels and runs \code{\link{continuous_compare}} across all of
#' \code{ContinuousCols}, writing one Excel file per grouping variable to
#' \code{Fpath} and returning the compiled results.
#'
#' @param dt A data.table containing the columns named in \code{CatTable} and
#'   \code{ContinuousCols}.
#' @param ContinuousCols Character vector of continuous column names to
#'   compare.
#' @param CatTable A data.table of binary grouping variables, as produced by
#'   \code{\link{CategoricalVariableCreate}} (columns
#'   \code{CategoricalVariable}, \code{Table1Subset}, \code{Table2Subset},
#'   \code{NameOrientation}).
#' @param Fpath Directory to write one \code{ContinuousComparisons_*.xlsx}
#'   file per grouping variable.
#' @return A data.table combining the \code{\link{continuous_compare}}
#'   results across every row of \code{CatTable}.
#' @examples
#' \dontrun{
#' load(system.file("extdata", "dfMigrant.rda", package = "dataMineR"))
#' CatTable <- CategoricalVariableCreate(dfMigrant, c("Migrant", "Homeless"))
#' continuous_compare_Compile(dfMigrant, ContinuousCols = "Age",
#'                             CatTable = CatTable, Fpath = tempdir())
#' }
#' @export
continuous_compare_Compile <- function(dt, ContinuousCols, CatTable, Fpath){
  CompiledAll <- data.table()
  if (!dir.exists(Fpath)) stop("Fpath does not exist: ", Fpath)
  for(i in seq_len(nrow(CatTable))){
    DTOne <- dt[as.character(dt[[CatTable$CategoricalVariable[i]]]) == CatTable$Table1Subset[i],]
    DTTwo <- dt[as.character(dt[[CatTable$CategoricalVariable[i]]]) == CatTable$Table2Subset[i],]
    CompsT <- continuous_compare(DT1 = DTOne, DT2 = DTTwo, compair = CatTable$NameOrientation[i], key = ContinuousCols)
    Fname <- paste("ContinuousComparisons", paste(gsub("-", "", CatTable$NameOrientation[i]), ".xlsx", sep = ""), sep = "_")
    openxlsx::write.xlsx(list(ContinuousComparisons = CompsT), file = file.path(Fpath, Fname), colWidths = "auto" )
    CompiledAll <- rbind(CompiledAll, CompsT)
    }
  return(CompiledAll)
}
