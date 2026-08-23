#' Compare categorical variables against a grouping column
#'
#' For each column named in \code{keys}, cross-tabulates it against
#' \code{compair_column} and runs a chi-squared or Fisher's exact test
#' (falling back to Monte Carlo simulation for larger contingency tables),
#' reporting odds ratios, risk ratios, Cramer's V, and effect-size
#' confidence intervals (via \pkg{effectsize}) alongside multiple-testing
#' corrected p-values (Bonferroni and FDR).
#'
#' @param DT A data.table containing \code{compair_column} and every column
#'   in \code{keys}.
#' @param compair Character label identifying this comparison, carried
#'   through to the output.
#' @param compair_column Column name of the (typically binary) grouping
#'   variable to compare each key against.
#' @param keys Character vector of categorical column names to test.
#' @param Stattest One of \code{"chi"} or \code{"fisher"}.
#' @param ci_level Confidence level used for effect-size intervals.
#' @param digits Number of digits used when formatting effect estimates.
#' @return A data.table with one row per element of \code{keys}, including
#'   odds ratio, risk ratio, Cramer's V, and (un)corrected p-values.
#' @note The multi-level (more than two columns/rows) comparison path can be
#'   unreliable for some contingency-table shapes; prefer binary
#'   (two-level) \code{keys} where possible.
#' @export
categorical_compare <- function(DT, compair, compair_column, keys, Stattest = "chi", ci_level = 0.95, digits = 4){
  #### Input validation ####
  if (!requireNamespace("effectsize", quietly = TRUE))
    stop("Package 'effectsize' is required. Please install it using install.packages('effectsize')")
  if (!compair_column %in% names(DT))
    stop(paste("compair_column not found in DT:", compair_column))
  missing_keys <- setdiff(keys, names(DT))
  if (length(missing_keys) > 0)
    stop(paste("Keys not found in DT:", paste(missing_keys, collapse = ", ")))
  if (!Stattest %in% c("chi", "fisher"))
    stop('Stattest must be "chi" or "fisher"')
  if (!is.numeric(ci_level) || ci_level <= 0 || ci_level >= 1)
    stop("ci_level must be a numeric value between 0 and 1 (exclusive)")
  if (!is.numeric(digits) || digits < 0)
    stop("digits must be a non-negative integer")
  if (length(keys) == 0L) return(data.table())
  fmt_effect <- function(val, low, high, d = digits){
    if (is.na(val)) return(NA_character_)
    sprintf(paste0("%.", d, "f [%.", d, "f, %.", d, "f]"), val, low, high)
  }
  test_dt <- NULL
  p_col <- if (Stattest == "chi") "chi_squared_p_value" else "fisher_exact_p_value"
  for(i in seq_along(keys)){
    mytab <- table(DT[[compair_column]], DT[[keys[i]]])
    is_2x2 <- nrow(mytab) == 2 && ncol(mytab) == 2
    col_marginals <- colSums(mytab)
    row_marginals <- rowSums(mytab)
    if(sum(col_marginals > 0) < 2 || sum(row_marginals > 0) < 2){
      message(sprintf("Key '%s': skipped - table has fewer than 2 non-zero marginals (no variation).", keys[i] ))
      base_row <- data.table( comparison = compair,
                              comparison_type = "summary",
                              parameter = keys[i],
                              n_rows = nrow(mytab),
                              n_cols = ncol(mytab),
                              n_rows_clean = NA_integer_,
                              n_cols_clean = NA_integer_,
                              cell_counts = NA_character_,
                              skipped = TRUE,
                              skip_reason = "fewer than 2 non-zero marginals",
                              CramersV = NA_real_, CramersV_CI_low   = NA_real_, CramersV_CI_high  = NA_real_, CramersV_str    = NA_character_,
                              TschuprowsT = NA_real_, TschuprowsT_CI_low= NA_real_, TschuprowsT_CI_high= NA_real_, TschuprowsT_str= NA_character_,
                              CohensW = NA_real_, CohensW_CI_low    = NA_real_, CohensW_CI_high   = NA_real_, CohensW_str     = NA_character_,
                              odds_ratio = NA_real_, odds_ratio_CI_low = NA_real_, odds_ratio_CI_high= NA_real_, odds_ratio_str  = NA_character_,
                              OR_pairwise = NA_character_, OR_onevsrest = NA_character_,
                              best_OR_pairwise = NA_character_, best_OR_onevsrest = NA_character_,
                              risk_ratio = NA_real_, risk_ratio_CI_low = NA_real_, risk_ratio_CI_high= NA_real_, risk_ratio_str  = NA_character_,
                              RR_pairwise = NA_character_, RR_onevsrest = NA_character_,
                              haldane_correction = FALSE)
      base_row[[p_col]] <- NA_real_
      test_dt <- rbind(test_dt, base_row, fill = TRUE)
      next
    }
    test_result <- tryCatch({
      if(Stattest == "chi"){
        test <- chisq.test(mytab)
        p_val <- test$p.value
      } else {
        if(is_2x2){
          test <- fisher.test(mytab, conf.level = ci_level)
        } else {
          message(sprintf("Key '%s': table is %dx%d - using Fisher's test with Monte Carlo simulation (B=2000).", keys[i], nrow(mytab), ncol(mytab) ))
          test <- fisher.test(mytab, simulate.p.value = TRUE, B = 2000)
        }
        p_val <- test$p.value
      }
      mytab_clean <- mytab[rowSums(mytab) > 0, colSums(mytab) > 0, drop = FALSE]
      is_2x2_clean <- nrow(mytab_clean) == 2 && ncol(mytab_clean) == 2
      cv_obj <- suppressWarnings(effectsize::cramers_v(mytab_clean, ci = ci_level))
      tt_obj <- suppressWarnings(effectsize::tschuprows_t(mytab_clean, ci = ci_level))
      cramers_v <- round(cv_obj$Cramers_v, digits)
      cramers_v_low <- round(cv_obj$CI_low, digits)
      cramers_v_high <- round(cv_obj$CI_high, digits)
      tschuprows_t <- round(tt_obj$Tschuprows_t, digits)
      tschuprows_t_low <- round(tt_obj$CI_low, digits)
      tschuprows_t_high <- round(tt_obj$CI_high, digits)
      k <- min(nrow(mytab_clean), ncol(mytab_clean))
      w_scale <- sqrt(k - 1)
      cohens_w <- round(cramers_v * w_scale, digits)
      cohens_w_low <- round(cramers_v_low * w_scale, digits)
      cohens_w_high <- round(cramers_v_high * w_scale, digits)
      or_rr_2x2 <- function(tab2x2) {
        tab_use <- if (any(tab2x2 == 0L)) tab2x2 + 0.5 else tab2x2
        or_o <- suppressWarnings(effectsize::oddsratio(tab_use, alternative = "two.sided", ci = ci_level))
        rr_o <- suppressWarnings(effectsize::riskratio(tab_use, ci = ci_level))
        list(or = round(or_o$Odds_ratio, digits),
                        or_low = round(or_o$CI_low,     digits),
                        or_high = round(or_o$CI_high,    digits),
                        rr = round(rr_o$Risk_ratio, digits),
                        rr_low = round(rr_o$CI_low,     digits),
                        rr_high = round(rr_o$CI_high,    digits),
                        haldane = any(tab2x2 == 0L)  )
      }
      fmt_multi <- function(results_list, stat){
        parts <- vapply(results_list, function(r){
          if (is.null(r)){ return(NA_character_) }
          val <- r[[stat]]
          low <- r[[paste0(stat, "_low")]]
          high <- r[[paste0(stat, "_high")]]
          p <- if(!is.null(r$p_val) && !is.na(r$p_val)){ formatC(r$p_val, format = "g", digits = 4) } else { "NA" }
          sprintf("%s: %s=%s [%s, %s], p=%s", r$label, toupper(stat),
                  formatC(val, format = "f", digits = digits),
                  formatC(low, format = "f", digits = digits),
                  formatC(high, format = "f", digits = digits), p)
        }, character(1))
        paste(na.omit(parts), collapse = "; ")
      }
      col_levels <- colnames(mytab_clean)
      pairwise_list <- list()
      onevsrest_list <- list()
      pval_2x2 <- function(tab2x2) {
        tryCatch({
          if (Stattest == "chi") chisq.test(tab2x2)$p.value
          else fisher.test(tab2x2, conf.level = ci_level)$p.value
        }, error = function(e) NA_real_)
      }
      if(is_2x2_clean){
        res_2x2 <- or_rr_2x2(mytab_clean)
        haldane_applied <- res_2x2$haldane
        odds_ratio <- res_2x2$or
        odds_ratio_low <- res_2x2$or_low
        odds_ratio_high <- res_2x2$or_high
        risk_ratio <- res_2x2$rr
        risk_ratio_low <- res_2x2$rr_low
        risk_ratio_high <- res_2x2$rr_high
        pairwise_label <- paste(col_levels[2], "vs", col_levels[1])
        p_2x2_str <- formatC(p_val, format = "g", digits = 4)
        or_pairwise_str <- sprintf("%s: OR=%s [%s, %s], p=%s", pairwise_label,
                                       formatC(odds_ratio, format = "f", digits = digits),
                                       formatC(odds_ratio_low, format = "f", digits = digits),
                                       formatC(odds_ratio_high, format = "f", digits = digits),
                                       p_2x2_str)
        rr_pairwise_str <- sprintf("%s: RR=%s [%s, %s], p=%s", pairwise_label,
                                       formatC(risk_ratio, format = "f", digits = digits),
                                       formatC(risk_ratio_low, format = "f", digits = digits),
                                       formatC(risk_ratio_high, format = "f", digits = digits),
                                       p_2x2_str)
        or_onevsrest_str <- or_pairwise_str
        rr_onevsrest_str <- rr_pairwise_str
        best_or_pairwise_str <- or_pairwise_str
        best_or_onevsrest_str <- or_pairwise_str
      } else {
        odds_ratio <- odds_ratio_low <- odds_ratio_high <- NA_real_
        risk_ratio <- risk_ratio_low <- risk_ratio_high <- NA_real_
        haldane_applied <- FALSE
        ref_col <- col_levels[1]
        pairwise_list <- lapply(col_levels[-1], function(lv) {
          sub_tab <- mytab_clean[, c(ref_col, lv), drop = FALSE]
          sub_tab <- sub_tab[rowSums(sub_tab) > 0, , drop = FALSE]
          if(nrow(sub_tab) != 2 || ncol(sub_tab) != 2){ return(NULL) }
          res <- or_rr_2x2(sub_tab)
          res$label <- paste(lv, "vs", ref_col)
          res$p_val <- pval_2x2(sub_tab)
          res$sub_tab <- sub_tab
          if(res$haldane){haldane_applied <<- TRUE}
          res
        })
        onevsrest_list <- lapply(col_levels, function(lv) {
          rest_cols <- col_levels[col_levels != lv]
          rest_counts <- rowSums(mytab_clean[, rest_cols, drop = FALSE])
          sub_tab <- cbind(mytab_clean[, lv, drop = FALSE], rest_counts)
          colnames(sub_tab) <- c(lv, "rest")
          sub_tab <- sub_tab[rowSums(sub_tab) > 0, , drop = FALSE]
          if(nrow(sub_tab) != 2 || ncol(sub_tab) != 2){ return(NULL) }
          res <- or_rr_2x2(sub_tab)
          res$label <- paste(lv, "vs rest")
          res$p_val <- pval_2x2(sub_tab)
          res$sub_tab <- sub_tab
          if(res$haldane){ haldane_applied <<- TRUE }
          res
        })
        or_pairwise_str <- fmt_multi(pairwise_list, "or")
        rr_pairwise_str <- fmt_multi(pairwise_list, "rr")
        or_onevsrest_str <- fmt_multi(onevsrest_list, "or")
        rr_onevsrest_str <- fmt_multi(onevsrest_list, "rr")
        best_or_fmt <- function(res_list){
          valid <- Filter(Negate(is.null), res_list)
          if(length(valid) == 0){ return(NA_character_) }
          p_vals <- vapply(valid, function(r) r$p_val, numeric(1))
          if (all(is.na(p_vals))){ return(NA_character_) }
          best <- valid[[which.min(p_vals)]]
          sprintf("(best p=%.4g) %s: OR=%s [%s, %s]", best$p_val, best$label,
                  formatC(best$or,      format = "f", digits = digits),
                  formatC(best$or_low,  format = "f", digits = digits),
                  formatC(best$or_high, format = "f", digits = digits))
        }
        best_or_pairwise_str <- best_or_fmt(pairwise_list)
        best_or_onevsrest_str <- best_or_fmt(onevsrest_list)
      }
      list(p_val = p_val,
           cramers_v = cramers_v,
           cramers_v_low = cramers_v_low,
           cramers_v_high = cramers_v_high,
           tschuprows_t = tschuprows_t,
           tschuprows_t_low = tschuprows_t_low, 
           tschuprows_t_high = tschuprows_t_high,
           cohens_w = cohens_w,
           cohens_w_low = cohens_w_low,
           cohens_w_high = cohens_w_high,
           odds_ratio = odds_ratio,
           odds_ratio_low = odds_ratio_low,
           odds_ratio_high = odds_ratio_high,
           risk_ratio = risk_ratio,
           risk_ratio_low = risk_ratio_low,
           risk_ratio_high = risk_ratio_high,
           or_pairwise_str = or_pairwise_str,
           rr_pairwise_str = rr_pairwise_str,
           or_onevsrest_str = or_onevsrest_str,
           rr_onevsrest_str = rr_onevsrest_str,
           best_or_pairwise_str = best_or_pairwise_str,
           best_or_onevsrest_str = best_or_onevsrest_str,
           haldane_applied = haldane_applied,
           n_rows_clean = nrow(mytab_clean),
           n_cols_clean = ncol(mytab_clean),
           pairwise_list = pairwise_list,
           onevsrest_list = onevsrest_list,
           error = NULL )
    }, error = function(e) {
      list(p_val = NA_real_,
           cramers_v = NA_real_,
           cramers_v_low = NA_real_, 
           cramers_v_high = NA_real_,
           tschuprows_t = NA_real_, 
           tschuprows_t_low = NA_real_,
           tschuprows_t_high = NA_real_,
           cohens_w = NA_real_, 
           cohens_w_low = NA_real_,
           cohens_w_high = NA_real_,
           odds_ratio = NA_real_,
           odds_ratio_low = NA_real_,
           odds_ratio_high = NA_real_,
           risk_ratio = NA_real_,
           risk_ratio_low = NA_real_,
           risk_ratio_high = NA_real_,
           or_pairwise_str = NA_character_,
           rr_pairwise_str = NA_character_,
           or_onevsrest_str = NA_character_,
           rr_onevsrest_str = NA_character_,
           best_or_pairwise_str = NA_character_,
           best_or_onevsrest_str = NA_character_,
           haldane_applied = FALSE,
           n_rows_clean = NA_integer_,
           n_cols_clean = NA_integer_,
           pairwise_list = list(),
           onevsrest_list = list(),
           error = conditionMessage(e))  })
    skipped <- !is.null(test_result$error)
    skip_reason <- if(skipped){ test_result$error } else { NA_character_ }
    if(skipped){ message(sprintf("Key '%s': test error - %s", keys[i], skip_reason)) }
    is_multilevel <- isTRUE( !is.na(test_result$n_cols_clean) && test_result$n_cols_clean > 2 ||
                             !is.na(test_result$n_rows_clean) && test_result$n_rows_clean > 2 )
    cell_names <- apply(expand.grid(rownames(mytab), colnames(mytab)), 1, paste, collapse = "|")
    cell_values <- as.vector(mytab)
    cell_counts <- paste(paste(cell_names, cell_values, sep = "="), collapse = ", ")
    base_row <- data.table(comparison = compair,
                           comparison_type = "summary",
                           parameter = keys[i],
                           n_rows = nrow(mytab),
                           n_cols = ncol(mytab),
                           n_rows_clean = test_result$n_rows_clean,
                           n_cols_clean = test_result$n_cols_clean,
                           cell_counts = cell_counts,
                           skipped = skipped,
                           skip_reason = skip_reason,
                           CramersV = test_result$cramers_v,
                           CramersV_CI_low = test_result$cramers_v_low,
                           CramersV_CI_high = test_result$cramers_v_high,
                           CramersV_str = fmt_effect(test_result$cramers_v, test_result$cramers_v_low, test_result$cramers_v_high),
                           TschuprowsT = test_result$tschuprows_t,
                           TschuprowsT_CI_low = test_result$tschuprows_t_low,
                           TschuprowsT_CI_high = test_result$tschuprows_t_high,
                           TschuprowsT_str = fmt_effect(test_result$tschuprows_t, test_result$tschuprows_t_low, test_result$tschuprows_t_high),
                           CohensW = test_result$cohens_w,
                           CohensW_CI_low = test_result$cohens_w_low,
                           CohensW_CI_high = test_result$cohens_w_high,
                           CohensW_str = fmt_effect(test_result$cohens_w, test_result$cohens_w_low, test_result$cohens_w_high),
                           odds_ratio = test_result$odds_ratio,
                           odds_ratio_CI_low = test_result$odds_ratio_low,
                           odds_ratio_CI_high = test_result$odds_ratio_high,
                           odds_ratio_str = fmt_effect(test_result$odds_ratio, test_result$odds_ratio_low, test_result$odds_ratio_high),
                           OR_pairwise = test_result$or_pairwise_str,
                           OR_onevsrest = test_result$or_onevsrest_str,
                           best_OR_pairwise = test_result$best_or_pairwise_str,
                           best_OR_onevsrest = test_result$best_or_onevsrest_str,
                           risk_ratio = test_result$risk_ratio,
                           risk_ratio_CI_low = test_result$risk_ratio_low,
                           risk_ratio_CI_high = test_result$risk_ratio_high,
                           risk_ratio_str = fmt_effect(test_result$risk_ratio, test_result$risk_ratio_low, test_result$risk_ratio_high),
                           RR_pairwise = test_result$rr_pairwise_str,
                           RR_onevsrest = test_result$rr_onevsrest_str,
                           haldane_correction = test_result$haldane_applied )
    base_row[[p_col]] <- test_result$p_val
    if(is_multilevel && !skipped){
      pw_valid <- Filter(Negate(is.null), test_result$pairwise_list)
      if (length(pw_valid) > 0) { pw_pvals <- vapply(pw_valid, function(r) r$p_val, numeric(1))
        if(!all(is.na(pw_pvals))){
          best_pw  <- pw_valid[[which.min(pw_pvals)]]
          base_row[, odds_ratio        := best_pw$or]
          base_row[, odds_ratio_CI_low := best_pw$or_low]
          base_row[, odds_ratio_CI_high:= best_pw$or_high]
          base_row[, odds_ratio_str    := fmt_effect(best_pw$or, best_pw$or_low, best_pw$or_high)]
          base_row[, risk_ratio        := best_pw$rr]
          base_row[, risk_ratio_CI_low := best_pw$rr_low]
          base_row[, risk_ratio_CI_high:= best_pw$rr_high]
          base_row[, risk_ratio_str    := fmt_effect(best_pw$rr, best_pw$rr_low, best_pw$rr_high)]
        }
      }
    }
    test_dt <- rbind(test_dt, base_row, fill = TRUE)
    
    if(is_multilevel && !skipped){
      make_expanded_row <- function(res, comp_type) {
        if(is.null(res)){ return(NULL) }
        sub_tab <- res$sub_tab
        if(!is.null(sub_tab)){
          sn <- apply(expand.grid(rownames(sub_tab), colnames(sub_tab)), 1, paste, collapse = "|")
          sv <- as.vector(sub_tab)
          sub_cell_counts <- paste(paste(sn, sv, sep = "="), collapse = ", ")
          sub_n_rows <- nrow(sub_tab)
          sub_n_cols <- ncol(sub_tab)
        } else {
          sub_cell_counts <- NA_character_
          sub_n_rows <- NA_integer_
          sub_n_cols <- NA_integer_
        }
        or_str_single <- fmt_effect(res$or, res$or_low, res$or_high)
        rr_str_single <- fmt_effect(res$rr, res$rr_low, res$rr_high)
        p_fmt <- if(!is.null(res$p_val) && !is.na(res$p_val)){ formatC(res$p_val, format = "g", digits = 4) } else {"NA"}
        or_str_full <- sprintf("%s: OR=%s [%s, %s], p=%s", res$label,
                               formatC(res$or, format = "f", digits = digits),
                               formatC(res$or_low, format = "f", digits = digits),
                               formatC(res$or_high, format = "f", digits = digits),
                               p_fmt)
        rr_str_full <- sprintf("%s: RR=%s [%s, %s], p=%s", res$label,
                               formatC(res$rr, format = "f", digits = digits),
                               formatC(res$rr_low, format = "f", digits = digits),
                               formatC(res$rr_high, format = "f", digits = digits),
                               p_fmt)
        row <- data.table(comparison = compair,
                          comparison_type = comp_type,
                          parameter = paste(keys[i], res$label, sep = "_"),
                          n_rows = sub_n_rows,
                          n_cols = sub_n_cols,
                          n_rows_clean = sub_n_rows,
                          n_cols_clean = sub_n_cols,
                          cell_counts = sub_cell_counts,
                          skipped = FALSE,
                          skip_reason = NA_character_,
                          CramersV = base_row$CramersV,
                          CramersV_CI_low = base_row$CramersV_CI_low,
                          CramersV_CI_high = base_row$CramersV_CI_high,
                          CramersV_str = base_row$CramersV_str,
                          TschuprowsT = base_row$TschuprowsT,
                          TschuprowsT_CI_low = base_row$TschuprowsT_CI_low,
                          TschuprowsT_CI_high = base_row$TschuprowsT_CI_high,
                          TschuprowsT_str = base_row$TschuprowsT_str,
                          CohensW = base_row$CohensW,
                          CohensW_CI_low = base_row$CohensW_CI_low,
                          CohensW_CI_high = base_row$CohensW_CI_high,
                          CohensW_str = base_row$CohensW_str,
                          odds_ratio = res$or,
                          odds_ratio_CI_low = res$or_low,
                          odds_ratio_CI_high = res$or_high,
                          odds_ratio_str = or_str_single,
                          OR_pairwise = if(comp_type == "pairwise"){ or_str_full } else { NA_character_ },
                          OR_onevsrest = if(comp_type == "onevsrest"){ or_str_full } else { NA_character_ },
                          best_OR_pairwise = NA_character_,
                          best_OR_onevsrest = NA_character_,
                          risk_ratio = res$rr,
                          risk_ratio_CI_low = res$rr_low,
                          risk_ratio_CI_high = res$rr_high,
                          risk_ratio_str = rr_str_single,
                          RR_pairwise = if(comp_type == "pairwise"){ rr_str_full } else { NA_character_ },
                          RR_onevsrest = if(comp_type == "onevsrest"){ rr_str_full } else { NA_character_ },
                          haldane_correction = res$haldane )
        row[[p_col]] <- res$p_val
        row
      }
      expanded_pairwise  <- rbindlist(lapply(test_result$pairwise_list,  make_expanded_row, comp_type = "pairwise"), fill = TRUE)
      expanded_onevsrest <- rbindlist(lapply(test_result$onevsrest_list, make_expanded_row, comp_type = "onevsrest"), fill = TRUE)
      expanded <- rbindlist( list(expanded_pairwise, expanded_onevsrest), fill = TRUE)
      if(!is.null(expanded) && nrow(expanded) > 0L){ test_dt <- rbind(test_dt, expanded, fill = TRUE) }
    }
  }
  summary_idx <- which(test_dt$comparison_type == "summary")
  all_p <- test_dt[[p_col]]
  summary_p <- all_p[summary_idx]
  valid_idx <- summary_idx[!is.na(summary_p)]
  adjusted_bonf <- rep(NA_real_, nrow(test_dt))
  adjusted_fdr <- rep(NA_real_, nrow(test_dt))
  adjusted_bonf[valid_idx] <- p.adjust(summary_p[!is.na(summary_p)], method = "bonferroni")
  adjusted_fdr[valid_idx] <- p.adjust(summary_p[!is.na(summary_p)], method = "BH")
  test_dt[,p_bonferroni := adjusted_bonf]
  test_dt[,p_fdr := adjusted_fdr]
  test_dt[,sig_unadjusted := all_p < 0.05]
  test_dt[,sig_bonferroni := p_bonferroni < 0.05]
  test_dt[,sig_fdr := p_fdr < 0.05]
  return(test_dt[])
}

#' Run \code{\link{categorical_compare}} across every binary grouping variable
#'
#' For each row of \code{CatTable} (as produced by
#' \code{\link{CategoricalVariableCreate}}), runs
#' \code{\link{categorical_compare}} of that grouping variable against every
#' column in \code{CategoricalCols}, writing one Excel file per grouping
#' variable to \code{Fpath} and returning the compiled results.
#'
#' @param dt A data.table containing the columns named in \code{CatTable} and
#'   \code{CategoricalCols}.
#' @param CategoricalCols Character vector of categorical column names to
#'   test against each grouping variable.
#' @param Type One of \code{"fisher"} or \code{"chi"}, passed through to
#'   \code{\link{categorical_compare}}.
#' @param CatTable A data.table of binary grouping variables, as produced by
#'   \code{\link{CategoricalVariableCreate}}.
#' @param Fpath Directory to write one \code{CategoricalComparisons_*.xlsx}
#'   file per grouping variable.
#' @return A data.table combining the \code{\link{categorical_compare}}
#'   results across every row of \code{CatTable}.
#' @export
Categorical_Compare_compile <- function(dt, CategoricalCols, Type = "fisher", CatTable, Fpath){
  CompiledAll <- data.table()
  if (!dir.exists(Fpath)) stop("Fpath does not exist: ", Fpath)
  for(i in seq_len(nrow(CatTable))){
    grouping_col <- CatTable$CategoricalVariable[i]
    keys <- setdiff(CategoricalCols, grouping_col)
    CatResults <- categorical_compare(DT=dt, compair=CatTable$NameOrientation[i], compair_column=grouping_col, keys=keys, Stattest = Type)
    Fname <- paste("CategoricalComparisons", Type, paste(gsub("-", "", CatTable$NameOrientation[i]), ".xlsx", sep = ""), sep = "_")
    openxlsx::write.xlsx(list(CategoricalComparisons = CatResults), file = file.path(Fpath, Fname), colWidths = "auto" )
    CompiledAll <- rbind(CompiledAll, CatResults)
    }
  return(CompiledAll)
}
