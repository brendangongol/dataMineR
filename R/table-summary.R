#' Build a grouped "Table 1" with SMD, p-values, and effect sizes
#'
#' Produces a grouped descriptive summary table (continuous variables as
#' N/mean/median, categorical as n/percent) from a data.frame or data.table,
#' with an overall column plus one column per level of \code{GroupBy}.
#' Optionally adds standardized mean differences (SMD), p-values (Wilcoxon/
#' Kruskal-Wallis for continuous, Fisher/chi-squared for categorical), and
#' effect estimates with confidence intervals (Hodges-Lehmann median
#' difference for continuous, odds ratio for categorical).
#'
#' @param TableDF A data.frame or data.table to summarize.
#' @param XLSXPath Output path for the Excel workbook (required if
#'   \code{WriteTable = TRUE}).
#' @param GroupBy Optional column name to stratify the table by. If
#'   \code{NULL} or not found in \code{TableDF}, an "Overall" column only is
#'   produced.
#' @param ExcludeSeparationGroup If \code{TRUE}, drop \code{GroupBy} itself
#'   from the row list.
#' @param RoundLevel Number of decimal places for summary statistics.
#' @param QuantileType Quantile algorithm passed to \code{\link[stats]{quantile}}
#'   (3 matches SAS; 7 is R's default).
#' @param group_order Optional character vector giving the display order of
#'   \code{GroupBy} levels (as columns, after "Overall").
#' @param weight_col Optional column name of survey/sampling weights to use
#'   for weighted means, variances, and quantiles. Note p-values and effect
#'   CIs remain unweighted base-R tests.
#' @param continuous_format One of \code{"n_mean_median"}, \code{"mean_sd"},
#'   or \code{"median_iqr"} controlling how continuous rows are summarized.
#' @param labels Optional named character vector mapping column names to
#'   display labels, e.g. \code{c(mortality30day = "30-day mortality")}.
#' @param show_missing One of \code{"inline"}, \code{"bottom"}, or
#'   \code{"none"} controlling how missing-value counts are reported.
#' @param smd Logical; add a standardized mean difference column comparing
#'   the first two non-overall groups.
#' @param smd_ref Optional group name to use as the SMD reference (defaults
#'   to the first non-overall group).
#' @param p_values Logical; add a p-value column.
#' @param p_value_digits Number of digits used when formatting p-values.
#' @param continuous_p_test One of \code{"wilcox"} or \code{"kruskal"}.
#' @param categorical_p_test One of \code{"fisher"} or \code{"chi"}.
#' @param effect_ci Logical; add an effect-estimate-with-95\%-CI column
#'   (Hodges-Lehmann difference for continuous rows, odds ratio for
#'   categorical rows).
#' @param NormalizeData Logical; if \code{TRUE}, relabel raw values (e.g.
#'   \code{"Checked"} -> \code{"Yes"}) before summarizing, per
#'   \code{NormalizeDTKey}.
#' @param NormalizeDTKey Optional data.table of relabeling rules used when
#'   \code{NormalizeData = TRUE}, with columns \code{ColToNormalize},
#'   \code{NormalizeColumn}, \code{NormalizeValue}, \code{ColToNormalizeValue}.
#' @param footnotes Logical; append explanatory footnotes to the table.
#' @param WriteTable Logical; if \code{TRUE}, write the table to
#'   \code{XLSXPath}.
#' @return A data.table holding the compiled summary table; also written to
#'   \code{XLSXPath} when \code{WriteTable = TRUE}.
#' @examples
#' \dontrun{
#' load(system.file("extdata", "dfMigrant.rda", package = "dataMineR"))
#' InitialTablebg(TableDF = dfMigrant, GroupBy = "Migrant", p_values = TRUE)
#' }
#' @export
InitialTablebg <- function(TableDF, XLSXPath = NULL, GroupBy = NULL, ExcludeSeparationGroup = FALSE,
                           RoundLevel = 1, QuantileType = 3, group_order = NULL, weight_col = NULL,
                           continuous_format = c("n_mean_median", "mean_sd", "median_iqr"),
                           labels = NULL, show_missing = c("inline", "bottom", "none"), smd = TRUE,
                           smd_ref = NULL, p_values = FALSE, p_value_digits = 4, continuous_p_test = c("wilcox", "kruskal"),
                           categorical_p_test = c("fisher", "chi"), effect_ci = FALSE, NormalizeData = FALSE,
                           NormalizeDTKey = NULL, footnotes = TRUE, WriteTable = FALSE){
  continuous_format <- match.arg(continuous_format)
  show_missing <- match.arg(show_missing)
  continuous_p_test <- match.arg(continuous_p_test)
  categorical_p_test <- match.arg(categorical_p_test)
  TableDF <- as.data.table(TableDF)
  #### Validate weight_col ####
  if(!is.null(weight_col) && !weight_col %in% names(TableDF)){ stop(sprintf("weight_col '%s' not found in TableDF.", weight_col)) }
  if (!is.null(weight_col)) {
    weights <- TableDF[[weight_col]]
    if (!is.numeric(weights) || any(!is.finite(weights[!is.na(weights)])) || any(weights[!is.na(weights)] < 0) || sum(weights, na.rm = TRUE) <= 0) {
      stop("weight_col must contain finite, non-negative numeric weights with a positive sum.")
    }
  }
  if(!is.null(weight_col) && (p_values || effect_ci)){
    warning(paste("InitialTablebg uses weight_col for descriptive summaries and SMD,",
                  "but p-values and effect confidence intervals are unweighted base-R tests.")) }
  ##########################################################################
  #### Stage 0: split into groups and record group N for column headers ####
  ##########################################################################
  if (!is.null(GroupBy) && GroupBy %in% names(TableDF)) {
    n_na <- sum(is.na(TableDF[[GroupBy]]))
    if(n_na > 0){
      message(sprintf("Excluding %d records with NA in '%s'.", n_na, GroupBy)) }
    TableDF <- TableDF[!is.na(TableDF[[GroupBy]]), ]
    Separate <- sort(unique(TableDF[[GroupBy]]))
    Separate <- Separate[!is.na(Separate)]
    if(!is.null(group_order)){
      missing_groups <- setdiff(group_order, as.character(Separate))
      if(length(missing_groups) > 0){
        warning(sprintf("group_order values not found in '%s': %s", GroupBy, paste(missing_groups, collapse = ", "))) }
      ordered_groups <- group_order[group_order %in% as.character(Separate)]
      remaining_groups <- as.character(Separate)[!as.character(Separate) %in% ordered_groups]
      Separate <- c(ordered_groups, remaining_groups)
    }
    TableList <- lapply(Separate, function(x) TableDF[TableDF[[GroupBy]] == x, ])
    names(TableList) <- as.character(Separate)
    TableList <- c(list(Overall = TableDF), TableList)
  } else {
    message("GroupBy is NULL or not found. Producing an Overall column only.")
    TableList <- list(Overall = TableDF)
  }
  #### Append group N to column header names ####
  group_n    <- vapply(TableList, nrow, integer(1))
  group_labels_hdr <- sprintf("%s (n=%s)", names(TableList), formatC(group_n, format = "d", big.mark = ","))
  names(group_labels_hdr) <- names(TableList)
  Classes <- vapply(TableList[["Overall"]], function(x) class(x)[1], character(1))
  if(ExcludeSeparationGroup && !is.null(GroupBy)){
    Classes <- Classes[names(Classes) != GroupBy] }
  if(!is.null(weight_col)){
    Classes <- Classes[names(Classes) != weight_col] }
  int_cols <- names(Classes)[Classes == "integer"]
  if(length(int_cols) > 0){
    warning(sprintf(paste("Integer-class columns will be summarized as CATEGORICAL",
                          "(one row per distinct value): %s. Convert to numeric",
                          "(e.g. ColClassConvert with ColClass = 'numeric') for",
                          "continuous treatment."),
                    paste(int_cols, collapse = ", "))) }
  #### SMD reference group ####
  non_overall <- names(TableList)[names(TableList) != "Overall"]
  if(smd && length(non_overall) >= 2) {
    smd_ref_grp <- if (!is.null(smd_ref) && smd_ref %in% non_overall)
      smd_ref else non_overall[1]
    smd_compare  <- non_overall[non_overall != smd_ref_grp]
  } else {
    smd <- FALSE
    smd_ref_grp <- NULL
    smd_compare <- NULL
  }
  ##########################
  #### Weighted helpers ####
  ##########################
  w_mean <- function(x, w){
    ok <- !is.na(x) & !is.na(w)
    if(sum(ok) == 0){ return(NA_real_) }
    weighted.mean(x[ok], w[ok])
  }
  w_var <- function(x, w){
    ok <- !is.na(x) & !is.na(w)
    if(sum(ok) < 2 || sum(w[ok]) <= 1){ return(NA_real_) }
    xw  <- w_mean(x, w)
    sum(w[ok] * (x[ok] - xw)^2) / (sum(w[ok]) - 1)
  }
  w_quantile <- function(x, w, probs){
    ok <- !is.na(x) & !is.na(w)
    if(sum(ok) == 0){ return(rep(NA_real_, length(probs))) }
    if(is.null(weight_col)){
      return(as.numeric(stats::quantile(x[ok], probs, type = QuantileType, names = FALSE))) }
    ord <- order(x[ok])
    xs <- x[ok][ord]
    ws <- w[ok][ord]
    cw <- cumsum(ws) / sum(ws)
    vapply(probs, function(p) xs[which(cw >= p)[1]], numeric(1))
  }
  get_weights <- function(tdt) {
    if(!is.null(weight_col)){ tdt[[weight_col]] } else { rep(1, nrow(tdt)) }
  }
  #####################
  #### SMD helpers ####
  #####################
  smd_continuous <- function(tdt_ref, tdt_cmp, col) {
    w_r <- get_weights(tdt_ref); x_r <- tdt_ref[[col]]
    w_c <- get_weights(tdt_cmp); x_c <- tdt_cmp[[col]]
    m_r <- w_mean(x_r, w_r); m_c <- w_mean(x_c, w_c)
    v_r <- w_var(x_r, w_r); v_c <- w_var(x_c, w_c)
    if(anyNA(c(m_r, m_c, v_r, v_c))){ return(NA_real_) }
    pooled_sd <- sqrt((v_r + v_c) / 2)
    if(pooled_sd == 0){ return(0) }
    round((m_r - m_c) / pooled_sd, RoundLevel + 2)
  }
  smd_categorical <- function(tdt_ref, tdt_cmp, col) {
    w_r <- get_weights(tdt_ref); x_r <- as.character(tdt_ref[[col]])
    w_c <- get_weights(tdt_cmp); x_c <- as.character(tdt_cmp[[col]])
    levs <- union(unique(x_r[!is.na(x_r)]), unique(x_c[!is.na(x_c)]))
    if(length(levs) == 0L){ return(NA_real_) }
    denom_r <- sum(w_r[!is.na(x_r)])
    denom_c <- sum(w_c[!is.na(x_c)])
    if(denom_r <= 0 || denom_c <= 0){ return(NA_real_) }
    smds <- vapply(levs, function(lv){
      p_r <- sum(w_r[!is.na(x_r) & x_r == lv])/denom_r
      p_c <- sum(w_c[!is.na(x_c) & x_c == lv])/denom_c
      pp<- (p_r + p_c)/2
      if(pp <= 0 || pp >= 1){ return(0) }
      abs((p_r - p_c)/sqrt(pp*(1-pp)))
    }, numeric(1))
    round(max(smds, na.rm = TRUE), RoundLevel + 2)
  }
  #########################
  #### P-value helpers ####
  #########################
  fmt_p_value <- function(p){
    if(is.na(p)){ return(NA_character_) }
    if(p < 10^(-p_value_digits)){
      return(paste0("<", formatC(10^(-p_value_digits), format = "f", digits = p_value_digits))) }
    formatC(p, format = "f", digits = p_value_digits)
  }
  p_value_continuous <- function(col){
    if(is.null(GroupBy) || length(non_overall) < 2){ return(NA_real_) }
    x <- suppressWarnings(as.numeric(TableList[[non_overall[1]]][[col]]))
    y <- suppressWarnings(as.numeric(TableList[[non_overall[2]]][[col]]))
    x <- x[!is.na(x)]
    y <- y[!is.na(y)]
    if(length(non_overall) == 2L && (length(x) < 2 || length(y) < 2)){ return(NA_real_) }
    if(length(non_overall) == 2L && length(unique(x)) == 1 && length(unique(y)) == 1){ return(NA_real_) }
    p <- tryCatch({
      if(continuous_p_test == "wilcox" && length(non_overall) == 2L) {
        stats::wilcox.test(x, y, conf.int = TRUE, paired = FALSE)$p.value
      } else {
        test_dt <- rbindlist(lapply(non_overall, function(grp) {
          data.table(group = grp,
                     value = suppressWarnings(as.numeric(TableList[[grp]][[col]])) )
        }), fill = TRUE)
        test_dt <- test_dt[!is.na(value) & !is.na(group)]
        if(test_dt[, uniqueN(group)] < 2L){ return(NA_real_) }
        stats::kruskal.test(value ~ group, data = test_dt)$p.value
      }
    }, error = function(e) NA_real_)
    p
  }
  p_value_categorical <- function(col){
    if(is.null(GroupBy) || length(non_overall) < 2){ return(NA_real_) }
    test_dt <- rbindlist(lapply(non_overall, function(grp){
      data.table(group = grp, value = as.character(TableList[[grp]][[col]]) )
    }), fill = TRUE)
    test_dt <- test_dt[!is.na(value) & !is.na(group)]
    if(test_dt[, uniqueN(group)] < 2L || test_dt[, uniqueN(value)] < 2L){ return(NA_real_) }
    tab <- table(test_dt$group, test_dt$value)
    p <- tryCatch({
      if(categorical_p_test == "fisher") {
        if(nrow(tab) == 2L && ncol(tab) == 2L) {
          stats::fisher.test(tab, conf.level = 0.95)$p.value
        } else {
          stats::fisher.test(tab, simulate.p.value = TRUE, B = 2000)$p.value
        }
      } else {
        stats::chisq.test(tab)$p.value
      }
    }, error = function(e) NA_real_)
    p
  }
  #############################################################################
  #### Effect estimate + exact 95% CI helpers (two-group comparisons only) ####
  #############################################################################
  fmt_effect <- function(est, lo, hi){
    if(anyNA(c(est, lo, hi))){ return(NA_character_) }
    sprintf("%s (%s, %s)", round(est, RoundLevel + 1), round(lo, RoundLevel + 1), round(hi, RoundLevel + 1))
  }
  effect_ci_continuous <- function(col){
    if(is.null(GroupBy) || length(non_overall) != 2L){ return(NA_character_) }
    x <- suppressWarnings(as.numeric(TableList[[non_overall[1]]][[col]])); x <- x[!is.na(x)]
    y <- suppressWarnings(as.numeric(TableList[[non_overall[2]]][[col]])); y <- y[!is.na(y)]
    if(length(x) < 2 || length(y) < 2){ return(NA_character_) }
    tryCatch({
      wt <- suppressWarnings(stats::wilcox.test(x, y, conf.int = TRUE, paired = FALSE))
      fmt_effect(unname(wt$estimate), wt$conf.int[1], wt$conf.int[2])
    }, error = function(e) NA_character_)
  }
  effect_ci_categorical <- function(col) {
    if(is.null(GroupBy) || length(non_overall) != 2L){ return(NA_character_) }
    test_dt <- rbindlist(lapply(non_overall, function(grp) {
      data.table(group = grp, value = as.character(TableList[[grp]][[col]])) }), fill = TRUE)
    test_dt <- test_dt[!is.na(value) & !is.na(group)]
    tab <- table(factor(test_dt$group, levels = non_overall), test_dt$value)
    if(nrow(tab) != 2L || ncol(tab) != 2L){ return(NA_character_) }
    tryCatch({
      ft <- stats::fisher.test(tab, conf.level = 0.95)
      fmt_effect(unname(ft$estimate), ft$conf.int[1], ft$conf.int[2])
    }, error = function(e) NA_character_)
  }
  ##############################################
  #### Compute per-group summary statistics ####
  ##############################################
  overall_levels <- lapply(names(Classes), function(col) {
    if(Classes[[col]] == "numeric"){ return(NULL) }
    x <- as.character(TableList[["Overall"]][[col]])
    sort(unique(x[!is.na(x)])) })
  names(overall_levels) <- names(Classes)
  CompleteTableList <- lapply(seq_along(TableList), function(a){
    TDT <- TableList[[a]]
    grp_name <- names(TableList)[a]
    W <- get_weights(TDT)
    CompleteTable <- data.table()
    for(b in seq_along(Classes)){
      col_name <- names(Classes)[b]
      col_class <- Classes[[b]]
      col_vals <- TDT[[col_name]]
      N_total <- nrow(TDT)
      N_na <- sum(is.na(col_vals))
      ok <- !is.na(col_vals)
      if(col_class == "numeric"){
        N_valid <- sum(ok)
        wok <- W[ok]; xok <- col_vals[ok]
        summary_str <- switch(continuous_format,
                              "n_mean_median" = {
                                Mea <- round(w_mean(xok, wok), RoundLevel)
                                Med <- round(w_quantile(xok, wok, 0.5), RoundLevel)
                                Q1 <- round(w_quantile(xok, wok, 0.25), RoundLevel)
                                Q3 <- round(w_quantile(xok, wok, 0.75), RoundLevel)
                                sprintf("%d, %s, %s (%s, %s)", N_valid, Mea, Med, Q1, Q3)
                              },
                              "mean_sd" = {
                                Mea <- round(w_mean(xok, wok), RoundLevel)
                                SD <- round(sqrt(w_var(xok, wok)), RoundLevel)
                                sprintf("%s +/- %s", Mea, SD)
                              },
                              "median_iqr" = {
                                Med <- round(w_quantile(xok, wok, 0.5), RoundLevel)
                                Q1 <- round(w_quantile(xok, wok, 0.25), RoundLevel)
                                Q3 <- round(w_quantile(xok, wok, 0.75), RoundLevel)
                                sprintf("%s (%s, %s)", Med, Q1, Q3)
                              } )
        na_row <- if(N_na > 0 && show_missing != "none"){
          data.table(Characteristic = "NNA", MasterColumn = col_name, Value = as.character(N_na), IsNARow = TRUE)
        } else { NULL }
        title_row <- data.table(Characteristic = paste0(col_name, ";Title"),
                                MasterColumn = col_name,
                                Value = "",
                                IsNARow = FALSE)
        body_row  <- data.table(Characteristic = col_name,
                                MasterColumn = col_name,
                                Value = summary_str,
                                IsNARow = FALSE)
        TempDT <- if(!is.null(na_row)){
          rbind(title_row, body_row, na_row)
        } else{ rbind(title_row, body_row) }
      } else {
        x_char <- as.character(col_vals)
        levs <- overall_levels[[col_name]]
        w_total <- sum(W[ok])
        counts <- vapply(levs, function(lv){ sum(W[ok & x_char == lv]) }, numeric(1))
        pcts <- if(w_total > 0){ round(counts/w_total*100, RoundLevel) } else { rep(0, length(levs)) }
        raw_n <- vapply(levs, function(lv) sum(ok & x_char == lv), integer(1))
        vp <- sprintf("%d (%s%%)", raw_n, pcts)
        body_rows <- data.table(Characteristic = c(paste0(col_name, ";Title"), levs),
                                MasterColumn = col_name, Value = c("", vp), IsNARow = FALSE)
        na_row <- if(N_na > 0 && show_missing != "none") {
          na_pct <- round(N_na/N_total*100, RoundLevel)
          data.table(Characteristic = "NNA", MasterColumn = col_name,
                     Value = sprintf("%d (%s%%)", N_na, na_pct),
                     IsNARow = TRUE)
        } else { NULL }
        TempDT <- if(!is.null(na_row)){
          rbind(body_rows, na_row)
        } else { body_rows }
      }
      setnames(TempDT, "Value", grp_name)
      CompleteTable <- rbind(CompleteTable, TempDT, fill = TRUE)
    }
    CompleteTable
  })
  names(CompleteTableList) <- names(TableList)
  ####################################
  #### Normalise selected columns ####
  ####################################
  if(NormalizeData && !is.null(NormalizeDTKey)) {
    CompleteTableList <- lapply(seq_along(CompleteTableList), function(a){
      tempTab  <- CompleteTableList[[a]]
      ColName  <- names(CompleteTableList)[a]
      for(b in seq_len(nrow(NormalizeDTKey))){
        norm_col <- NormalizeDTKey$NormalizeColumn[b]
        norm_val <- NormalizeDTKey$NormalizeValue[b]
        target_col <- NormalizeDTKey$ColToNormalize[b]
        target_val <- NormalizeDTKey$ColToNormalizeValue[b]
        norm_row <- tempTab[MasterColumn == norm_col & Characteristic == norm_val]
        if(nrow(norm_row) == 0L){ next }
        NormValue  <- suppressWarnings(as.numeric(gsub(" .+", "", norm_row[[ColName]][1])))
        if(is.na(NormValue) || NormValue == 0){ next }
        idx <- which(tempTab$MasterColumn == target_col & tempTab$Characteristic == target_val)
        if(length(idx) == 0L){ next }
        raw_n <- suppressWarnings(as.numeric(gsub(" .+", "", tempTab[[ColName]][idx[1]])))
        if(is.na(raw_n)){ next }
        normalized_str <- sprintf("%s Normalized: %d/%d %s%%",
                                  tempTab[[ColName]][idx[1]], raw_n, NormValue,
                                  round(raw_n / NormValue * 100, RoundLevel))
        tempTab[idx[1], (ColName) := normalized_str]
      }
      tempTab
    })
    names(CompleteTableList) <- names(TableList)
  }
  ##############################################
  #### Merge all groups into one wide table ####
  ##############################################
  CompleteTable <- CompleteTableList[["Overall"]]
  CompleteTable[, RowOrder := seq_len(.N)]
  CompleteTable[, MergeColumn := paste(MasterColumn, Characteristic, sep = ";")]
  for(a in seq_along(CompleteTableList)[-1]){
    grp_name <- names(CompleteTableList)[a]
    AddTable <- copy(CompleteTableList[[a]])
    AddTable[, MergeColumn := paste(MasterColumn, Characteristic, sep = ";")]
    missing_keys <- setdiff(CompleteTable$MergeColumn, AddTable$MergeColumn)
    if(length(missing_keys) > 0){
      Missing1 <- CompleteTable[MergeColumn %in% missing_keys, .(MergeColumn, MasterColumn, Characteristic, IsNARow)]
      Missing1[, (grp_name) := "0 (0%)"]
      AddTable <- rbind(AddTable, Missing1, fill = TRUE)
    }
    CompleteTable <- merge(CompleteTable, AddTable[, !c("Characteristic", "MasterColumn", "IsNARow"), with = FALSE], by = "MergeColumn", all.x = TRUE )
    CompleteTable <- CompleteTable[order(RowOrder)]
  }
  #####################
  #### SMD columns ####
  #####################
  if(smd){
    ref_data <- TableList[[smd_ref_grp]]
    for(cmp_grp in smd_compare){
      cmp_data <- TableList[[cmp_grp]]
      smd_col <- sprintf("SMD_%s_vs_%s", cmp_grp, smd_ref_grp)
      smd_vals <- vapply(names(Classes), function(col){
        if(!is.null(GroupBy) && col == GroupBy){ return(NA_real_) }
        cls <- Classes[[col]]
        if(cls == "numeric"){ smd_continuous(ref_data, cmp_data, col) 
          } else{ smd_categorical(ref_data, cmp_data, col) }
      }, numeric(1))
      smd_dt <- data.table( MasterColumn = names(Classes), smd_val = smd_vals )
      CompleteTable[, (smd_col) := NA_real_]
      for(col in names(Classes)){
        idx <- which(CompleteTable$MasterColumn == col & !grepl(";Title$", CompleteTable$MergeColumn) & CompleteTable$Characteristic != "NNA")[1]
        if(!is.na(idx)){ set(CompleteTable, idx, smd_col, smd_dt[MasterColumn == col, smd_val]) }
      }
    }
    smd_cols <- grep("^SMD_", names(CompleteTable), value = TRUE)
    CompleteTable[, SMD_imbalanced := apply( .SD, 1, function(r) any(!is.na(r) & abs(r) > 0.1)), .SDcols = smd_cols]
  }
  ########################
  #### P-value column ####
  ########################
  if(p_values){
    p_vals <- vapply(names(Classes), function(col) {
      if(!is.null(GroupBy) && col == GroupBy){ return(NA_real_) }
      cls <- Classes[[col]]
      if(cls == "numeric"){ p_value_continuous(col)
      } else { p_value_categorical(col) }
    }, numeric(1))
    p_dt <- data.table(
      MasterColumn = names(Classes),
      p_value = vapply(p_vals, fmt_p_value, character(1))
    )
    CompleteTable[, P_value := NA_character_]
    for(col in names(Classes)){
      idx <- which(CompleteTable$MasterColumn == col & !grepl(";Title$", CompleteTable$MergeColumn) & CompleteTable$Characteristic != "NNA")[1]
      if(!is.na(idx)){ set(CompleteTable, idx, "P_value", p_dt[MasterColumn == col, p_value]) }
    }
  }
  ##################################################
  #### Effect estimate with exact 95% CI column ####
  ##################################################
  if(effect_ci){
    eff_vals <- vapply(names(Classes), function(col) {
      if(!is.null(GroupBy) && col == GroupBy){ return(NA_character_) }
      if(Classes[[col]] == "numeric"){ effect_ci_continuous(col)
      } else { effect_ci_categorical(col) }
    }, character(1))

    CompleteTable[, Effect_95CI := NA_character_]
    for(col in names(Classes)){
      idx <- which(CompleteTable$MasterColumn == col & !grepl(";Title$", CompleteTable$MergeColumn) & CompleteTable$Characteristic != "NNA")[1]
      if(!is.na(idx)){ set(CompleteTable, idx, "Effect_95CI", eff_vals[[col]]) }
    }
  }
  #########################################################
  #### Apply variable labels and rename column headers ####
  #########################################################
  if(!is.null(labels)){
    for(col in names(labels)){
      idx <- which(CompleteTable$MasterColumn == col & grepl(";Title$", CompleteTable$MergeColumn))
      if(length(idx) > 0){
        set(CompleteTable, idx, "Characteristic", paste0(labels[[col]], ";Title")) }
    }
  }
  title_merge_keys <- CompleteTable$MergeColumn[grepl(";Title$", CompleteTable$MergeColumn)]
  CompleteTable[, Characteristic := gsub(";Title$|;.+", "", Characteristic)]
  for(grp in names(group_labels_hdr)){
    if(grp %in% names(CompleteTable)){setnames(CompleteTable, grp, group_labels_hdr[[grp]]) }
  }
  lead_cols <- c("RowOrder", "MergeColumn", "MasterColumn", "Characteristic", "IsNARow")
  lead_cols <- intersect(lead_cols, names(CompleteTable))
  group_cols <- unname(group_labels_hdr[names(group_labels_hdr) %in% names(TableList)])
  group_cols <- intersect(group_cols, names(CompleteTable))
  p_cols <- intersect(c("P_value", "Effect_95CI"), names(CompleteTable))
  smd_cols <- grep("^SMD_", names(CompleteTable), value = TRUE)
  other_cols <- setdiff(names(CompleteTable), c(lead_cols, group_cols, p_cols, smd_cols))
  setcolorder(CompleteTable, c(lead_cols, group_cols, p_cols, smd_cols, other_cols))
  ################################################################
  #### show_missing = "bottom". move NNA rows to end of table ####
  ################################################################
  if(show_missing == "bottom" && any(CompleteTable$IsNARow)){
    na_block <- CompleteTable[IsNARow == TRUE]
    na_block[, Characteristic := sprintf("Missing: %s", ifelse(MasterColumn %in% names(labels), labels[MasterColumn], MasterColumn))]
    CompleteTable <- rbind(CompleteTable[IsNARow == FALSE], na_block)
    CompleteTable[, RowOrder := seq_len(.N)]
  }
  #################################
  #### Auto-generate footnotes ####
  #################################
  footnote_str <- NULL
  if(footnotes){
    fn_parts <- character(0)
    fn_parts <- c(fn_parts, switch(continuous_format,
                                   "n_mean_median" = "Continuous variables: N, mean, median (Q1, Q3).",
                                   "mean_sd" = "Continuous variables: mean +/- SD.",
                                   "median_iqr" = "Continuous variables: median (Q1, Q3)." ))
    fn_parts <- c(fn_parts, "Categorical variables: N (%).")
    if (!is.null(weight_col)) {
      fn_parts <- c(fn_parts, sprintf("Descriptive statistics and SMD are survey-weighted using '%s'.", weight_col))
      if(p_values || effect_ci){
        fn_parts <- c(fn_parts, "P-values and effect confidence intervals are unweighted base-R tests.") }
    }
    if(show_missing == "inline"){ fn_parts <- c(fn_parts, "NNA = number of records with missing values.") }
    if(show_missing == "bottom"){ fn_parts <- c(fn_parts, "Missing value counts shown at the bottom of the table.") }
    if(smd){ 
      fn_parts <- c(fn_parts, sprintf("SMD = standardized mean difference (reference group: %s). |SMD| > 0.1 indicates meaningful imbalance.", smd_ref_grp)) }
    if(p_values){
      fn_parts <- c(fn_parts,
                    sprintf("P-values compare non-overall '%s' groups. Continuous variables: %s%s. Categorical variables: %s test.",
                            GroupBy, continuous_p_test, ifelse(continuous_p_test == "wilcox", " for two groups and Kruskal-Wallis when more than two groups are present", ""), categorical_p_test))  }
    if(effect_ci){
      eff_groups <- if(length(non_overall) == 2L){
        sprintf(" (%s minus %s / %s vs %s)", non_overall[1], non_overall[2], non_overall[1], non_overall[2])
      } else { ""
      fn_parts <- c(fn_parts, sprintf(paste("Effect (95%% CI): Hodges-Lehmann median difference with exact CI",
                                            "(Wilcoxon) for continuous variables; odds ratio with exact CI (Fisher)",
                                            "for binary variables%s; blank for variables with more than two levels."), eff_groups))}
    }
    if(is.null(weight_col)){
      fn_parts <- c(fn_parts, sprintf("Quartiles computed using type=%d%s.", QuantileType, ifelse(QuantileType == 3, " (SAS-compatible)", "")))
    } else { fn_parts <- c(fn_parts, "Quartiles computed via weighted cumulative-distribution method.") }
    footnote_str <- paste(fn_parts, collapse = " ")
    attr(CompleteTable, "footnotes") <- footnote_str 
  }
  ##############################
  #### Write Excel workbook ####
  ##############################
  if(WriteTable){
    if(is.null(XLSXPath)){ stop("XLSXPath must be provided when WriteTable = TRUE.") }
    export_cols <- setdiff(names(CompleteTable), c("RowOrder", "MergeColumn", "MasterColumn", "IsNARow", "SMD_imbalanced"))
    export_dt <- CompleteTable[, export_cols, with = FALSE]
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "Summary Table")
    openxlsx::writeData(wb, "Summary Table", export_dt, startRow = 1, startCol = 1)
    n_data_rows <- nrow(export_dt)
    n_cols <- ncol(export_dt)
    #### Header style ####
    openxlsx::addStyle(wb, "Summary Table",
             openxlsx::createStyle(fontName = "Calibri", fontSize = 11, textDecoration = "bold", halign = "center", border = "Bottom"),
             rows = 1, cols = seq_len(n_cols), gridExpand = TRUE)
    #### Variable title rows: grey background ####
    title_rows <- which(CompleteTable$MergeColumn %in% title_merge_keys) + 1L
    if(length(title_rows) > 0){
      openxlsx::addStyle(wb, "Summary Table", openxlsx::createStyle(fgFill = "#D3D3D3"), rows = title_rows, cols = seq_len(n_cols), gridExpand = TRUE)
      openxlsx::addStyle(wb, "Summary Table", openxlsx::createStyle(fgFill = "#D3D3D3", textDecoration = "bold"), rows = title_rows, cols = 1L, gridExpand = TRUE)
    }
    if(smd){
      smd_export_cols <- grep("^SMD_", names(export_dt), value = TRUE)
      for(sc in smd_export_cols) {
        sc_idx <- which(names(export_dt) == sc)
        smd_vals_col <- CompleteTable[[sc]]
        bold_rows <- which(!is.na(smd_vals_col) & abs(smd_vals_col) > 0.1) + 1L
        if (length(bold_rows) > 0){
          openxlsx::addStyle(wb, "Summary Table", openxlsx::createStyle(textDecoration = "bold"), rows = bold_rows, cols = sc_idx, gridExpand = FALSE, stack = TRUE) }
      }
    }
    openxlsx::setColWidths(wb, "Summary Table", cols = seq_len(n_cols), widths = "auto")
    #### Footnote row ####
    if (footnotes && !is.null(footnote_str)) {
      fn_row <- n_data_rows + 2L
      openxlsx::writeData(wb, "Summary Table", footnote_str, startRow = fn_row, startCol = 1)
      openxlsx::addStyle(wb, "Summary Table", openxlsx::createStyle(fontSize = 9, fontColour = "#555555", wrapText = TRUE), rows = fn_row, cols = 1L)
      openxlsx::mergeCells(wb, "Summary Table", cols = seq_len(n_cols), rows = fn_row)
    }
    openxlsx::saveWorkbook(wb, XLSXPath, overwrite = TRUE)
    message(sprintf("Workbook saved: %s", XLSXPath))
  }
  return(CompleteTable[])
}
