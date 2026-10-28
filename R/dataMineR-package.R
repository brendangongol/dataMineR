#' dataMineR: Clinical Data Mining and Comparison Tools
#'
#' Functions for grouped Table 1 summaries, regression screening, continuous
#' and categorical group comparisons, correlation and agreement analyses,
#' survival analysis, mediation analysis, dimensionality reduction, power
#' analysis, bootstrap predictor-stability screening, survey-weighted LASSO
#' variable selection (\code{\link{LassoCompile}}), and propensity score
#' matching, inverse probability weighting and entropy balancing
#' (\code{\link{PSRegression}}, \code{\link{IPTWRegression}},
#' \code{\link{EBRegression}}).
#'
#' See \code{demo("dataMineRWorkflow", package = "dataMineR")} for an
#' end-to-end example built on the package's synthetic example datasets
#' (\code{\link{dfMigrant}}, \code{\link{copd_data}}, \code{\link{dream_data}},
#' \code{\link{peru_data}}, \code{\link{egs_misc_1621}}, and others).
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import data.table
#' @import ggplot2
#' @importFrom ggpubr ggdensity theme_pubr
#' @importFrom stats as.dist as.formula anova binomial chisq.test coef
#'   complete.cases cor dist fisher.test glm hclust lm median
#'   p.adjust pchisq predict prcomp quantile relevel setNames
#'   shapiro.test wilcox.test weighted.mean na.omit ave confint
#'   delete.response fitted gaussian lm.wfit model.frame model.matrix
#'   model.response na.fail plogis pt qlogis qnorm qt quasibinomial
#'   quasipoisson reformulate sd terms vcov
#' @importFrom survival coxph survdiff survfit Surv
#' @importFrom utils combn txtProgressBar setTxtProgressBar head tail
#' @importFrom pROC roc coords
#' @importFrom scales hue_pal
#' @importFrom svyVarSel wlasso.plot
## usethis namespace: end
NULL

utils::globalVariables(c(
  ".", ".N", ".SD", ".SDcols", ":=", ".data", ".dep_source", ".ind_source", ".SigFlag",
  "dependent", "independent", "Significance", "p_fdr", "names", "formula",
  "comparison", "parameter", "PlotValue", "Row", "Col", "Value", "group",
  "model", "value", "Feature", "CategoricalVariable", "binomial",
  "Pr(>|t|)", "Pr(>|z|)", "Pr(>F)",
  "RowOrder", "MergeColumn", "MasterColumn", "Characteristic", "IsNARow",
  "smd_val", "SMD_imbalanced", "P_value", "p_value", "Effect_95CI",
  "N", "N_pair", "N_x", "N_y", "X", "Y", "medN", "minN", "maxN",
  "PC1", "PC2", "Group", "Label", "UMAP1", "UMAP2", "cutoff", "AUC",
  "odds_ratio_CI_low", "odds_ratio_CI_high", "odds_ratio_str",
  "risk_ratio_CI_low", "risk_ratio_CI_high", "risk_ratio_str",
  "p_bonferroni", "sig_unadjusted", "sig_bonferroni", "sig_fdr",
  "p_fdr_overall", "p_bonferroni_overall", "sig_fdr_overall",
  "StabilityFraction", "TimesSignificant", "SurvTime", "SurvStatus",
  "BonfSignificance", "BHSignificance", "FDR", "p.adjust",
  ".column", ".ps_fit_weight", ".row_label", "ModelID", "analysis_id", "coefficient",
  "column", "correlation", "covariate", "display_label", "joint_duplicate", "label",
  "lambda", "logit_ps", "near_duplicate", "predictor", "reason", "role", "size",
  "std_coef", "variable", "w"
))
