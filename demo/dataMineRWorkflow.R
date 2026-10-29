################################################################################
################################################################################
#### dataMineRWorkflow.R #######################################################
################################################################################
################################################################################

#################################
#### Environmental variables ####
#################################
demo_packages <- c(
  "data.table", "ggplot2", "ggpubr", "openxlsx", "effectsize", "ggfortify",
  "umap", "GGally", "pwr", "lsr", "pwrss", "blandr", "mediation", "patchwork",
  "Hmisc", "vtable", "seminr", "svyVarSel", "glmnet", "survey", "MatchIt", "WeightIt",
  "cobalt", "sandwich"
)
missing_packages <- demo_packages[!vapply(demo_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) {
  stop("Install the packages required by the full demo: ",
       paste(missing_packages, collapse = ", "))
}
suppressPackageStartupMessages({
  library(dataMineR)
  library(data.table)
  library(ggplot2)
  library(ggpubr)
})
#### Override these with DATAMINER_OUTPUT_DIR / DATAMINER_WORKFLOW_NAME. ####
WritePath <- Sys.getenv("DATAMINER_OUTPUT_DIR", unset = getwd())
DirName   <- Sys.getenv("DATAMINER_WORKFLOW_NAME", unset = "dataMineRWorkflow-results")
WorkflowSeed <- suppressWarnings(as.integer(Sys.getenv("DATAMINER_SEED", unset = "2026")))
if (length(WorkflowSeed) != 1L || is.na(WorkflowSeed)) {
  stop("DATAMINER_SEED must be a valid integer.")
}
set.seed(WorkflowSeed)

ensure_dir <- function(path) {
  if (!dir.exists(path) && !dir.create(path, recursive = TRUE, showWarnings = FALSE)) {
    stop("Unable to create output directory: ", path)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}
OutputDir <- ensure_dir(file.path(WritePath, DirName))

##############################
#### Load example dataset ####
##############################
data(dfMigrant, package = "dataMineR")
dt <- as.data.table(copy(dfMigrant))
GroupBy <- "Migrant"  #### primary two-group comparison variable for the Table 1 / bubble-plot pipeline ####
NonAnalysisCols <- c("MRN", "CSN", "DateofDeath", "DateofAdmission", "DateofDischarge",
                     "MigrationDetails", "TimeofSepsis")

#############################################################
#### remove columns with only one value or all NA values ####
#############################################################
NValuesPerCol <- vapply(dt, data.table::uniqueN, integer(1), na.rm = TRUE)
dt <- dt[, NValuesPerCol > 1, with = FALSE]

#########################################
#### Remove spaces from column names ####
#########################################
colnames(dt) <- gsub(" ", "_", colnames(dt))

#################################################################################################################################################
#################################################################################################################################################
#### analysis ###################################################################################################################################
#################################################################################################################################################
#################################################################################################################################################

######################################################
#### Create a table with all counts or statistics ####
######################################################
TableSummary <- InitialTablebg(TableDF = dt,
                               XLSXPath = file.path(OutputDir, "Table1.xlsx"),
                               GroupBy = GroupBy,
                               ExcludeSeparationGroup = FALSE,
                               RoundLevel = 1,
                               QuantileType = 3,  
                               NormalizeData = FALSE,  
                               WriteTable = TRUE,
                               weight_col = NULL,
                               continuous_format = "n_mean_median",
                               labels = c(InhospitalMortality = "In-hospital mortality",
                                         mortality30day = "30-day mortality",
                                         mortality90day = "90-day mortality"),
                               show_missing = "inline",
                               smd = TRUE,
                               smd_ref = NULL,
                               p_values = TRUE,
                               p_value_digits = 4,
                               continuous_p_test = "wilcox",
                               categorical_p_test = "fisher",
                               effect_ci = TRUE,
                               footnotes = TRUE)

#############################
#### Define column types ####
#############################
#### Adjust these overrides when adapting the workflow to another dataset. ####
ForcedCategoricalCols <- intersect(c("ZIP", "quintiles", "quartiles"), names(dt))
ForcedContinuousCols <- character()
#### Continuous columns ####
numeric_cols <- vapply(dt, is.numeric, logical(1))
if (length(intersect(ForcedCategoricalCols, ForcedContinuousCols)) > 0L) {
  stop("A column cannot appear in both ForcedCategoricalCols and ForcedContinuousCols.")
}
if (any(!ForcedContinuousCols %in% names(dt)[numeric_cols])) {
  stop("ForcedContinuousCols must contain existing numeric columns.")
}
AnalysisCols <- setdiff(names(dt), NonAnalysisCols)
CatTable <- CategoricalVariableCreate(dt = dt, CategoricalDependents = AnalysisCols)
CategoricalDependents <- intersect(colnames(dt), CatTable$CategoricalVariable)
ContinuousCols <- union(
  setdiff(names(dt)[numeric_cols], c(NonAnalysisCols, CategoricalDependents, ForcedCategoricalCols)),
  intersect(ForcedContinuousCols, AnalysisCols)
)
ContinuousDependents <- ContinuousCols
#### Categorical columns: binary columns, non-numeric keys, and explicit overrides ####
CategoricalCols <- unique(c(
  CategoricalDependents,
  setdiff(names(dt)[!numeric_cols], c(NonAnalysisCols, ForcedContinuousCols)),
  ForcedCategoricalCols
))
CategoricalKeys <- CategoricalCols

##############################
#### Check data normality ####
##############################
#### Shapiro-Wilk normality test ####
(DTNorm <- normalityfindR(DT = dt, cols = ContinuousCols))
DTNorm$variable <- gsub(";", " ", DTNorm$variable)
openxlsx::write.xlsx(list(DataNormality = DTNorm), file = file.path(OutputDir, "NormalDistributionCheck.xlsx"), colWidths = "auto")
#### q-q plots ####
qq_plots <- lapply(ContinuousCols, function(column) ggqqplot(dt[[column]]) + ylab(column))
QQPlotPath <- ensure_dir(file.path(OutputDir, "QQPlots"))
for(i in seq_along(ContinuousCols)) {
  ggsave(file.path(QQPlotPath, paste0(gsub(";", "_", ContinuousCols[i]), ".png")),
         plot = qq_plots[[i]], width = 5, height = 5, dpi = 100)
}
#### Density plots ####
DensityPlotPath <- ensure_dir(file.path(OutputDir, "DensityPlots"))
DensityContinuousCols <- intersect(c("Age", "IndexSOFAscoreattimeofsepsis", "Lactateattimeofsepsis", "COMORBIDITYINDEX"), ContinuousCols)
DensityCategoricalDependents <- intersect(c("Migrant", "InhospitalMortality", "ICUadmission", "Mechanicalventilation"), CategoricalDependents)
DensityPlotCompile(PlotDT = dt, ContinuousCols = DensityContinuousCols,
                   CategoricalDependents = DensityCategoricalDependents,
                   SavePlots = TRUE, OutPath = DensityPlotPath)

########################################################
#### Compare continuous variables between two groups ####
########################################################
ContinuousComparisonPath <- ensure_dir(file.path(OutputDir, "ContinuousComparisonResults"))
CompiledContinuous <- continuous_compare_Compile(dt = dt, ContinuousCols = ContinuousCols, CatTable = CatTable,
                                                 Fpath = ContinuousComparisonPath)
CompiledContinuous <- CompiledContinuous[order(CompiledContinuous$mann_whitney_u_pvalue),]
openxlsx::write.xlsx(list(ContinuousAnalysis = CompiledContinuous), file = file.path(OutputDir, "CompiledContinuous.xlsx"), colWidths = "auto")
#### Visualize results with an overview bubble plot ####
BubblePlot <- ComparisonOverviewBubblePlot(Fpath = NULL,
                                           XColumn = "comparison",
                                           YColumn = "parameter",
                                           SignificantDigits = 4,
                                           RemoveAllEmpty = TRUE,
                                           SizeVariable = "mann_whitney_u_pvalue",
                                           log2SizeVariable = TRUE,
                                           SizeSubsetLevel = "<0.05",
                                           ColorVariable = "CliffsDelta",
                                           log2ColorVariable = FALSE,
                                           ColorSubsetLevel = NA,
                                           SourceData = dt,
                                           ShowMarginalBars = TRUE,
                                           MarginalAxes = "y",
                                           MarginalPanelRatio = 12,
                                           ColorLimitQuantile = 0.95,
                                           ClusterMethod = c("correlation"),
                                           ClusterAxes  = c("both"),
                                           SigColumn = NULL,
                                           SigThreshold = 0.05,
                                           Xlab = NULL,
                                           Ylab = NULL,
                                           NAFill = "grey85",
                                           ResultsDT = CompiledContinuous)
ggsave(file.path(OutputDir, "ContinuousComparisonPlotOverview.png"),
       plot = BubblePlot[["Plot"]], width = 25, height = 25, units = "in", dpi = 100,
       limitsize = FALSE)
openxlsx::write.xlsx(list(OverviewTable = BubblePlot[["Table"]]), file = file.path(OutputDir, "ContinuousComparisonPlotOverview.xlsx"), colWidths = "auto")

##################################################################################################
#### Compare categorical variables between each binary dependent and all other categoricals   ####
##################################################################################################
CategoricalComparisonPath <- ensure_dir(file.path(OutputDir, "CategoricalComparisonResultsFisher"))
CompiledCategorical <- Categorical_Compare_compile(dt = dt, CategoricalCols = CategoricalKeys, CatTable = CatTable, Type = "fisher",
                                                   Fpath = CategoricalComparisonPath)
CompiledCategorical <- CompiledCategorical[order(CompiledCategorical$fisher_exact_p_value),]
CategoricalSummary <- CompiledCategorical[comparison_type == "summary"]
openxlsx::write.xlsx(list(CategoricalAnalysis = CompiledCategorical,
                          CategoricalSummary = CategoricalSummary),
                     file = file.path(OutputDir, "CompiledCategorical.xlsx"), colWidths = "auto")
#### Visualize results with an overview bubble plot ####
BubblePlot <- ComparisonOverviewBubblePlot(Fpath = NULL,
                                           XColumn = "comparison",
                                           YColumn = "parameter",
                                           SignificantDigits = 4,
                                           RemoveAllEmpty = TRUE,
                                           SizeVariable = "fisher_exact_p_value",
                                           log2SizeVariable = TRUE,
                                           SizeSubsetLevel = "<0.05",
                                           ColorVariable = "odds_ratio",
                                           log2ColorVariable = FALSE,
                                           ColorSubsetLevel = NA,
                                           SourceData = dt,
                                           ShowMarginalBars = TRUE,
                                           MarginalAxes = "y",
                                           MarginalPanelRatio = 12,
                                           ColorLimitQuantile = 0.95,
                                           ClusterMethod = c("correlation"),
                                           ClusterAxes  = c("both"),
                                           SigColumn = NULL,
                                           SigThreshold = 0.05,
                                           Xlab = NULL,
                                           Ylab = NULL,
                                           NAFill = "grey85",
                                           ResultsDT = CategoricalSummary)
ggsave(file.path(OutputDir, "CategoricalComparisonPlotOverview.png"),
       plot = BubblePlot[["Plot"]], width = 30, height = 20, units = "in", dpi = 100,
       limitsize = FALSE)
openxlsx::write.xlsx(list(OverviewTable = BubblePlot[["Table"]]), file = file.path(OutputDir, "CategoricalComparisonPlotOverview.xlsx"), colWidths = "auto")

#######################################################################
#### Perform linear regression analyses on each dependent variable ####
#######################################################################
LinearRegressionPath <- ensure_dir(file.path(OutputDir, "LinearRegressionResults"))
CompiledLinear <- RegressionCompile(DT = dt,
                                    dependents = ContinuousDependents,
                                    independents = ContinuousCols,
                                    Fpath = LinearRegressionPath,
                                    removeAdditionalDependents = FALSE,
                                    ReturnSignificant = FALSE,
                                    Parametric = TRUE,
                                    type = "Linear")
linear_p_cols <- grep("Pr(>", colnames(CompiledLinear), fixed = TRUE, value = TRUE)
if (nrow(CompiledLinear) > 0 && length(linear_p_cols) > 0) {
  CompiledLinear <- CompiledLinear[order(CompiledLinear[[linear_p_cols[1]]]),]
} else {
  warning("CompiledLinear is empty or has no model p-value column; skipping p-value sort.")
}
openxlsx::write.xlsx(list(RegressionAnalysis = CompiledLinear), file = file.path(OutputDir, "CompiledLinear.xlsx"), colWidths = "auto")
#### Visualize results with an overview bubble plot ####
BubblePlot <- ComparisonOverviewBubblePlot(Fpath = NULL,
                                           XColumn = "dependent",
                                           YColumn = "independent",
                                           SignificantDigits = 4,
                                           RemoveAllEmpty = TRUE,
                                           SizeVariable = "Pr(>|t|)",
                                           log2SizeVariable = TRUE,
                                           SizeSubsetLevel = "<0.05",
                                           ColorVariable = "Estimate",
                                           log2ColorVariable = FALSE,
                                           ColorSubsetLevel = NA,
                                           SourceData = dt,
                                           ShowMarginalBars = TRUE,
                                           MarginalAxes = "y",
                                           MarginalPanelRatio = 12,
                                           ColorLimitQuantile = 0.95,
                                           ClusterMethod = "correlation",
                                           ClusterAxes = "both",
                                           SigColumn = NULL,
                                           SigThreshold = 0.05,
                                           Xlab = NULL,
                                           Ylab = NULL,
                                           NAFill = "grey85",
                                           ResultsDT = CompiledLinear)
ggsave(file.path(OutputDir, "LinearRegressionPlotOverview.png"),
       plot = BubblePlot[["Plot"]], width = 30, height = 20, units = "in", dpi = 100,
       limitsize = FALSE)
openxlsx::write.xlsx(list(OverviewTable = BubblePlot[["Table"]]), file = file.path(OutputDir, "LinearRegressionPlotOverview.xlsx"), colWidths = "auto")

###########################################################################
#### Perform logistic regression analyses on each categorical variable ####
###########################################################################
LogisticRegressionPath <- ensure_dir(file.path(OutputDir, "LogisticRegressionResults"))
CompiledLogistic <- RegressionCompile(DT = dt,
                                      dependents = CategoricalDependents,
                                      independents = ContinuousCols,
                                      Fpath = LogisticRegressionPath,
                                      removeAdditionalDependents = TRUE,
                                      ReturnSignificant = FALSE,
                                      Parametric = TRUE,
                                      #### Reference group for each binary dependent comes straight out of the auto-built CatTable ####
                                      CatTable = CatTable,
                                      CatReferenceCol = "Table2Subset",
                                      type = "Logistic")
logistic_p_cols <- grep("Pr(>", colnames(CompiledLogistic), fixed = TRUE, value = TRUE)
if (nrow(CompiledLogistic) > 0 && length(logistic_p_cols) > 0) {
  CompiledLogistic <- CompiledLogistic[order(CompiledLogistic[[logistic_p_cols[1]]]),]
} else {
  warning("CompiledLogistic is empty or has no model p-value column; skipping p-value sort.")
}
openxlsx::write.xlsx(list(RegressionAnalysis = CompiledLogistic), file = file.path(OutputDir, "CompiledLogistic.xlsx"), colWidths = "auto")
#### Visualize results with an overview bubble plot ####
BubblePlot <- ComparisonOverviewBubblePlot(Fpath = NULL,
                                           XColumn = "dependent",
                                           YColumn = "independent",
                                           SignificantDigits = 4,
                                           RemoveAllEmpty = TRUE,
                                           SizeVariable = "Pr(>|z|)",
                                           log2SizeVariable = TRUE,
                                           SizeSubsetLevel = "<0.05",
                                           ColorVariable = "Estimate",
                                           log2ColorVariable = FALSE,
                                           ColorSubsetLevel = NA,
                                           SourceData = dt,
                                           ShowMarginalBars = TRUE,
                                           MarginalAxes = "y",
                                           ClusterMethod = "correlation",
                                           ClusterAxes = "both",
                                           SigColumn = "p_fdr",
                                           SigThreshold = 0.05,
                                           ResultsDT = CompiledLogistic)
ggsave(file.path(OutputDir, "LogisticRegressionPlotOverview.png"),
       plot = BubblePlot[["Plot"]], width = 20, height = 20, units = "in", dpi = 100,
       limitsize = FALSE)
openxlsx::write.xlsx(list(OverviewTable = BubblePlot[["Table"]]), file = file.path(OutputDir, "LogisticRegressionPlotOverview.xlsx"), colWidths = "auto")

################################################
#### Heatmap alternative to the bubble plot ####
################################################
HeatmapResults <- SignificanceHeatmapR(ResultsDT = CompiledLogistic, RowVar = "dependent", ColVar = "independent",
                                       ValueVar = "Pr(>|z|)", Transform = "neglog10")
HeatmapResults[["Plot"]]
ggsave(file.path(OutputDir, "LogisticSignificanceHeatmap.png"), plot = HeatmapResults[["Plot"]], width = 10, height = 8, dpi = 150)

####################################################
#### Create summary table of continuous columns ####
####################################################
vtable::st(dt)
vtable::st(dt, file = file.path(OutputDir, "SummaryTable.html"))


################################################################################
################################################################################
################################################################################
#### LASSO REGRESSION ##########################################################
################################################################################
################################################################################
################################################################################

##########################
#### Define variables ####
##########################
library(svyVarSel); library(data.table); library(glmnet)
outcome_vars <- c("Death_01", "LOS", "complic_respiratory", "complic_cardiac", "complic_renal")   # the LASSO is fitted on the first; the Matching section uses all
exposure_vars <- c("sud_any_01", "elix_alcohol", "elix_drug")
pred_vars   <- c("sud_any_01", "white_01", "hispanic_01", "black_01", "other_race_01", "asian_pi_01",
                 "native_american_01", "medicare_01", "private_01", "medicaid_01", "other_payor_01",
                 "first_quartile_01", "second_quartile_01", "third_quartile_01", "fourth_quartile_01",
                 "PRDAY_first_egs_proc_positive", "sex_01", "AGE", paste0("elix_", c("chf", "carit",
                                                                                     "valv", "pcd", "pvd", "hypunc", "hypc", "para", "ond", "cpd", "diabunc", "diabc",
                                                                                     "hypothy", "rf", "ld", "pud", "aids", "lymph", "metacanc", "solidtum", "rheumd",
                                                                                     "coag", "obes", "wloss", "fed", "blane", "dane", "alcohol", "drug", "psycho", "depre")))

data(egs_misc_1621, package = "dataMineR")   #### Synthetic stand-in bundled with the package (see ?egs_misc_1621)
model_dt <- as.data.table(egs_misc_1621)[, c(outcome_vars, pred_vars, "NIS_STRATUM", "HOSP_NIS", "DISCWT", "YEAR"), with = FALSE]

#######################
#### Data cleaning ####
#######################
# Inline numeric coercion & complete case filtering
num_cols <- c(outcome_vars, pred_vars, "DISCWT")
model_dt[, (num_cols) := lapply(.SD, function(x) suppressWarnings(as.numeric(as.character(x)))), .SDcols = num_cols]
# Complete cases on the LASSO columns only: a missing extra outcome does not drop the record (each outcome model drops its own)
model_dt <- model_dt[complete.cases(model_dt[, c(outcome_vars[1], pred_vars, "NIS_STRATUM", "HOSP_NIS", "DISCWT", "YEAR"), with = FALSE])]
# Drop zero-variance predictors
(zero_var  <- pred_vars[sapply(pred_vars, function(v) uniqueN(model_dt[[v]]) < 2)])
if(length(zero_var) > 0){ pred_vars <- setdiff(pred_vars, zero_var) }

################################
#### Check Reference values ####
################################
dummy_groups <- list(race = c("white_01", "hispanic_01", "black_01", "other_race_01", "asian_pi_01", "native_american_01"),
                     payer = c("medicare_01", "private_01", "medicaid_01", "other_payor_01"),
                     income = c("first_quartile_01", "second_quartile_01", "third_quartile_01", "fourth_quartile_01"))
for (g in names(dummy_groups)){ cat(g, "N variables per row:\n"); print(table(rowSums(model_dt[, dummy_groups[[g]], with = FALSE]))) }
#### Each dummy group must has one "1" per row, so remove reference level variables that can be derived from the rest ####
pred_vars <- setdiff(pred_vars, c(race = "white_01", payer = "private_01", income = "fourth_quartile_01") )
#### Pool unique multi-year info ####
model_dt[, strata_id := paste(YEAR, NIS_STRATUM, sep = "_")]
model_dt[, cluster_id := paste(YEAR, HOSP_NIS, sep = "_")]
#### Hold out 20% of hospitals within strata for validation, before any fitting, tuning or exploration ####
#### lasso_ready_df (training) and validation_df (held out) come from the same prepared data; no hospital is in both. ####
holdout <- ClusterStratifiedHoldoutSplit(model_dt, cluster = "cluster_id", strata = "strata_id", share = 0.2, outcome = outcome_vars[1], seed = 1621)
#### Convert to data frame ####
lasso_ready_df <- as.data.frame(holdout$training)
validation_df <- as.data.frame(holdout$validation)

###################
#### Fit model ####
###################
set.seed(1621)
lasso_Fit <- svyVarSel_wlasso_adjusted(data = lasso_ready_df,
                                       col.y = outcome_vars[1],
                                       col.x = pred_vars,
                                       family = "binomial",
                                       strata = "strata_id",
                                       weights = "DISCWT",
                                       cluster = "cluster_id",
                                       method = "dCV", k = 5, R = 1)

########################
#### Review results ####
########################
#### lambda_rule: "min" = lowest CV error; "1se" = largest lambda within 1 SE of the minimum (sparser model) ####
#### Prints position, lambda, N vars and CV error for both rules, and flags lambda.min at the end of the grid ####
lambda_select <- svyVarSel_lambda_select(fit = lasso_Fit, lambda_rule = "min")
grid <- lambda_select$grid; j <- lambda_select$j; j_1se <- lambda_select$j_1se   # used by the plot and checks below

#### Plot results ####
wlasso.plot(lasso_Fit); abline(v = log(grid[j_1se]), lty = 3, col = "blue")
#### Coefficients at the chosen lambda: $coefs (all), $selected (nonzero), $excluded (names of the zero ones) ####
lasso_coefs <- svyVarSel_selected_coefs(fit = lasso_Fit, j = j)
(selected_vars <- lasso_coefs$selected)
lasso_coefs[["excluded"]]

##########################################################
#### Test the reason excluded predictors are left out ####
##########################################################
excluded_reasons <- svyVarSel_excluded_reasons(fit = lasso_Fit, data = lasso_ready_df, outcome = outcome_vars[1], weights = "DISCWT", j = j, PrintSummary = FALSE)
excluded_reasons$summary
excluded_reasons$details

#########################################################################################
#### Build optimized model by adjusting lambda settings so that a minimum is reached ####
#########################################################################################
lambda_refine <- svyVarSel_refine_lambda_grid(fit = lasso_Fit)
lasso_optimized <- svyVarSel_wlasso_tune(data = lasso_ready_df,
                                         col.y = outcome_vars[1],
                                         col.x = pred_vars,
                                         family = "binomial",
                                         strata = "strata_id",
                                         weights = "DISCWT",
                                         cluster = "cluster_id",
                                         method = "dCV",
                                         k = 5,
                                         R = 1,
                                         lambda.grid = lambda_refine$grid,
                                         max_rounds = 3,
                                         seed = 1621)
# Use minimum CV error for the prediction.
wlasso.plot(lasso_optimized)
svyVarSel_lambda_select(fit = lasso_optimized, lambda_rule = "min")

optimized_coefs <- as.matrix(lasso_optimized$model$min$beta)
optimized_coefs[optimized_coefs[, 1] != 0, , drop = FALSE]

##########################
#### Optimize K and R ####
##########################
kr_stability <- svyVarSel_kr_stability(data = lasso_ready_df,
                                       col.y = outcome_vars[1],
                                       col.x = pred_vars,
                                       family = "binomial",
                                       strata = "strata_id",
                                       weights = "DISCWT",
                                       cluster = "cluster_id",
                                       method = "dCV",
                                       K = c(5, 10),
                                       R = c(1, 3, 5),
                                       seeds = c(1621, 1622),
                                       lambda.grid = lasso_optimized$lambda$grid)
kr_stability$suggested             # K and R to use; refit lasso_optimized with them if they differ from the values above
kr_stability$selection_frequency   # share of seeds selecting each predictor, per K x R

#########################################
#### coefficient stability bootstrap ####
#########################################
coef_bootstrap <- svyVarSel_coef_bootstrap(fit = lasso_optimized, data = lasso_ready_df, outcome = outcome_vars[1], weights = "DISCWT",
                                           cluster = "cluster_id", strata = "strata_id", B = 200, seed = 1621)
coefficient_stability_table <- coef_bootstrap$table

####################################################################
#### Unpenalized association analysis with Confidence intervals ####
####################################################################
library(survey)
confidence_refit <- svyVarSel_unpenalized_refit(fit = lasso_optimized, data = lasso_ready_df, outcome = outcome_vars[1],
                                                weights = "DISCWT", cluster = "cluster_id", strata = "strata_id", PrintSummary = FALSE)
names(confidence_refit)
confidence_refit[["comparison"]]
confidence_refit[["table"]]
confidence_table <- confidence_refit$table

###################################
#### Evaluate prediction model ####
###################################
prediction_evaluation <- evaluate_prediction_model(
  model = lasso_optimized, validation_data = validation_df,
  training_clusters = unique(lasso_ready_df$cluster_id), outcome = outcome_vars[1],
  weights = "DISCWT", cluster = "cluster_id", strata = "strata_id", B = 500
)
prediction_evaluation$summary       # AUC, Brier, log loss, calibration O/E + CIs
prediction_evaluation$calibration   # Weighted observed vs predicted mortality


#####################################################################################
#### LASSO for every dependent: the outcomes and the exposure (double selection) #### #### Runs the steps above (tune, lambda, excluded reasons, bootstrap, unpenalized refit, validation) for each dependent. ####
##################################################################################### #### A LASSO on the exposure shows what predicts it; with the outcome LASSOs that gives both halves of double selection. ####
LassoDependents <- unique(c(outcome_vars, exposure_vars))   # outcomes, then every exposure used below
LassoPath <- file.path(WritePath, DirName, "LassoResults")
LassoResults <- LassoCompile(DT = lasso_ready_df,
                             dependents = LassoDependents,
                             predictors = pred_vars,             # a dependent is never its own predictor
                             family = c(LOS = "gaussian"),       # match the Linear LOS outcome model below; Poisson outcomes must be explicit
                             weights = "DISCWT", cluster = "cluster_id", strata = "strata_id",
                             validation = validation_df,         # held-out hospitals; binomial dependents only for now
                             method = "svyglm",                  # unpenalized refit: "glm", "svyglm" or "lm"; lm requires gaussian outcomes
                             cv_method = "dCV", k = 5, R = 1,    # replicate-weight cross-validation
                             lambda_rule = "min",                # "min" = lowest CV error; "1se" = sparser model
                             max_rounds = 3,                     # lambda-grid refinements until lambda.min is inside the grid
                             seed = 1621,
                             excluded_reasons = TRUE,            # why each left-out predictor was excluded
                             B = 200,                            # coefficient stability bootstrap refits (0 skips)
                             B_validation = 500,                 # bootstrap draws for the validation intervals
                             Fpath = LassoPath,                  # one LASSO_Sub_DataPredictions_<dependent>.xlsx per dependent + LASSO_README.xlsx
                             checkpoint = TRUE,                  # save each dependent as it finishes; a rerun resumes where it stopped
                             singleton = "fixed",                # "fail" stops the bootstrap step for singleton strata
                             exclude_predictors = NULL)          # optional: leave predictors out of one dependent's LASSO after a first run,
# e.g. list(sud_any_01 = c("elix_drug", "elix_alcohol")) once the overview
# plot flags them as near-duplicates of the exposure
saveRDS(LassoResults, file.path(LassoPath, "LassoResults.rds"))      # keep the result: LassoCompile() clears its checkpoints once it finishes
# LassoResults <- readRDS(file.path(LassoPath, "LassoResults.rds"))  # in a later session, run this instead of LassoCompile() (data and settings above must be unchanged)
LassoResults$models                                        # one row per dependent: lambda, stability, refit comparison, validation
LassoResults$table[LassoResults$table$selected, ]            # selected predictors of every dependent
wlasso.plot(LassoResults$fits[["LOS"]]$fit)                  # CV error curve for one dependent

####################
#### LASSO plot ####
####################
LassoPlotPath <- file.path(LassoPath, "Plots")
lasso_paths <- list(); lasso_correlations <- list()
for (dep in names(LassoResults$fits)) {
  result <- LassoResults$fits[[dep]]
  if (is.null(result$fit) || is.null(result$selection)) {
    message("Skipping LASSO plots for ", dep, ": the model did not fit. See LassoResults$models$error.")
    next
  }
  training_rows <- lasso_ready_df[result$rows, , drop = FALSE]
  lasso_paths[[dep]] <- LassoCoefficientPathPlot(result$fit, training_rows,
                                                        weights = "DISCWT", j = result$selection$j, Fpath = LassoPlotPath, PrintPlot = FALSE,
                                                        FileName = paste0(dep, "_CoefficientPaths.pdf"))
  lasso_correlations[[dep]] <- LassoCorrelationPlot(result$fit, training_rows,
                                                          weights = "DISCWT", j = result$selection$j, Fpath = LassoPlotPath, PrintPlot = FALSE,
                                                          FileName = paste0(dep, "_PredictorCorrelations.pdf"))
}

##############################
#### Bubble overview plot ####
##############################
LassoOverview <- LassoOverviewPlot(results = LassoResults,
                                   data = lasso_ready_df,             # the DT given to LassoCompile (weighted SDs)
                                   exposures = exposure_vars,         # every exposure used by the matching section
                                   weights = "DISCWT",
                                   SigColumn = "p_fdr",               # used only when ShowRefitMarks = TRUE
                                   SigThreshold = 0.05,
                                   ShowRefitMarks = FALSE,             # TRUE: outline significant and X unreliable unpenalized refits
                                   WeakThreshold = 0.1,               # confounder strong if |standardized exposure-LASSO coefficient| >= this
                                   phi_threshold = 0.8,               # near-duplicate of the exposure: |r| >= this (MatchingAssumptionCheck)
                                   joint_min_cor = 0.3,               # codes that together nearly are the exposure (MatchingAssumptionCheck)
                                   share_threshold = 0.01,            # extreme exposed share (MatchingAssumptionCheck)
                                   ReferenceLine = 0.7,               # pane 1 dashed line at |r| (display only; NULL = no line)
                                   SizeBreaks = c(25, 50, 75, 100),   # bootstrap % shown in the size legend
                                   ShowNotSelected = FALSE,           # TRUE adds predictors selected for no dependent
                                   ColorLimitQuantile = 0.95,         # colour scale capped at this quantile of |coefficient|
                                   ClusterMethod = "correlation",     # row order within each role: "correlation", "signed" or "none"
                                   Fpath = LassoPlotPath,
                                   FileName = "LASSO_SelectionOverview.pdf",
                                   ShowRoles = FALSE,
                                   ShowRefitMarks = FALSE,
                                   PrintPlot = TRUE)
LassoOverview$roles   # each predictor's role and the dependents it was selected for

#########################################
#### Recommended models for matching ####
#########################################
RecommendedModels <- do.call(rbind, Filter(Negate(is.null), lapply(exposure_vars, function(exposure_name) {
  tryCatch({
    recommended <- MatchingSetupFromLasso(results = LassoResults,
                                          roles = LassoOverview$roles_by_exposure,
                                          exposure = exposure_name,
                                          DT = model_dt,
                                          dependents = outcome_vars,
                                          method = "svyglm",                      # "glm", "svyglm" or "lm"
                                          type = c(LOS = "Linear"),               # aligned with the Gaussian LOS LASSO above
                                          include_weak = TRUE,                    # FALSE: match on strong confounders only
                                          outcome_adjust = TRUE,                  # FALSE: no outcome-model covariates
                                          exclude = "PRDAY_first_egs_proc_positive")   # measured after the exposure (possible mediator)
    recommended$ModelID <- paste0("LASSO_", exposure_name, "_", seq_len(nrow(recommended)))
    recommended
  }, error = function(e) {
    message("No LASSO-recommended matching model for ", exposure_name, ": ", conditionMessage(e))
    NULL
  })
})))
RecommendedModels

################################################################################
################################################################################
################################################################################
#### Matching ##################################################################
################################################################################
################################################################################
################################################################################

########################
#### Matching Setup ####
########################
library(MatchIt); library(cobalt); library(survey); library(sandwich)
exposure_var <- exposure_vars
CandidateCovariates <- setdiff(pred_vars, exposure_var)
BalanceCovariates <- unique(c(CandidateCovariates, unlist(dummy_groups), "YEAR"))
BalanceCovariates <- intersect(BalanceCovariates, names(model_dt))
MatchingSetupPath <- file.path(WritePath, DirName)   # flag and setup files; results go in subfolders below
# dir.create(MatchingSetupPath, recursive = TRUE, showWarnings = FALSE)

#### Flag overlap and timing concerns; preserve every candidate ####
Positivity <- MatchingAssumptionCheck(model_dt, exposure = exposure_var, covariates = CandidateCovariates,
                                      weights = "DISCWT", post_exposure = NULL,  # Supply column names here when timing is known.
                                      estimand = c("ATT", "ATE"),                # ATT for matching / EB, ATE for IPTW
                                      IPTWSampleSizeMinimum = 0.10)
openxlsx::write.xlsx(list(README = Positivity$readme, Flags = Positivity$table, Checks = Positivity$checks), file.path(MatchingSetupPath, "Matching_CovariateFlags.xlsx"))

#################################################
#### Set up matching regression combinations ####
#################################################
RegressionSetup <- MatchingRegressionSetup(DT = model_dt,
                                           dependents = outcome_vars,
                                           independents = exposure_var,
                                           matching = unique(Positivity$table$variable),
                                           single = TRUE,
                                           joint = TRUE,
                                           method = "svyglm", # "glm", "svyglm" or "lm"
                                           type = c(LOS = "Linear"),
                                           flags = Positivity$table,          # flags annotate; never remove
                                           examples = TRUE,                   # adds example rows (see Notes) at the end of the same sheet; delete any you do not want to run
                                           recommended = RecommendedModels)   # LASSO double-selection rows (LASSO_1, ...), placed first
openxlsx::write.xlsx(list(Setup = RegressionSetup), file.path(MatchingSetupPath, "Matching_RegressionSetup.xlsx"), colWidths = "auto")   # edit this file, then read it back below; rerunning this line overwrites your edits

#########################################
#### Propensity score matching setup ####
#########################################
MatchingPath <- file.path(WritePath, DirName, "PropensityScore")
dir.create(MatchingPath, recursive = TRUE, showWarnings = FALSE)
RegressionSetup <- openxlsx::read.xlsx(file.path(MatchingSetupPath, "Matching_RegressionSetup.xlsx"), sheet = "Setup")   # rerun from here after editing

#############################
#### perform PS matching ####
#############################
PSResults <- PSRegression(DT = model_dt, setup = RegressionSetup,
                          weights = "DISCWT", cluster = "cluster_id", strata = "strata_id",
                          Fpath = MatchingPath, balance_covariates = BalanceCovariates,
                          flags = Positivity$table,
                          checkpoint = TRUE,           # save each matching as it finishes; a rerun resumes where it stopped
                          detail = FALSE)              # TRUE adds every per-covariate balance column and writes the full Balance + OutcomeBalance sheets
PSResults
#### Visualize PS matching results with an overview bubble plot ####
PSBubblePlot <- MatchedComparisonOverviewBubblePlot(Fpath = MatchingPath,
                                                    XColumn = "dependent",
                                                    YColumn = "independent",
                                                    #### Added to the row label only when they differ for the same exposure x outcome
                                                    LabelColumns = c("MatchingVariable", "Method", "Type", "estimand", "trim", "moments",
                                                                     "IndependentVariable", "Ratio", "PSCaliperSD", "CovariateCalipers"),
                                                    LabelSeparator = " | ",
                                                    SignificantDigits = 4,
                                                    RemoveAllEmpty = TRUE,
                                                    SizeVariable = "Pr(>|t|)",   # the p-value column for every Type
                                                    log2SizeVariable = TRUE,
                                                    SizeSubsetLevel = "<0.05",
                                                    ColorVariable = "Estimate",  # log odds / log ratio, or difference for Linear
                                                    log2ColorVariable = FALSE,
                                                    ColorSubsetLevel = NA,
                                                    NColumn = "NcountsNoNA",     # marginal bars: patients in each matched analysis
                                                    ShowMarginalBars = TRUE,
                                                    MarginalAxes = "y",
                                                    MarginalPanelRatio = 8,
                                                    ColorLimitQuantile = 0.95,
                                                    ClusterMethod = "correlation",
                                                    ClusterAxes = "both",
                                                    #### Outline bubbles surviving per-dependent FDR
                                                    SigColumn = "p_fdr",
                                                    SigThreshold = 0.05,
                                                    Xlab = NULL,
                                                    Ylab = NULL,
                                                    NAFill = "grey85"); PSBubblePlot[["Plot"]]
png(file = file.path(WritePath, DirName, "PSPlotOverview.png"), width = 2000, height = 2000, units = "px", res = 100); print(PSBubblePlot[["Plot"]]); dev.off()
openxlsx::write.xlsx(list(OverviewTable = PSBubblePlot[["Table"]], Labels = PSBubblePlot[["Labels"]]), file = file.path(WritePath, DirName, "PSPlotOverview.xlsx"), colWidths = "auto")


#######################################################
#### Inverse Probability Treatment Weighting setup ####
#######################################################
IPTWPath <- file.path(WritePath, DirName, "IPTW")
dir.create(IPTWPath, recursive = TRUE, showWarnings = FALSE)
IPTWSetup <- openxlsx::read.xlsx(file.path(MatchingSetupPath, "Matching_RegressionSetup.xlsx"), sheet = "Setup")

#############################################
#### perform IPTW followed by regression ####
#############################################
IPTWResults <- IPTWRegression(DT = model_dt, setup = IPTWSetup,
                              weights = "DISCWT", cluster = "cluster_id", strata = "strata_id",
                              estimand = "ATE",            # "ATE" (everyone), "ATT" (exposed) or "ATO" (overlap population)
                              trim = NULL,                 # NULL, or a weight quantile to cap at, e.g. 0.99
                              Fpath = IPTWPath, balance_covariates = BalanceCovariates,
                              flags = Positivity$table,
                              checkpoint = TRUE,           # save each weighting as it finishes; a rerun resumes where it stopped
                              detail = FALSE)              # TRUE adds every per-covariate balance column and writes the full Balance + OutcomeBalance sheets
IPTWResults
#### Visualize IPTW results with an overview bubble plot ####
IPTWBubblePlot <- MatchedComparisonOverviewBubblePlot(Fpath = IPTWPath,
                                                      XColumn = "dependent",
                                                      YColumn = "independent",
                                                      LabelColumns = c("MatchingVariable", "Method", "Type", "estimand", "trim", "moments",
                                                                       "IndependentVariable", "Ratio", "PSCaliperSD", "CovariateCalipers"),
                                                      LabelSeparator = " | ",
                                                      SignificantDigits = 4,
                                                      RemoveAllEmpty = TRUE,
                                                      SizeVariable = "Pr(>|t|)",
                                                      log2SizeVariable = TRUE,
                                                      SizeSubsetLevel = "<0.05",
                                                      ColorVariable = "Estimate",
                                                      log2ColorVariable = FALSE,
                                                      ColorSubsetLevel = NA,
                                                      NColumn = "NcountsNoNA",     # marginal bars: patients in each weighted analysis
                                                      ShowMarginalBars = TRUE,
                                                      MarginalAxes = "y",
                                                      MarginalPanelRatio = 8,
                                                      ColorLimitQuantile = 0.95,
                                                      ClusterMethod = "correlation",
                                                      ClusterAxes = "both",
                                                      SigColumn = "p_fdr",
                                                      SigThreshold = 0.05,
                                                      Xlab = NULL,
                                                      Ylab = NULL,
                                                      NAFill = "grey85"); IPTWBubblePlot[["Plot"]]
png(file = file.path(WritePath, DirName, "IPTWPlotOverview.png"), width = 2000, height = 2000, units = "px", res = 100); print(IPTWBubblePlot[["Plot"]]); dev.off()
openxlsx::write.xlsx(list(OverviewTable = IPTWBubblePlot[["Table"]], Labels = IPTWBubblePlot[["Labels"]]), file = file.path(WritePath, DirName, "IPTWPlotOverview.xlsx"), colWidths = "auto")


###########################
#### Entropy Balancing ####
###########################
EBPath <- file.path(WritePath, DirName, "EntropyBalancing")
dir.create(EBPath, recursive = TRUE, showWarnings = FALSE)
EBSetup <- openxlsx::read.xlsx(file.path(MatchingSetupPath, "Matching_RegressionSetup.xlsx"), sheet = "Setup")

#########################################################
#### perform entropy balancing followed by regression ####
#########################################################
EBResults <- EBRegression(DT = model_dt, setup = EBSetup,
                          weights = "DISCWT", cluster = "cluster_id", strata = "strata_id",
                          estimand = "ATT",            # "ATT" (controls reweighted to the exposed) or "ATE" (everyone)
                          moments = 1,                 # 1 balances means; 2 also balances variances
                          Fpath = EBPath, balance_covariates = BalanceCovariates,
                          flags = Positivity$table,
                          checkpoint = TRUE,           # save each balancing as it finishes; a rerun resumes where it stopped
                          detail = FALSE)              # TRUE adds every per-covariate balance column and writes the full Balance + OutcomeBalance sheets
EBResults
#### Visualize entropy balancing results with an overview bubble plot ####
EBBubblePlot <- MatchedComparisonOverviewBubblePlot(Fpath = EBPath,
                                                    XColumn = "dependent",
                                                    YColumn = "independent",
                                                    LabelColumns = c("MatchingVariable", "Method", "Type", "estimand", "trim", "moments",
                                                                     "IndependentVariable", "Ratio", "PSCaliperSD", "CovariateCalipers"),
                                                    LabelSeparator = " | ",
                                                    SignificantDigits = 4,
                                                    RemoveAllEmpty = TRUE,
                                                    SizeVariable = "Pr(>|t|)",
                                                    log2SizeVariable = TRUE,
                                                    SizeSubsetLevel = "<0.05",
                                                    ColorVariable = "Estimate",
                                                    log2ColorVariable = FALSE,
                                                    ColorSubsetLevel = NA,
                                                    NColumn = "NcountsNoNA",     # marginal bars: patients in each balanced analysis
                                                    ShowMarginalBars = TRUE,
                                                    MarginalAxes = "y",
                                                    MarginalPanelRatio = 8,
                                                    ColorLimitQuantile = 0.95,
                                                    ClusterMethod = "correlation",
                                                    ClusterAxes = "both",
                                                    SigColumn = "p_fdr",
                                                    SigThreshold = 0.05,
                                                    Xlab = NULL,
                                                    Ylab = NULL,
                                                    NAFill = "grey85"); EBBubblePlot[["Plot"]]
png(file = file.path(WritePath, DirName, "EBPlotOverview.png"), width = 2000, height = 2000, units = "px", res = 100); print(EBBubblePlot[["Plot"]]); dev.off()
openxlsx::write.xlsx(list(OverviewTable = EBBubblePlot[["Table"]], Labels = EBBubblePlot[["Labels"]]), file = file.path(WritePath, DirName, "EBPlotOverview.xlsx"), colWidths = "auto")



################################################################################
################################################################################
#### Misc ######################################################################
################################################################################
################################################################################
MiscPath <- ensure_dir(file.path(OutputDir, "Misc"))

#####################################################################################
#### All-pairs Spearman correlation, with scatterplots of the significant pairs  ####                                           ####
#####################################################################################
res2 <- Hmisc::rcorr(as.matrix(dt[, ContinuousCols, with = FALSE]), type = "spearman")
flat <- flattenCorrMatrix(res2$r, res2$P, res2$n)
flat <- flat[order(abs(flat$cor), decreasing = TRUE),]
openxlsx::write.xlsx(list(AllCorrelations = flat), file = file.path(MiscPath, "AllCorrelations.xlsx"), colWidths = "auto")
CorrelationPlotPath <- ensure_dir(file.path(MiscPath, "CorrelationScatterplots"))
flatSig <- flat[flat$p < 0.05,]
for(i in seq_len(min(nrow(flatSig), 10))){
  p <- scatterplotgraph(dataset = dt, xaxis = flatSig$row[i], yaxis = flatSig$column[i],
                        xlabel = flatSig$row[i], ylabel = flatSig$column[i],
                        title = paste("correlation between", flatSig$row[i], "and", flatSig$column[i]))
  ggsave(file.path(CorrelationPlotPath,
                   paste0(gsub("[^A-Za-z0-9_]", "_", paste(flatSig$row[i], flatSig$column[i], sep = "_")), ".png")),
        plot = p, width = 6, height = 5, dpi = 150)
}

##########################################################################################################
#### Multi-variable regression on significant single-variable predictors, plus confounding screening. ####
#### which significant predictors lose significance when co-included with another predictor? #############
##########################################################################################################
MultiRegressionPath <- ensure_dir(file.path(MiscPath, "MultiVariableRegression"))
CompiledSignificantLinear <- RegressionCompile(DT = dt, dependents = ContinuousDependents, independents = ContinuousCols,
                                               Fpath = MultiRegressionPath,
                                               ReturnSignificant = TRUE, removeAdditionalDependents = FALSE,
                                               Parametric = TRUE, type = "Linear")
MultiLinearModels <- MultiModelRegression(DataTable = dt, deps = ContinuousDependents,
                                          IndependentPredictorDT = CompiledSignificantLinear,
                                          UseSignificantPredictorsOnly = TRUE, parametric = TRUE,
                                          CalculateAUC = FALSE, Type = "Linear", ReturnSignificantOnly = "F")
openxlsx::write.xlsx(list(MultiLinearModels = MultiLinearModels), file = file.path(MiscPath, "MultiLinearModels.xlsx"), colWidths = "auto")
PairedVarResults <- PairedMultiVariableSummary(DT = dt, CompiledSignificant = CompiledSignificantLinear, AddBack = NA)
openxlsx::write.xlsx(list(PairedRegression = PairedVarResults[["PairedRegressionDT"]],
                          InhibitionTable  = PairedVarResults[["InhibitionTable"]],
                          SignificantPredsInhibition = PairedVarResults[["SignificantPredsInhibitionDT"]]),
                     file = file.path(MiscPath, "PredictorConfounding.xlsx"), colWidths = "auto")

#####################################
#### Simple two-group comparison ####
#####################################
(fisherMigrant <- categorical_compare(DT = dt, compair = "Migrant-Native", compair_column = GroupBy,
                                       keys = setdiff(CategoricalKeys, GroupBy), Stattest = "fisher"))
openxlsx::write.xlsx(list(CategoricalMigrantNative = fisherMigrant), file = file.path(MiscPath, "CategoricalMigrantNative.xlsx"), colWidths = "auto")
DTOne <- dt[dt$Migrant == 1,]; DTTwo <- dt[dt$Migrant == 0,]
(continuousMigrant <- continuous_compare(DT1 = DTOne, DT2 = DTTwo, compair = "Migrant-Native", key = ContinuousCols))
openxlsx::write.xlsx(list(ContinuousMigrantNative = continuousMigrant), file = file.path(MiscPath, "ContinuousMigrantNative.xlsx"), colWidths = "auto")

##################################################
#### Survival analysis: Kaplan-Meier + Cox PH ####
##################################################
SurvivalCovariates <- setdiff(c(CategoricalDependents, ContinuousCols), c("InhospitalMortality", "mortality30day", "mortality90day"))
SurvivalResults <- SurvivalCompileR(dt = dt, AdmitCol = "DateofAdmission", DeathCol = "DateofDeath",
                                    DischargeCol = "DateofDischarge", StratifyBy = "Migrant",
                                    DeathFormat = "%m/%d/%Y",
                                    Covariates = SurvivalCovariates, TimeUnits = "days", Fpath = MiscPath)
SurvivalResults[["KMOverallPlot"]]
SurvivalResults[["KMStratifiedPlot"]]
SurvivalResults[["LogRankP"]]
SurvivalResults[["UnivariateCox"]]
ggsave(file.path(MiscPath, "KaplanMeier_Overall.png"), plot = SurvivalResults[["KMOverallPlot"]], width = 6, height = 5, dpi = 150)
ggsave(file.path(MiscPath, "KaplanMeier_ByMigrant.png"), plot = SurvivalResults[["KMStratifiedPlot"]], width = 6, height = 5, dpi = 150)

############################
#### Mediation analysis ####
############################
MediationTable <- data.table(DependentVariable = "ICUadmission", IndependentVariable = "Migrant", PutativeMediator = "hpi2score")
MediationResults <- MediationTestR(DT = copy(dt), MediationTable = MediationTable)
openxlsx::write.xlsx(list(IndependentModels = MediationResults[["IndependentModels"]],
                          MediationAnalysis  = MediationResults[["MediationAnalysis"]]),
                     file = file.path(MiscPath, "MediationAnalysis.xlsx"), colWidths = "auto")
#### Multi-mediator PLS-SEM: do several socioeconomic factors jointly mediate Migrant -> ICU admission? ####
MediationVariables <- c("Migrant", "economic", "education", "housing", "ICUadmission")
dtMediation <- dt[stats::complete.cases(dt[, MediationVariables, with = FALSE])]
PLSMediation <- MultiMediationPLSR(dt = dtMediation, Treatment = "Migrant",
                                   Mediators = c("economic", "education", "housing"),
                                   Outcome = "ICUadmission", Nboot = 100, Seed = WorkflowSeed)
PLSMediation[["IndirectEffects"]]

########################################################
#### Predictor "stability" via bootstrap resampling ####
########################################################
RegressionStabilityPath <- ensure_dir(file.path(MiscPath, "RegressionStability"))
StabilityDependents  <- intersect(c("Lactateattimeofsepsis", "IndexSOFAscoreattimeofsepsis"), ContinuousDependents)
StabilityIndependents <- intersect(c("Age", "COMORBIDITYINDEX", "hpi2score", "economic"), ContinuousCols)
StabilityResults <- RegressionStabilityR(dt = dt, dependents = StabilityDependents,
                                         independents = StabilityIndependents,
                                         Nboots = 15, Fraction = 0.5, Seed = WorkflowSeed)
openxlsx::write.xlsx(list(RegressionStability = StabilityResults), file = file.path(RegressionStabilityPath, "RegressionStability.xlsx"), colWidths = "auto")

########################################
#### ANCOVA / interaction screening ####
########################################
AncovaResults <- AncovaInteractionR(dt = dt, Dependent = "IndexSOFAscoreattimeofsepsis",
                                    Covariates = c("Age", "Migrant", "Homeless"),
                                    InteractionPairs = list(c("Age", "Migrant")))
AncovaResults[["MainEffects"]]
AncovaResults[["InteractionEffects"]]

######################################
#### Effect size / power analysis ####
######################################
PowerAnalysisR(test = "t", x = dt[Migrant == 1,]$Age, y = dt[Migrant == 0,]$Age, power = 0.80)
PowerAnalysisR(test = "summary_means", mu1 = 3.4, mu2 = 2.8, sd1 = 3.5, sd2 = 3.1, power = 0.80)

#########################################################################################
#### ROC threshold sweep, DeLong ROC comparison, and Bland-Altman agreement analysis ####
#########################################################################################
data(dream_data, peru_data, copd_data, package = "dataMineR")
dtDream <- as.data.table(dream_data)
ThresholdSweep <- ThresholdSweepAUCR(dt = dtDream, Outcome = "opioids", Predictor = "ahi", Direction = "ge")
ThresholdSweep[["Plot"]]
DeLongResults <- ROCDeLongTest(DataList = list(DreamData = dtDream), DepVariable = "opioids",
                               IndVariable1 = "ahi", IndVariable2 = "rdi", paired = TRUE, Direction = "<")
DeLongResults[["DreamData"]][["ROCDT"]]

dtPeru <- as.data.table(peru_data)
BAResults <- BlandAltman(DataList = list(PeruData = dtPeru), Var1 = "hcto_1", Var2 = "hcto_2",
                         SigLevel = 0.95, ColorFeature = FALSE, ColorBy = NA, ScaleData = FALSE)
BAResults[["BARaw"]][["PeruData"]]
BAResults[["BAPlot"]][["PeruData"]]

#####################################################################
#### Dimensionality reduction, pairs plots, and grouped boxplots ####
#####################################################################
dtCopd <- as.data.table(copd_data)
CopdFeatures <- c("FEV1", "FVC", "sao2", "AHI", "CRP", "hemoglobin", "leukocytes", "cholesterol")
PCAResults <- PCAPlotR(dt = dtCopd, features = CopdFeatures, GroupCol = "group")
PCAResults[["Plot"]]
UMAPResults <- UMAPPlotR(dt = dtCopd, features = CopdFeatures, GroupCol = "group", seed = WorkflowSeed)
UMAPResults[["Plot"]]
PairsPlotResults <- PairsPlotR(dt = dtCopd, features = CopdFeatures, GroupCol = "group")
PairsPlotResults
BoxplotPath <- ensure_dir(file.path(MiscPath, "Boxplots"))
BoxplotResults <- GroupedBoxplotStatR(dt = dtCopd, ValueCols = CopdFeatures, GroupCol = "group",
                                      OutPath = BoxplotPath, SavePlots = TRUE)
