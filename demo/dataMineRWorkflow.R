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
  "Hmisc", "vtable", "seminr"
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
#### Misc #######################################################################
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
