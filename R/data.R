#' Synthetic COPD/OSA overlap-syndrome biomarker dataset
#'
#' Synthetic patient-level data for demonstrating group-comparison,
#' dimensionality-reduction, and boxplot workflows across four cohorts:
#' COPD, OSA (obstructive sleep apnea), their overlap ("OVERLAP"), and
#' healthy controls.
#'
#' @format A \code{data.table} (see \pkg{data.table}) with 180 rows and 20 variables, including
#'   \code{group} (\code{"OVERLAP"}, \code{"COPD"}, \code{"CONTROL"},
#'   \code{"OSA_BASAL"}), demographics (\code{age}, \code{sex}), spirometry
#'   (\code{FEV1}, \code{FVC}, \code{FEV1_FVC}), sleep/oxygenation
#'   (\code{AHI}, \code{sao2}), labs (\code{CRP}, \code{hemoglobin},
#'   \code{leukocytes}, \code{cholesterol}), and patient-reported/functional
#'   measures (\code{CAT}, \code{distance_6min}, \code{active_smoker}).
#' @source Randomly generated; see \code{data-raw/create_example_data.R}.
#' @examples
#' data(copd_data)
#' PCAPlotR(copd_data, features = c("FEV1", "FVC", "sao2", "AHI"),
#'          GroupCol = "group")
"copd_data"

#' Synthetic migrant-status sepsis cohort dataset
#'
#' Synthetic patient-level data modeled on a hospital sepsis registry,
#' comparing migrant vs. non-migrant patients across demographics,
#' comorbidities, sepsis severity, geographic/social-deprivation indices,
#' and outcomes (mortality, ICU admission, length of stay).
#'
#' @format A \code{data.table} (see \pkg{data.table}) with 155 rows and 83 variables, including the
#'   binary grouping variable \code{Migrant}; admission/discharge/death
#'   date-time columns (\code{DateofAdmission}, \code{DateofDischarge},
#'   \code{DateofDeath}); sepsis severity (\code{IndexSOFAscoreattimeofsepsis},
#'   \code{Lactateattimeofsepsis}); ~20 binary comorbidity flags (\code{MI},
#'   \code{CHF}, \code{Diabetes}, ...); outcomes
#'   (\code{InhospitalMortality}, \code{mortality30day}, \code{mortality90day},
#'   \code{ICUadmission}); and neighborhood social-deprivation indices
#'   (\code{hpi2score}, \code{economic}, \code{education}, \code{housing}, ...).
#' @source Randomly generated; see \code{data-raw/create_example_data.R}.
#' @examples
#' data(dfMigrant)
#' InitialTablebg(dfMigrant, GroupBy = "Migrant", p_values = TRUE)
"dfMigrant"

#' Synthetic cardiac arrest resuscitation dataset
#'
#' Synthetic patient-level data modeled on a cardiac-arrest resuscitation
#' registry, with CPR-quality metrics, arrest etiology, and outcomes.
#' Several columns include realistic missingness patterns.
#'
#' @format A \code{data.table} (see \pkg{data.table}) with 3455 rows and 16 variables, including
#'   end-tidal CO2 (\code{ETCO2_init}, \code{ETCO2_time}), CPR-quality
#'   metrics (\code{CPR_time}, \code{CCF}, \code{Depth}, \code{Rate},
#'   \code{defib_attempts}), arrest characteristics (\code{Etiology},
#'   \code{OHCA}, \code{IHCA_subsequent}), and outcomes (\code{ROSC},
#'   \code{Survival to DC}, \code{CPC score}).
#' @source Randomly generated; see \code{data-raw/create_example_data.R}.
"cardiac_data"

#' Synthetic opioid-use polysomnography (PSG) dataset
#'
#' Synthetic sleep-study data comparing patients with and without opioid
#' use, including standard PSG metrics (respiratory events, sleep staging,
#' oxygenation) and PAP-therapy adherence data.
#'
#' @format A \code{data.table} (see \pkg{data.table}) with 676 rows and 117 variables, including the
#'   binary grouping variable \code{opioids}; respiratory-event indices
#'   (\code{rdi}, \code{ahi}, \code{cai_2}); sleep staging
#'   (\code{n1}, \code{n2}, \code{n3}, \code{rem}); oxygenation
#'   (\code{spo2_mean}, \code{spo2_nadir}); daytime sleepiness (\code{eds},
#'   \code{ess_ini}); and PAP-therapy adherence metrics
#'   (\code{days_fu_cpap}, \code{percent_days_4h}, \code{pressure_median_resmed}, ...).
#' @source Randomly generated; see \code{data-raw/create_example_data.R}.
#' @examples
#' data(dream_data)
#' ThresholdSweepAUCR(dream_data, Outcome = "opioids", Predictor = "ahi",
#'                     Direction = "ge")
"dream_data"

#' Synthetic Peru high-altitude cognition dataset
#'
#' Synthetic data modeled on a high-altitude acclimatization study,
#' combining demographics, pulmonary function, baseline/overnight sleep
#' studies, and a battery of cognitive-performance tasks (psychomotor
#' vigilance, Trail Making, Attention Network Test, balloon analog risk
#' task). Many task-derived columns are missing-at-random per participant
#' (task not completed).
#'
#' @format A \code{data.table} (see \pkg{data.table}) with 48 rows and 89 variables, including
#'   demographics/altitude (\code{age}, \code{altitude}, \code{cms_score});
#'   duplicate hematocrit measurements (\code{hcto_1}, \code{hcto_2}) useful
#'   for a Bland-Altman agreement demo; pulmonary function
#'   (\code{fvc_actual}, \code{fev1_actual}); sleep-study metrics
#'   (\code{ahi_b}, \code{nadir_spo2}, \code{totalahi}); and cognitive-task
#'   outputs (\code{pvt_medianrt}, \code{tmta_sec}, \code{accuracy}, ...).
#' @source Randomly generated; see \code{data-raw/create_example_data.R}.
#' @examples
#' data(peru_data)
#' BlandAltman(DataList = list(Peru = peru_data), Var1 = "hcto_1",
#'             Var2 = "hcto_2", SigLevel = 0.95, ColorFeature = FALSE,
#'             ColorBy = NA, ScaleData = FALSE)
"peru_data"

#' SmartPill gastric-motility dataset with derived group/demographic labels
#'
#' A lightly reshaped copy of \code{medicaldata::smartpill}, with
#' human-readable \code{group} (Trauma/Healthy), \code{sex}, and \code{race}
#' factor labels derived from their coded originals, a derived \code{bmi}
#' column, and \code{.} replaced with \code{_} in column names.
#'
#' @format A \code{data.table} (see \pkg{data.table}) with 95 rows and 26 variables. See
#'   \code{\link[medicaldata]{smartpill}} for the source column
#'   definitions; \code{group}, \code{sex}, \code{race}, and \code{bmi} are
#'   added here.
#' @source \pkg{medicaldata}; see \code{data-raw/create_example_data.R}.
"MData"

#' NHANES adults subset with a derived diabetes/smoking group
#'
#' One row per adult (age >= 18) participant from the \pkg{NHANES} package's
#' 2009-2012 survey cycle, with a derived four-level \code{group} combining
#' diabetes and current-smoking status.
#'
#' @format A \code{data.table} (see \pkg{data.table}) with 4381 rows and 77 variables. See
#'   \code{\link[NHANES]{NHANES}} for the source column definitions;
#'   \code{group} (one of \code{"DIABETIC_SMOKER"}, \code{"DIABETIC_NONSMOKER"},
#'   \code{"HEALTHY"}, \code{"SMOKER_ONLY"}) is added here.
#' @source \pkg{NHANES}; see \code{data-raw/create_example_data.R}.
"nhanes_like"
