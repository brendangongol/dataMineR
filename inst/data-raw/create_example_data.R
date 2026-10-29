################################################################################
#### inst/data-raw/create_example_data.R #######################################
#### Builds the package's example datasets and saves them to                ####
#### inst/extdata/<name>.rda (one object per file, named like the file).    ####
#### Re-run this script from the package root (e.g. with devtools::load_all ####
#### or source()) whenever the generation logic below changes; it is not    ####
#### run automatically at install/build time. After installing, load a      ####
#### dataset with                                                           ####
####   load(system.file("extdata", "dfMigrant.rda", package = "dataMineR")) ####
################################################################################
################################################################################
#### Create Synthetic COPD/OSA biomarker dataset for package demonstration #####
################################################################################
set.seed(42)
N <- 180
groups <- rep(c("OVERLAP", "COPD", "CONTROL", "OSA_BASAL"), each = 45)
is_patient <- groups %in% c("OVERLAP", "COPD")
copd_data <- data.frame(id = sample(10000:999999, N),group = groups,
                        sex = ifelse(groups %in% c("OVERLAP","COPD"),
                                     sample(c("Male","Female"), N, replace=TRUE, prob=c(0.87,0.13)),
                                     sample(c("Male","Female"), N, replace=TRUE, prob=c(0.55,0.45))),
                        age = round(ifelse(is_patient, rnorm(N, 69, 7), rnorm(N, 44, 9))),
                        FEV1 = round(ifelse(is_patient, pmax(0.3, rnorm(N, 1.7, 0.6)), pmax(2.0, rnorm(N, 3.3, 0.8))), 2),
                        FVC = round(ifelse(is_patient, pmax(0.5, rnorm(N, 3.1, 0.8)), pmax(2.5, rnorm(N, 4.0, 0.9))), 2),
                        sao2 = round(ifelse(is_patient, pmax(70, rnorm(N, 93, 3)), pmax(90, rnorm(N, 97, 1)))),
                        AHI = round(pmax(0, ifelse(groups == "OVERLAP", rnorm(N, 38, 20),
                                                   ifelse(groups == "OSA_BASAL", rnorm(N, 40, 18),
                                                          ifelse(groups == "COPD", rnorm(N,  2,  3), rnorm(N, 4, 3)))))),
                        CRP = round(pmax(0.01, rgamma(N, shape = ifelse(is_patient, 1.5, 0.8), rate = 0.9)), 2),
                        hemoglobin = round(pmax(8, rnorm(N, 15.0, 1.8)), 1),
                        leukocytes = round(pmax(2, rnorm(N, 7.5, 2.2)), 1),
                        cholesterol = round(pmax(80, rnorm(N, 200,  45))),
                        CAT = ifelse(is_patient, round(pmax(0, pmin(40, rnorm(N, 13, 7)))), NA),
                        distance_6min = ifelse(is_patient, round(pmax(60, rnorm(N, 385, 100))), NA),
                        active_smoker = ifelse(is_patient, rbinom(N, 1, 0.28), NA),
                        PLATE = sample(1:7, N, replace = TRUE),
                        stringsAsFactors = FALSE)
copd_data$FEV1_FVC <- round(copd_data$FEV1 / copd_data$FVC, 3)
copd_data$FEV1_FVC[is_patient] <- pmin(copd_data$FEV1_FVC[is_patient], 0.69)
copd_data$osa <- sample(0:1, nrow(copd_data), replace = TRUE)
copd_data$copd <- sample(0:1, nrow(copd_data), replace = TRUE)
copd_data$ovs <- sample(0:1, nrow(copd_data), replace = TRUE)

####################################################################
#### Create Synthetic Migrant dataset for package demonstration ####
####################################################################
set.seed(42)
N <- 155
has_geo <- rbinom(N, 1, 0.40)
make_z <- function(mn = 0, sd = 0.5) ifelse(has_geo, round(rnorm(N, mn, sd), 6), NA)
make_p <- function() ifelse(has_geo, round(pmax(0, pmin(100, rnorm(N, 50, 28))), 6), NA)
admit_date <- as.POSIXct("2019-01-01") + runif(N, 0, as.numeric(as.POSIXct("2023-06-01") - as.POSIXct("2019-01-01")))
discharge_date <- admit_date + round(pmax(1, rnorm(N, 9, 7))) * 86400
sepsis_time    <- admit_date + runif(N, -3600, 14400)
county_fips_v  <- ifelse(has_geo, sample(c(6073, 6025), N, replace = TRUE, prob = c(0.92, 0.08)), NA)
hpi2_pctile_v  <- ifelse(has_geo, round(pmax(0, pmin(100, rnorm(N, 42, 28))), 6), NA)
dfMigrant <- data.frame(
  MRN = sample(10000000:29999999, N),
  CSN = sample(60000000000:69999999999, N),
  DateofDeath = ifelse(rbinom(N, 1, 0.065) == 1, format(discharge_date + round(runif(N, 0, 30)) * 86400, "%m/%d/%Y"), NA),
  Migrant  = rbinom(N, 1, 0.49),
  Homeless = rbinom(N, 1, 0.35),
  MigrationDetails = paste0(
    sample(c("SW note - undocumented",
             "SW note - Restricted Medical r/t Immigration status",
             "SW note - undocumented immigrant",
             "OSH note - Restricted M-cal due to immigration status"),
           N, replace = TRUE), " - ", format(as.Date("2005-01-01") + sample(0:6000, N, replace = TRUE), "%m/%d/%Y")),
  ModeofArrival = sample(
    c("Paramedic Unit", "Automobile", "Walk in", "Basic/BLS Ambulance",
      "Critical Care Transport", "Helicopter", "Other", "Taxi"), N, replace = TRUE, prob = c(0.38, 0.28, 0.14, 0.10, 0.05, 0.02, 0.02, 0.01)),
  ICDCode = sample(
    c("A41.9", "A41.59", "A41.9,R65.21", "A41.51", "A41.89,A41.51",
      "A41.81", "A41.9,R65.20", "A41.01"), N, replace = TRUE, prob = c(0.35, 0.20, 0.18, 0.12, 0.06, 0.04, 0.03, 0.02)),
  DateofAdmission      = format(admit_date,      "%m/%d/%Y %H:%M"),
  DateofDischarge      = format(discharge_date,  "%m/%d/%Y %H:%M"),
  DischargeDisposition = sample(
    c("Home Routine", "Home Healthcare IP Admit Related", "Expired",
      "Long Term Care Facility", "Acute Care Facility (not this hospital)",
      "Hospice Home", "Hospice in Place", "Against Medical Advice"),
    N, replace = TRUE, prob = c(0.38, 0.32, 0.065, 0.08, 0.06, 0.025, 0.025, 0.025)),
  Age = round(pmax(18, pmin(90, rnorm(N, 47, 13))), 6),
  Sex = sample(c("Male", "Female"), N, replace = TRUE, prob = c(0.42, 0.58)),
  Race = sample(
    c("Other Race or Mixed Race", "White", "Black or African American",
      "Asian", "Undetermined", "Unknown (Patient cannot or refuses to declare race)"),
    N, replace = TRUE, prob = c(0.58, 0.18, 0.10, 0.06, 0.05, 0.03)),
  Ethnicity = sample(
    c("Other Hispanic, Latino(a) or Spanish origin", "Hispanic/Latino",
      "Not Hispanic, Latino/a, or Spanish origin",
      "Unreported (Patient refused or chose not to disclose)"),
    N, replace = TRUE, prob = c(0.60, 0.22, 0.12, 0.06)),
  Language = sample(
    c("Spanish", "English", "Arabic", "Other", "Portuguese", "Russian"),
    N, replace = TRUE, prob = c(0.60, 0.26, 0.05, 0.05, 0.02, 0.02)),
  InsuranceStatus = sample(
    c("Medicaid - California", "Medicaid Managed Care",
      "County Medical Services", "Medicaid - Out of State"),
    N, replace = TRUE, prob = c(0.65, 0.20, 0.10, 0.05)),
  TimeofSepsis = format(sepsis_time, "%m/%d/%Y %H:%M"),
  IndexSOFAscoreattimeofsepsis = round(pmax(0, pmin(16, rnorm(N, 3.4, 3.5)))),
  Lactateattimeofsepsis = round(pmax(0.4, rnorm(N, 2.95, 3.3)), 1),
  Administrationofvasoactivemedications = rbinom(N, 1, 0.24),
  Mechanicalventilation = rbinom(N, 1, 0.18),
  ICUadmission = rbinom(N, 1, 0.29),
  Suspectedsourceofinfection = sample( c("Unknown", "GU", "Lung", "GI", "Bacteremia", "SSTI", "Bone", "GI/GU", "Lung/GU", "Bone/Bacteremia"),
    N, replace = TRUE, prob = c(0.20, 0.18, 0.13, 0.12, 0.10, 0.08, 0.07, 0.05, 0.04, 0.03)),
  ComfortcareHospiceatdischarge = rbinom(N, 1, 0.08),
  MI = rbinom(N, 1, 0.11),
  CHF = rbinom(N, 1, 0.31),
  Peripheralvasculardisease = rbinom(N, 1, 0.12),
  CVAorTIA = rbinom(N, 1, 0.10),
  Dementia = rbinom(N, 1, 0.05),
  Chronicpulmonarydisease = rbinom(N, 1, 0.26),
  Connectivetissuedisease = rbinom(N, 1, 0.04),
  PUD = rbinom(N, 1, 0.14),
  Liverdisease = rbinom(N, 1, 0.49),
  Diabetes = rbinom(N, 1, 0.50),
  HemiplegiaParaplegia = rbinom(N, 1, 0.05),
  CKD = rbinom(N, 1, 0.31),
  Solidtumor = rbinom(N, 1, 0.60),
  Leukemia = rbinom(N, 1, 0.05),
  Lymphoma = rbinom(N, 1, 0.03),
  AIDS = rbinom(N, 1, 0.47),
  COMORBIDITYINDEX = round(pmax(0, pmin(12, rnorm(N, 3.6, 3.1)))),
  InhospitalMortality = rbinom(N, 1, 0.065),
  mortality30day = rbinom(N, 1, 0.02),
  mortality90day = rbinom(N, 1, 0.007),
  SERVICE = sample(
    c("Medicine Hillcrest", "General Medicine", "Pulmonary / Critical Care (HC)",
      "Emergency General Surgery (HC)", "Surgery Colorectal (LJ)",
      "Bone Marrow Transplant", "Family Medicine", "ED Medicine"),
    N, replace = TRUE, prob = c(0.35, 0.20, 0.15, 0.10, 0.07, 0.06, 0.04, 0.03)),
  CensusTract = ifelse(has_geo, round(runif(N, 6.07e9, 6.08e9)), NA),
  pop2010 = ifelse(has_geo, round(rnorm(N, 5134, 1583)), NA),
  City = NA_character_,
  ZIP = ifelse(has_geo, sample(c(91911, 91950, 92103, 92113, 92116, 92123, 92283, 92102), N, replace = TRUE), NA),
  County_FIPS = county_fips_v,
  County_Name = ifelse(county_fips_v == 6073, "San Diego", ifelse(county_fips_v == 6025, "Imperial", NA)),
  UrbanType = ifelse(has_geo, sample(c("urban_area", "urban_cluster", "rural"), N, replace = TRUE, prob = c(0.75, 0.15, 0.10)), NA),
  hpi2score = ifelse(has_geo, round(rnorm(N, -0.15, 0.60), 6), NA),
  hpi2_pctile_pos = hpi2_pctile_v,
  hpi_least_healthy_25pct = ifelse(has_geo, ifelse(hpi2_pctile_v <= 25, "Yes", "No"), NA),
  quintiles = ifelse(has_geo, sample(1:5, N, replace = TRUE), NA),
  quartiles = ifelse(has_geo, sample(1:4, N, replace = TRUE), NA),
  economic = make_z(), economic_pctile = make_p(),
  education = make_z(), education_pctile = make_p(),
  housing = make_z(), houseing_pctile = make_p(),
  healthcareaccess = make_z(), healthcareaccess_pctile = make_p(),
  neighborhood = make_z(), neighborhood_pctile = make_p(),
  pollution = make_z(), pollution_pctile = make_p(),
  transportation = make_z(), transportation_pctile = make_p(),
  social = make_z(), social_pctile = make_p(),
  insured = ifelse(has_geo, round(pmax(50, pmin(100, rnorm(N, 83, 12))), 1), NA),
  insured_pctile = make_p(),
  LEB = ifelse(has_geo, round(pmax(200, rnorm(N, 250, 60)), 2), NA),
  LEB_pctile = make_p(),
  white_pct = ifelse(has_geo, round(pmax(0, pmin(100, rnorm(N, 33.9, 25.0))), 6), NA),
  black_pct = ifelse(has_geo, round(pmax(0, pmin(100, rnorm(N,  8.0,  7.0))), 6), NA),
  asian_pct = ifelse(has_geo, round(pmax(0, pmin(100, rnorm(N,  8.5,  8.0))), 6), NA),
  latino_pct = ifelse(has_geo, round(pmax(0, pmin(100, rnorm(N, 46.8, 28.7))), 6), NA),
  multiple_pct = ifelse(has_geo, round(pmax(0, pmin(10,  rnorm(N,  3.0,  2.0))), 6), NA),
  stringsAsFactors = FALSE )

##########################################
#### Synthetic cardiac arrest dataset ####
##########################################
set.seed(42)
N <- 3455
mask_na <- function(x, n_present) { idx <- sample(N, N - n_present); x[idx] <- NA; x }
cardiac_data <- data.frame(
  ETCO2_init = mask_na(round(pmax(0, pmin(98, rnorm(N, 20.1, 15.8))), 1), 758),
  ETCO2_time = mask_na(round(pmax(0, pmin(36, rexp(N, rate = 0.32))), 1), 689),
  defib_attempts = mask_na(round(pmax(0, pmin(39, rexp(N, rate = 0.69)))), 1114),
  CPR_time = mask_na(round(pmax(0, pmin(180, rexp(N, rate = 0.068)))), 1015),
  CCF = mask_na(round(pmax(0.39, pmin(0.99, rnorm(N, 0.93, 0.04))), 2), 929),
  Depth = mask_na( round(pmax(0.73, pmin(5.0, rnorm(N, 2.84, 0.45))), 2), 931),
  Rate = mask_na(round(pmax(3.5, pmin(149, rnorm(N, 114.3, 9.5))), 1), 932),
  Age = mask_na(round(pmax(0, pmin(145, rnorm(N, 59.7, 18.0)))), 3435),
  Sex = mask_na(sample(c(1, 0), N, replace = TRUE, prob = c(0.665, 0.335)), 3438),
  Etiology = mask_na( sample(c("PEA", "VF", "VT", "Asystole", "asystole", "PEA/asystole", "Vagal", "VF/VT"),
           N, replace = TRUE, prob = c(0.745, 0.111, 0.087, 0.019, 0.018, 0.004, 0.002, 0.001)), 3421),
  Matrix = mask_na(sample(c("D1", "R2", "C1", "U", "D2", "C4", "C2", "R1", "R4", "C3", "N1", "Unknown"), N, replace = TRUE,
           prob = c(0.198, 0.146, 0.131, 0.130, 0.080, 0.079, 0.060, 0.039, 0.036, 0.026, 0.011, 0.003)), 3422),
  ROSC = mask_na(rbinom(N, 1, prob = 0.723), 3428),
  `Survival to DC` = mask_na(rbinom(N, 1, prob = 0.295), 3373),
  `CPC score` = mask_na(sample(c("5", "1", "2", "3", "4"), N, replace = TRUE, prob = c(0.695, 0.156, 0.070, 0.055, 0.011)), 3363),
  OHCA = ifelse(rbinom(N, 1, prob = 941 / N) == 1, 1, NA),
  IHCA_subsequent = ifelse(rbinom(N, 1, prob = 433 / N) == 1, 1, NA),
  stringsAsFactors = FALSE, check.names      = FALSE )

##################################################
#### Synthetic Opioid PSG/sleep study dataset ####
##################################################
set.seed(42)
N <- 676
is_opio <- c(rep(TRUE, 169), rep(FALSE, 507))  # opioid cohort flag
rdi_v <- round(pmax(5,  rnorm(N, ifelse(is_opio, 40, 33), 28)), 1)
ahi_v <- round(pmax(0,  rnorm(N, ifelse(is_opio, 30, 25), 28)), 1)
spo2_v <- round(pmax(81, rnorm(N, 93.6, 2.4)), 1)
nadir_v <- round(pmax(35, rnorm(N, 80.8, 9.3)))
tst_v <- round(pmax(10, rnorm(N, 262, 109)), 1)
n1_v <- round(pmax(0,  rnorm(N, 13.4, 13.2)), 1)
n2_v <- round(pmax(7.9,pmin(97.1, rnorm(N, 61.8, 15.2))), 1)
n3_v <- round(pmax(0,  rnorm(N, 12.9, 14.3)), 1)
rem_v <- pmax(0, 100 - n1_v - n2_v - n3_v)   # constrain to 100%
rem_v <- round(pmin(rem_v, 56.2), 1)
supine_v <- round(pmax(0, pmin(100, rnorm(N, 58.7, 35.6))), 2)
dream_data <- data.frame(
  study_id = sample(1:999, N),
  cohort = ifelse(is_opio, "Opioids", "No opioids"),
  redcap_repeat_instrument = NA_real_,
  redcap_repeat_instance = NA_real_,
  sex = sample(c("Female", "Male"), N, replace = TRUE, prob = c(0.57, 0.43)),
  age = round(pmax(20, pmin(89, rnorm(N, 57.3, 12.2)))),
  bmi = round(pmax(15.2, rnorm(N, 34.7, 8.4)), 1),
  date_psg = format( as.Date("2016-01-01") + sample(0:800, N, replace = TRUE), "%Y-%m-%d"),
  rdi = rdi_v,
  ahi = ahi_v,
  cai_2 = round(pmax(0, ifelse(is_opio, rnorm(N, 4.0, 10), rnorm(N, 1.8, 6))), 1),
  cai_percent_of_rdi = round(pmax(0, rnorm(N, 5.5, 12)), 1),
  no_oas = round(pmax(0, rnorm(N, 70, 55))),
  no_cas = round(pmax(0, ifelse(is_opio, rnorm(N, 6, 15), rnorm(N, 1, 4)))),
  no_mas = round(pmax(0, rnorm(N, 5, 12))),
  no_hypos = round(pmax(0, rnorm(N, 38, 35))),
  no_reras = round(pmax(0, rnorm(N, 10, 18))),
  percent_hypos = round(pmax(0, pmin(100, rnorm(N, 48, 30))), 1),
  calc_rdi_ahi = round(pmax(0, rnorm(N, 9, 12)), 1),
  rdi_nrem = round(pmax(0, rnorm(N, 35, 33)), 1),
  rdi_rem = round(pmax(0, rnorm(N, 28, 35)), 1),
  rdi_supine = round(pmax(0, rnorm(N, 57, 40)), 1),
  rdi_lateral_calc = round(pmax(0, rnorm(N, 18, 26)), 1),
  spo2_mean = spo2_v,
  spo2_nadir = nadir_v,
  calc_spo2_mean_nadir = round(spo2_v - nadir_v, 1),
  spo2_lt89 = round(pmax(0, rgamma(N, shape = 0.8, rate = 0.11)), 2),
  bicarb_1 = round(pmax(16, pmin(39, rnorm(N, 25.9, 2.9)))),
  bicarb_2 = round(pmax(16, pmin(39, rnorm(N, 26.1, 2.9)))),
  bicarb_3 = round(pmax(16, pmin(39, rnorm(N, 26.0, 2.9)))),
  date_initial_visit = format( as.Date("2015-06-01") + sample(0:600, N, replace = TRUE), "%Y-%m-%dT08:00:00Z"),
  ess_ini = round(pmax(0, pmin(24, rnorm(N, 10, 5)))),
  co2_monitoring = sample(c("No", "Yes"), N, replace = TRUE, prob = c(0.82, 0.18)),
  tst_abstracted = tst_v,
  se = round(pmax(3.4,  pmin(99, rnorm(N, 72.3, 17.6))), 1),
  sl = round(pmax(0,    rnorm(N, 27.4, 33.5)), 1),
  reml = round(pmax(0,    rnorm(N, 134, 83.3)), 1),
  awakenings = round(pmax(0,    rnorm(N, 30, 22))),
  stage_shifts = round(pmax(0,    rnorm(N, 115, 65))),
  n1 = n1_v,
  n2 = n2_v,
  n3 = n3_v,
  rem = rem_v,
  supine = supine_v,
  lateral = round(pmax(0, pmin(100, 100 - supine_v)), 2),
  mean_event_duration = round(pmax(5, rnorm(N, 20, 8)), 1),
  max_event_duration  = round(pmax(10, rnorm(N, 48, 20)), 1),
  ari_spont = round(pmax(0, rnorm(N, 9.5, 9)), 1),
  ari_plms = round(pmax(0, rnorm(N, 4.5, 8)), 1),
  plmi = round(pmax(0, rnorm(N, 14.1, 25)), 1),
  ari_total = round(pmax(2.5, rnorm(N, 28.3, 23.7)), 1),
  ari_resp = round(pmax(0,   rnorm(N, 14.1, 14.5)), 1),
  smoke_py = NA_real_,
  days_fu_cpap = ifelse(rbinom(N, 1, 0.75), round(pmax(1, rnorm(N, 80, 55))), NA),
  days_initial_psg = round(pmax(1, rnorm(N, 60, 40))),
  percent_days_4h = ifelse(rbinom(N, 1, 0.70), round(pmax(0, pmin(100, rnorm(N, 60, 28))), 1), NA),
  pressure_median_resmed = ifelse(rbinom(N, 1, 0.65), round(pmax(4, rnorm(N, 8.5, 2.5)), 1), NA),
  pressure_95_resmed = ifelse(rbinom(N, 1, 0.65), round(pmax(4, rnorm(N, 11.2, 3.1)), 1), NA),
  leak_median_resmed = ifelse(rbinom(N, 1, 0.65), round(pmax(0, rnorm(N, 5.5, 9)), 1), NA),
  leak_95_resmed = ifelse(rbinom(N, 1, 0.65), round(pmax(0, rnorm(N, 28, 22)), 1), NA),
  ai_flow_resmed = ifelse(rbinom(N, 1, 0.65), round(pmax(0, rnorm(N, 5, 8)), 1), NA),
  ahi_flow_resmed = ifelse(rbinom(N, 1, 0.65), round(pmax(0, rnorm(N, 6, 9)), 1), NA),
  cai_flow_resmed = ifelse(rbinom(N, 1, 0.65), round(pmax(0, rnorm(N, 1, 3)), 1), NA),
  oai_flow_resmed = ifelse(rbinom(N, 1, 0.65), round(pmax(0, rnorm(N, 2, 4)), 1), NA),
  loop_gain_1 = round(pmax(0.14, pmin(1.46, rnorm(N, 0.60, 0.19))), 6),
  loop_gain_n = round(pmax(0.14, pmin(0.91, rnorm(N, 0.41, 0.11))), 6),
  delay = round(pmax(5.12, pmin(23.28, rnorm(N, 10.36, 2.57))), 6),
  vra = round(pmax(0, pmin(100, rnorm(N, 37.4, 24.8))), 6),
  ar_threshold_pup = round(pmax(69.39, pmin(331.67, rnorm(N, 125.4, 32.2))), 6),
  vpassive = round(pmax(0, rnorm(N, 85, 25)), 6),
  vactive = round(pmax(0, rnorm(N, 90, 30)), 6),
  vcomp = round(rnorm(N, 5, 20), 6),
  vmin = round(pmax(0, rnorm(N, 25, 25)), 6),
  ahi_pup = round(pmax(0, rnorm(N, 27.4, 20.2)), 6),
  t90 = round(pmax(0, pmin(100, rgamma(N, shape = 0.6, rate = 0.08))), 6),
  severity = sample(c("Severe", "Mild", "Moderate"), N, replace = TRUE, prob = c(0.451, 0.300, 0.249)),
  cai_gt5 = sample(c("CAI<5/hr",  "CAI>5/hr"),  N, replace = TRUE, prob = c(0.903, 0.097)),
  cai_gt15 = sample(c("CAI<15/hr", "CAI>15/hr"), N, replace = TRUE, prob = c(0.958, 0.042)),
  obesity = sample(c("Morbidly obese", "Obese", "Overweight", "Normal"), N, replace = TRUE,
                    prob = c(0.463, 0.234, 0.188, 0.115)),
  eds = sample(c("No", "Yes"), N, replace = TRUE, prob = c(0.570, 0.430)),
  opioids_ini = ifelse(is_opio, "Yes", "No"),
  opioids_psg = ifelse(is_opio, "Yes", "No"),
  opioids_fu = ifelse(is_opio, "Yes", "No"),
  opioids_any = ifelse(is_opio, "Yes", "No"),
  opioids = as.integer(is_opio),
  COV_NREM = round(pmax(0.09, pmin(0.79, rnorm(N, 0.25, 0.10))), 6),
  COV_REM = round(pmax(0.00, pmin(0.66, rnorm(N, 0.27, 0.09))), 6),
  COV_SLEEP = round(pmax(0.09, pmin(0.70, rnorm(N, 0.25, 0.09))), 6),
  COV_WAKE = round(pmax(0.09, pmin(0.80, rnorm(N, 0.26, 0.09))), 6),
  COV_NREM_noEvts = round(pmax(0.09, pmin(0.79, rnorm(N, 0.24, 0.10))), 6),
  COV_REM_noEvts = round(pmax(0.00, pmin(0.66, rnorm(N, 0.26, 0.09))), 6),
  COV_SLEEP_noEvts = round(pmax(0.09, pmin(0.70, rnorm(N, 0.24, 0.09))), 6),
  COV_WAKE_noEvts = round(pmax(0.09, pmin(0.80, rnorm(N, 0.25, 0.09))), 6),
  COV_WAKE_first = round(pmax(0.09, pmin(0.80, rnorm(N, 0.25, 0.09))), 6),
  RR_NREM = round(pmax(1.79, pmin(7.44, rnorm(N, 3.82, 0.74))), 6),
  RR_REM = round(pmax(1.79, pmin(7.44, rnorm(N, 3.50, 0.70))), 6),
  RR_SLEEP = round(pmax(1.79, pmin(7.44, rnorm(N, 3.78, 0.73))), 6),
  RR_WAKE = round(pmax(1.79, pmin(7.44, rnorm(N, 3.20, 0.65))), 6),
  RR_NREM_noEvts = round(pmax(1.79, pmin(7.44, rnorm(N, 3.83, 0.74))), 6),
  RR_REM_noEvts = round(pmax(1.79, pmin(7.44, rnorm(N, 3.51, 0.70))), 6),
  RR_SLEEP_noEvts = round(pmax(1.79, pmin(7.44, rnorm(N, 3.79, 0.73))), 6),
  RR_WAKE_noEvts = round(pmax(1.79, pmin(7.44, rnorm(N, 3.21, 0.65))), 6),
  Ti_NREM = round(pmax(0.83, pmin(3.81, rnorm(N, 1.81, 0.41))), 6),
  Ti_REM = round(pmax(0.83, pmin(3.81, rnorm(N, 1.65, 0.38))), 6),
  Ti_SLEEP = round(pmax(0.83, pmin(3.81, rnorm(N, 1.79, 0.40))), 6),
  Ti_WAKE = round(pmax(0.83, pmin(3.81, rnorm(N, 1.55, 0.38))), 6),
  Ti_NREM_noEvts = round(pmax(0.83, pmin(3.81, rnorm(N, 1.82, 0.41))), 6),
  Ti_REM_noEvts = round(pmax(0.83, pmin(3.81, rnorm(N, 1.66, 0.38))), 6),
  Ti_SLEEP_noEvts = round(pmax(0.83, pmin(3.81, rnorm(N, 1.80, 0.40))), 6),
  Ti_WAKE_noEvts = round(pmax(0.83, pmin(3.81, rnorm(N, 1.56, 0.38))), 6),
  HFrEF_pmh = rbinom(N, 1, 0.08),
  HFpEF_pmh = rbinom(N, 1, 0.10),
  chf = sample(c("No CHF", "HFrEF", "HFpEF"), N, replace = TRUE, prob = c(0.84, 0.08, 0.08)),
  afib_pmh = rbinom(N, 1, 0.14),
  afib = sample(c("No A fib", "A fib"), N, replace = TRUE, prob = c(0.86, 0.14)),
  copd_pmh = rbinom(N, 1, 0.18),
  copd = sample(c("No COPD", "COPD"), N, replace = TRUE, prob = c(0.82, 0.18)),
  stringsAsFactors = FALSE )

########################################################
#### Synthetic Peru high-altitude cognition dataset ####
########################################################
set.seed(42)
N <- 48
peru_data <- data.frame(
  study_id = paste0(sample(c("CdP","AGS","BST","CCA","CLE","DRE","EHS"), N, replace = TRUE), "_",
                          sprintf("%03d", sample(70:200, N, replace = FALSE))),
  biological_sex = rep(c("F", "M"), each = N / 2),
  age = round(pmax(19, pmin(65, rnorm(N, 44.4, 13.2)))),
  dob = format(as.Date("1960-01-01") + sample(0:14600, N, replace = TRUE), "%m/%d/%y"),
  altitude = sample(c(4300, 4350, 4200, 4338, 3600, 3500), N, replace = TRUE,
                          prob = c(0.42, 0.21, 0.13, 0.13, 0.06, 0.05)),
  occupation = sample(c("Student", "house wife", "farmer", "shoe keeper", "teacher", "merchant"),
                          N, replace = TRUE, prob = c(0.35, 0.20, 0.18, 0.10, 0.10, 0.07)),
  treatment = sample(c("No", "no"), N, replace = TRUE, prob = c(0.91, 0.09)),
  height = round(pmax(1.40, pmin(1.80, rnorm(N, 1.54, 0.07))), 2),
  weight = round(pmax(46.5, pmin(89.0, rnorm(N, 64.4, 10.7))), 1),
  bmi = round(pmax(20.4, pmin(37.7, rnorm(N, 27.1, 4.2))), 2),
  hcto_1 = round(pmax(38, pmin(75, rnorm(N, 53.5, 7.3)))),
  hcto_2 = round(pmax(38, pmin(75, rnorm(N, 53.5, 7.3)))),
  hcto_av = round(pmax(40, pmin(72, rnorm(N, 53.5, 7.3))), 1),
  cms_score = sample(c(0,1,2,3,4,5,8,9,12,23), N, replace = TRUE,
                    prob = c(0.44, 0.21, 0.06, 0.15, 0.02, 0.02, 0.02, 0.02, 0.04, 0.02)),
  hr = round(pmax(56, pmin(97, rnorm(N, 71.1, 9.2)))),
  spo2 = round(pmax(70, pmin(90.5, rnorm(N, 84.5, 4.7))), 1),
  sbp_av = round(pmax(57.5, pmin(154, rnorm(N, 112.1, 18.7))), 1),
  dbp_av = round(pmax(50.5, pmin(97, rnorm(N, 73.1, 11.0))), 1),
  fvc_pred = round(pmax(2.0, rnorm(N, 3.1, 0.5)), 2),
  fev1_pred = round(pmax(1.5, rnorm(N, 2.7, 0.4)), 2),
  fev1fvc_pred = round(pmax(75, pmin(100, rnorm(N, 88, 5)))),
  fvc_actual = round(pmax(1.8, rnorm(N, 3.1, 0.6)), 2),
  fev1_actual = round(pmax(1.4, rnorm(N, 2.7, 0.5)), 2),
  fev1fvc_actual = round(pmax(70, pmin(100, rnorm(N, 86, 6)))),
  pef_actual = round(pmax(3, rnorm(N, 7.5, 2.2)), 2),
  sleep_efficiency_b = ifelse(rbinom(N,1,0.83), round(pmax(39,pmin(100,rnorm(N,75.5,10))),1), NA),
  n1_b = ifelse(rbinom(N,1,0.83), round(pmax(0, rnorm(N,18,14)),1), NA),
  n2_b = ifelse(rbinom(N,1,0.83), round(pmax(10,pmin(95,rnorm(N,62,16))),1), NA),
  n3_b = ifelse(rbinom(N,1,0.83), round(pmax(0, rnorm(N,10,10)),1), NA),
  rem_b = ifelse(rbinom(N,1,0.83), round(pmax(0, rnorm(N,14,9)),1), NA),
  odi_total_b = ifelse(rbinom(N,1,0.83), round(pmax(0, rnorm(N,18,15)),1), NA),
  nadir_saturation_b = ifelse(rbinom(N,1,0.67), round(pmax(62,pmin(90,rnorm(N,72,5)))), NA),
  mean_saturation_b = ifelse(rbinom(N,1,0.75), round(pmax(75,pmin(91,rnorm(N,81,3.3))),1), NA),
  ahi_b = ifelse(rbinom(N,1,0.83), round(pmax(0, rnorm(N,20.3,17.2)),1), NA),
  arousal_index_b = ifelse(rbinom(N,1,0.83), round(pmax(0, rnorm(N,19,12)),1), NA),
  totalahi = ifelse(rbinom(N,1,0.92), round(pmax(0, rnorm(N,18.4,17.2)),1), NA),
  totalari = ifelse(rbinom(N,1,0.92), round(pmax(5.8,pmin(51.9,rnorm(N,21.4,11.8))),1), NA),
  sleepefficiency = ifelse(rbinom(N,1,0.92), round(pmax(0.39,pmin(0.86,rnorm(N,0.68,0.13))),2), NA),
  duration_sleep_min = ifelse(rbinom(N,1,0.92), round(pmax(50,rnorm(N,370,90)),1), NA),
  nremahi = ifelse(rbinom(N,1,0.92), round(pmax(0, rnorm(N,16,17)),1), NA),
  remahi = ifelse(rbinom(N,1,0.92), round(pmax(0, rnorm(N,24,20)),1), NA),
  spo2mean_sleep = ifelse(rbinom(N,1,0.92), round(pmax(50,pmin(95,rnorm(N,81,4))),2), NA),
  nadir_spo2 = ifelse(rbinom(N,1,0.92), round(pmax(50,pmin(90,rnorm(N,72,5)))), NA),
  spo2below90p_sleep_prct = ifelse(rbinom(N,1,0.92), round(pmax(0,pmin(100,rnorm(N,30,30))),2), NA),
  hvr = ifelse(rbinom(N,1,0.92), round(rnorm(N, 0.25, 0.19), 6), NA),
  hcvr = ifelse(rbinom(N,1,0.92), round(pmax(0, rnorm(N, 1.38, 0.5)), 6), NA),
  endtidal_co2 = ifelse(rbinom(N,1,0.92), round(pmax(25,pmin(50,rnorm(N,36,4))),2), NA),
  sao2 = ifelse(rbinom(N,1,0.92), round(pmax(70,pmin(98,rnorm(N,84,5))),2), NA),
  ams_score = ifelse(rbinom(N,1,0.40), round(pmax(0, rnorm(N,4,3))), NA),
  pvt_resp = ifelse(rbinom(N,1,0.94), round(pmax(50,  rnorm(N,700,200))), NA),
  pvt_lapses = ifelse(rbinom(N,1,0.94), round(pmax(0,   rnorm(N,15.9,14.8))), NA),
  pvt_errors = ifelse(rbinom(N,1,0.94), round(pmax(0,   rnorm(N,2,3))), NA),
  pvt_medianrt = ifelse(rbinom(N,1,0.94), round(pmax(205, rnorm(N,384,247)),1), NA),
  pvt_meanrt = ifelse(rbinom(N,1,0.94), round(pmax(229, rnorm(N,440,270)),1), NA),
  pvt_minrt = ifelse(rbinom(N,1,0.94), round(pmax(140, rnorm(N,200,40)),1), NA),
  pvt_maxrt = ifelse(rbinom(N,1,0.94), round(pmax(500, rnorm(N,950,400)),1), NA),
  pvt_stddevrt = ifelse(rbinom(N,1,0.94), round(pmax(30,  rnorm(N,350,200)),3), NA),
  pvt_efficiency = ifelse(rbinom(N,1,0.94), round(pmax(0,   pmin(100,rnorm(N,82,15)))), NA),
  tmta_sec = ifelse(rbinom(N,1,0.98), round(pmax(22.6,pmin(150,rnorm(N,63.8,32.9))),2), NA),
  tmta_mistakes = ifelse(rbinom(N,1,0.98), round(pmax(0,   rnorm(N,1.5,2))), NA),
  tmtb_sec = ifelse(rbinom(N,1,0.98), round(pmax(47.9,pmin(650,rnorm(N,145.5,128.4))),2), NA),
  tmtb_mistakes = ifelse(rbinom(N,1,0.98), round(pmax(0,   rnorm(N,2,3))), NA),
  vpa_night = ifelse(rbinom(N,1,0.92), round(pmax(0,pmin(16,rnorm(N,5,4))),1), NA),
  vpa_morn = ifelse(rbinom(N,1,0.92), round(pmax(0,pmin(16,rnorm(N,5,4))),1), NA),
  vpa_mem = ifelse(rbinom(N,1,0.92), round(pmax(0,pmin(16,rnorm(N,4,4))),1), NA),
  ert_corr = ifelse(rbinom(N,1,0.94), round(pmax(6, pmin(28,rnorm(N,20.3,5.9)))), NA),
  erc_incorr = ifelse(rbinom(N,1,0.94), round(pmax(0, rnorm(N,3.7,2.5))), NA),
  ert_meanrt = ifelse(rbinom(N,1,0.94), round(pmax(2099,rnorm(N,5902,2903)),1), NA),
  ert_efficiency = ifelse(rbinom(N,1,0.94), round(pmax(0,pmin(1,rnorm(N,0.84,0.12))),4), NA),
  balloonscoll = ifelse(rbinom(N,1,0.92), round(pmax(0,  rnorm(N,9,5))), NA),
  balloonspop = ifelse(rbinom(N,1,0.92), round(pmax(0,  rnorm(N,7,5))), NA),
  totalpumps = ifelse(rbinom(N,1,0.92), round(pmax(50, rnorm(N,130,60))), NA),
  meanpumps = ifelse(rbinom(N,1,0.92), round(pmax(2,  rnorm(N,4.2,1.5)),2), NA),
  bartaccuracy = ifelse(rbinom(N,1,0.92), round(pmax(0,pmin(1,rnorm(N,0.56,0.25))),4), NA),
  psqi = ifelse(rbinom(N,1,0.71), round(pmax(1,pmin(12,rnorm(N,5.7,2.8)))), NA),
  grunt1 = ifelse(rbinom(N,1,0.92), round(pmax(0,pmin(10,rnorm(N,4,2)))), NA),
  grunt2 = ifelse(rbinom(N,1,0.92), round(pmax(0,pmin(10,rnorm(N,4,2)))), NA),
  grunt3 = ifelse(rbinom(N,1,0.92), round(pmax(0,pmin(10,rnorm(N,5,2)))), NA),
  grunt4 = ifelse(rbinom(N,1,0.92), round(pmax(0,pmin(10,rnorm(N,4,2)))), NA),
  alerteffect = ifelse(rbinom(N,1,0.90), round(rnorm(N,  30, 60), 1), NA),
  orientingeffect = ifelse(rbinom(N,1,0.90), round(rnorm(N,  55, 35), 1), NA),
  conflict_effect = ifelse(rbinom(N,1,0.90), round(rnorm(N, 110, 55), 1), NA),
  grand_mean_effect = ifelse(rbinom(N,1,0.90), round(rnorm(N, 680,200), 1), NA),
  accuracy = ifelse(rbinom(N,1,0.90), round(pmax(50,pmin(100,rnorm(N,90,10))),1), NA),
  nocue_congruent_median_rt = ifelse(rbinom(N,1,0.90), round(pmax(400,rnorm(N,660,140)),1), NA),
  nocue_incongruent_median_rt = ifelse(rbinom(N,1,0.90), round(pmax(400,rnorm(N,740,180)),1), NA),
  doublecue_congruent_median_rt = ifelse(rbinom(N,1,0.90), round(pmax(400,rnorm(N,630,120)),1), NA),
  spatialcue_congruent_median_rt = ifelse(rbinom(N,1,0.90), round(pmax(350,rnorm(N,590,120)),1), NA),
  spatialcue_incongruent_median_rt= ifelse(rbinom(N,1,0.90), round(pmax(400,rnorm(N,680,150)),1), NA),
  stringsAsFactors = FALSE)

###############################################################################
#### Use medical data for demonstration #######################################
###############################################################################
if (!requireNamespace("medicaldata", quietly = TRUE))
  stop("Package 'medicaldata' is required to rebuild the example data.")
MData <- medicaldata::smartpill
#### Make group labels human-readable ####
MData$group <- factor(MData$Group, levels = c(0, 1), labels = c("Trauma", "Healthy"))
MData$sex <- factor(MData$Gender, levels = c(0, 1), labels = c("Female", "Male"))
MData$race <- factor(MData$Race, levels = 1:5, labels = c("White", "Black", "Asian_PI", "Hispanic", "Other"))
#### Derive BMI after converting height from centimeters to meters ####
MData$bmi <- round(MData$Weight / (MData$Height / 100)^2, 1)
#### Select and rename to match target structure ####
colnames(MData) <- gsub("\\.", "_", colnames(MData))  # replace dots with underscores

#################################
#### NHANES clinical dataset #### 2009-2012 (10,000 rows, 76 variables)
#################################
if (!requireNamespace("NHANES", quietly = TRUE))
  stop("Package 'NHANES' is required to rebuild the example data.")
#### Adults only, one row per person, complete key variables ####
nhanes_source <- NHANES::NHANES
nhanes_like <- nhanes_source[nhanes_source$Age >= 18, ]
nhanes_like <- nhanes_like[!duplicated(nhanes_like$ID), ]
#### Derive a 4-level group from diabetes + smoking status ####
# DIABETIC_SMOKER: heaviest disease burden
# DIABETIC_NONSMOKER 
# HEALTHY: no diabetes, never smoked
# SMOKER_ONLY: smoker, no diabetes
nhanes_like$group <- with(nhanes_like, ifelse( Diabetes == "Yes" & SmokeNow == "Yes",  "DIABETIC_SMOKER",
                      ifelse(Diabetes == "Yes" & SmokeNow != "Yes", "DIABETIC_NONSMOKER",
                      ifelse( Diabetes != "Yes" & Smoke100 == "No", "HEALTHY", "SMOKER_ONLY" )) ))
#### Drop rows where group can't be determined (NAs in key variables) ####
nhanes_like <- nhanes_like[!is.na(nhanes_like$group), ]

#########################################################################
#### Synthetic NIS emergency general surgery dataset (egs_misc_1621) ####
#########################################################################
#### Stands in for egs_misc_1621.csv (NIS 2016-2021): the LASSO section's outcome, predictors and survey design   ####
#### columns. Prevalences and effects mirror the real run: SUD always carries an Elixhauser drug or alcohol code, ####
#### mortality is about 3.7%, hospitals are nested in year-specific strata, and a few strata hold one hospital.   ####
#### Strata hold about 12 hospitals each, so the workflow can hold out hospitals within strata for validation.      ####
make_egs <- function(N, hospitals_per_year, n_strata, n_singleton_strata, hospital_id_start) {
  years <- 2016:2021
  codes <- sample(1000:9999, n_strata + n_singleton_strata)            # stratum codes reused every year, as in the NIS
  strata <- c(rep(codes[seq_len(n_strata)], length.out = hospitals_per_year - n_singleton_strata),
              codes[n_strata + seq_len(n_singleton_strata)])
  hosp <- data.frame(YEAR = rep(years, each = hospitals_per_year),
                     HOSP_NIS = hospital_id_start + seq_len(length(years) * hospitals_per_year),
                     NIS_STRATUM = as.vector(replicate(length(years), sample(strata))),
                     hospital_effect = rnorm(length(years) * hospitals_per_year, 0, 0.25))
  d <- hosp[sample(nrow(hosp), N, replace = TRUE), ]
  one_hot <- function(levels, prob) {
    g <- sample(levels, N, replace = TRUE, prob = prob)
    setNames(as.data.frame(sapply(levels, function(l) as.integer(g == l))), levels)
  }
  race <- one_hot(c("white_01", "hispanic_01", "black_01", "other_race_01", "asian_pi_01", "native_american_01"),
                  c(0.64, 0.14, 0.13, 0.04, 0.045, 0.005))
  payer <- one_hot(c("medicare_01", "private_01", "medicaid_01", "other_payor_01"), c(0.40, 0.31, 0.19, 0.10))
  income <- one_hot(c("first_quartile_01", "second_quartile_01", "third_quartile_01", "fourth_quartile_01"), c(0.29, 0.26, 0.24, 0.21))
  sud <- rbinom(N, 1, 0.062)
  age <- round(pmax(18, pmin(90, rnorm(N, ifelse(sud == 1, 43, 53), 18))))
  #### SUD patients always have an Elixhauser drug or alcohol code; a few non-SUD patients also have one ####
  drug <- ifelse(sud == 1, rbinom(N, 1, 0.61), rbinom(N, 1, 0.0004))
  alcohol <- ifelse(sud == 1, rbinom(N, 1, 0.47), rbinom(N, 1, 0.0022))
  neither <- sud == 1 & drug == 0 & alcohol == 0
  drug[neither] <- rbinom(sum(neither), 1, 0.55); alcohol[neither] <- 1 - drug[neither]
  elix_prev <- c(chf = 0.08, carit = 0.12, valv = 0.04, pcd = 0.02, pvd = 0.04, hypunc = 0.33, hypc = 0.16, para = 0.02, ond = 0.05,
                 cpd = 0.15, diabunc = 0.11, diabc = 0.10, hypothy = 0.10, rf = 0.11, ld = 0.06, pud = 0.02, aids = 0.003, lymph = 0.01,
                 metacanc = 0.03, solidtum = 0.03, rheumd = 0.03, coag = 0.05, obes = 0.18, wloss = 0.06, fed = 0.25, blane = 0.02,
                 dane = 0.03, psycho = 0.02, depre = 0.12)
  elix <- sapply(elix_prev, function(p) rbinom(N, 1, pmin(0.95, p * exp(0.02 * (age - 53)))))   # comorbidities rise with age
  colnames(elix) <- paste0("elix_", names(elix_prev))
  out <- data.frame(d[, c("YEAR", "HOSP_NIS", "NIS_STRATUM")], sud_any_01 = sud, race, payer, income,
                    PRDAY_first_egs_proc_positive = rnbinom(N, size = 1, mu = 1.3), sex_01 = rbinom(N, 1, 0.50), AGE = age,
                    elix, elix_alcohol = alcohol, elix_drug = drug,
                    DISCWT = round(runif(N, 4.85, 5.05), 5), row.names = NULL)
  #### Log-odds of death, using coefficients close to the real LASSO fit (SUD acts only through elix_drug/elix_alcohol) ####
  beta <- c(hispanic_01 = -0.34, black_01 = 0.16, other_race_01 = 0.04, asian_pi_01 = -0.41, medicare_01 = 0.24, medicaid_01 = 0.22,
            other_payor_01 = 0.31, first_quartile_01 = 0.33, second_quartile_01 = 0.19, third_quartile_01 = 0.08,
            PRDAY_first_egs_proc_positive = 0.023, sex_01 = 0.09, AGE = 0.034, elix_chf = 0.52, elix_carit = 0.73, elix_valv = -0.36,
            elix_pcd = 0.48, elix_pvd = 0.98, elix_hypunc = -0.62, elix_hypc = -0.42, elix_para = 0.98, elix_ond = 1.52, elix_cpd = 0.32,
            elix_diabunc = -0.07, elix_diabc = 0.08, elix_hypothy = -0.21, elix_rf = 0.32, elix_ld = 0.92, elix_pud = -0.23,
            elix_aids = 0.16, elix_lymph = 0.28, elix_metacanc = 0.51, elix_solidtum = 0.10, elix_rheumd = 0.21, elix_coag = 1.40,
            elix_obes = -0.34, elix_wloss = 0.66, elix_fed = 1.45, elix_blane = -0.20, elix_dane = -0.48, elix_alcohol = -0.03,
            elix_drug = -0.27, elix_psycho = 0.29, elix_depre = -0.31)
  eta <- drop(as.matrix(out[, names(beta)]) %*% beta) + d$hospital_effect
  intercept <- uniroot(function(a) mean(plogis(a + eta)) - 0.037, c(-20, 5))$root   # overall mortality about 3.7%
  out$Death_01 <- rbinom(N, 1, plogis(intercept + eta))
  #### Further outcomes for the propensity-score section, with known direct SUD effects (conditional on the covariates): ####
  #### LOS rate ratio exp(0.15) = 1.16; respiratory OR exp(0.4) = 1.49; cardiac OR exp(0.2) = 1.22; renal none (null).     ####
  #### Age raises every outcome and SUD patients are younger, so unadjusted comparisons understate the SUD effects.        ####
  age_c <- out$AGE - 53
  out$LOS <- rnbinom(N, size = 2, mu = exp(1.5 + 0.01 * age_c + 0.3 * out$elix_fed + 0.2 * out$elix_chf + 0.15 * sud + 0.5 * d$hospital_effect))
  out$complic_respiratory <- rbinom(N, 1, plogis(-3.0 + 0.02 * age_c + 0.5 * out$elix_cpd + 0.4 * sud))
  out$complic_cardiac <- rbinom(N, 1, plogis(-3.2 + 0.03 * age_c + 0.6 * out$elix_chf + 0.6 * out$elix_carit + 0.2 * sud))
  out$complic_renal <- rbinom(N, 1, plogis(-3.5 + 0.02 * age_c + 0.9 * out$elix_rf))
  out[, c("Death_01", setdiff(names(out), "Death_01"))]
}
set.seed(42)
egs_misc_1621 <- make_egs(N = 40000, hospitals_per_year = 150, n_strata = 12, n_singleton_strata = 4, hospital_id_start = 10000)

################################################################################
#### Save every example dataset into inst/extdata/*.rda ####
#### Save as data.table so the demo can use data.table syntax directly. Public ####
#### package functions also coerce ordinary data.frames at their boundaries.   ####
################################################################################
data.table::setDT(copd_data)
data.table::setDT(dfMigrant)
data.table::setDT(cardiac_data)
data.table::setDT(dream_data)
data.table::setDT(peru_data)
data.table::setDT(MData)
data.table::setDT(nhanes_like)
data.table::setDT(egs_misc_1621)

dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)
for (nm in c("copd_data", "dfMigrant", "cardiac_data", "dream_data", "peru_data", "MData", "nhanes_like", "egs_misc_1621")) {
  save(list = nm, file = file.path("inst/extdata", paste0(nm, ".rda")), compress = "bzip2")
}
