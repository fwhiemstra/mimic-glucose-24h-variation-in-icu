# Run-Time Environment: 	R version 4.3.1
# Authors:				        Floor Hiemstra & Laura Kervezee
# Date created:           January 2024
# Project:      			    24-hour glucose variation in ICU patients
# Filename:               mimic-iv-24h-glucose-variation_01_preprocessing
# Purpose:  			        Extraction and preprocessing for the time-stamped dataframe including:
#                         - Glucose data and info
#                         - Nutrition data
#                         - Patient & ICU stay info 
#                         - Dextrose, insulin, glucocorticoids administrations
#                         - RASS scores
#                         - Ventilation mode 
# Datafiles used:			    MIMIC-IV SQL database


########################################################################################################################
##### Packages -----
library(RPostgreSQL)
library(tidyverse)
library(stringr)
library(lubridate)
library(dplyr)
library(gtsummary)
library(ivs)

########################################################################################################################
##### STEP 1: Connect with SQL database -----
mimic_icu <- dbConnect(RPostgres::Postgres(),
                       host     = "localhost",
                       dbname   = "mimic",
                       user     = "postgres",
                       password = "postgres",
                       bigint   = "integer",
                       port="5432",
                       options  = "-c search_path=mimiciv_icu")

mimic_hosp <- dbConnect(RPostgres::Postgres(),
                        host     = "localhost",
                        dbname   = "mimic",
                        user     = "postgres",
                        password = "postgres",
                        bigint   = "integer",
                        port="5432",
                        options  = "-c search_path=mimiciv_hosp")

mimic_derived <- dbConnect(RPostgres::Postgres(),
                           host     = "localhost",
                           dbname   = "mimic",
                           user     = "postgres",
                           password = "postgres",
                           bigint   = "integer",
                           port="5432",
                           options  = "-c search_path=mimiciv_derived")


########################################################################################################################
##### STEP 2: Subject selection -----

## Extract all ICUstays and all patient info 
df_icustays_all <- tbl(mimic_icu, "icustays") %>% collect()

# Inclusion flowchart numbers
df_flowchart_subjects <- data.frame(total_icustays_all = nrow(df_icustays_all), 
                                    total_patients_all = length(unique(df_icustays_all$subject_id)))

## Exclude stays with length of stay (LOS) < 4 days
df_icustays_subject_selection <- df_icustays_all %>% filter(los >= 4) %>% collect()

# Inclusion flowchart numbers
df_flowchart_subjects$total_incl_icustays_los4 <- nrow(df_icustays_subject_selection)
df_flowchart_subjects$total_incl_patients_los4 <- length(unique(df_icustays_subject_selection$subject_id))
df_flowchart_subjects$total_excl_icustays_los4_from_all <- df_flowchart_subjects$total_icustays_all - df_flowchart_subjects$total_incl_icustays_los4
df_flowchart_subjects$total_excl_patients_los4_from_all <- df_flowchart_subjects$total_patients_all - df_flowchart_subjects$total_incl_patients_los4
## Exclude readmission stays 
df_icustays_subject_selection <- df_icustays_subject_selection %>%  group_by(subject_id) %>% filter(intime == min(intime)) %>% ungroup()

# Inclusion flowchart numbers
df_flowchart_subjects$total_incl_icustays_los4_readm <- nrow(df_icustays_subject_selection)
df_flowchart_subjects$total_incl_patients_los4_readm <- length(unique(df_icustays_subject_selection$subject_id))
df_flowchart_subjects$total_excl_icustays_los4_readm <- df_flowchart_subjects$total_incl_icustays_los4 - df_flowchart_subjects$total_incl_icustays_los4_readm
df_flowchart_subjects$total_excl_patients_los4_readm <- df_flowchart_subjects$total_incl_patients_los4 - df_flowchart_subjects$total_incl_patients_los4_readm

## Create arrays of included stayids en patientids
stayids_incl_subject_selection <- df_icustays_subject_selection$stay_id
subjectids_incl_subject_selection <- df_icustays_subject_selection$subject_id

## WRAP-UP:
# df_icustays_all: dataframe with all ICU stays
# df_icustays_subject_selection: dataframe with selected ICU stays
# stayids_incl_subject_selection: array with stay_ids of included ICU stays
# subjectids_incl_subject_selection: array with subject_ids of included ICU stays


########################################################################################################################
##### STEP 3: Extract nutrition data -----
## Collect itemids of nutrition events 
dfItems <- tbl(mimic_icu, "d_items") %>% collect() # Table with definition of itemids (e.g. types of nutrition)
dfitems_nut <-  dfItems %>% filter(grepl("nutrition", tolower(category)) | label == "PO Intake") %>% select(itemid, label, category)
itemid_nut <- dfitems_nut %>% pull(itemid)

## Obtain all feeding events in selected ICU stays from 'inputevents' table 
df_nut_all <- tbl(mimic_icu, "inputevents") %>% 
  filter(stay_id %in% stayids_incl_subject_selection & itemid %in% itemid_nut & amount > 0) %>%  collect() %>% 
  # Add feeding label
  left_join(dfitems_nut %>% select(itemid, category, label), by="itemid") %>% 
  # Add icustay information
  left_join(df_icustays_subject_selection %>% select(stay_id, intime, outtime, first_careunit, last_careunit, los), by=c("stay_id")) %>% 
  group_by(stay_id) %>% 
  # Starttime and endtime as seq time (time since midnight of first day), so time of day info is preserved
  mutate(starttime_seq = as.numeric(difftime(starttime, min(as.Date(intime)), units="hours")),
         endtime_seq = as.numeric(difftime(endtime, min(as.Date(intime)), units="hours")),
         outtime_seq = as.numeric(difftime(outtime, min(as.Date(intime)), units="hours")),
         intime_seq = as.numeric(difftime(intime, min(as.Date(intime)), units="hours"))) %>% 
  # Add index to indicate if feeding periods are adjacent
  ungroup() %>% group_by(stay_id) %>% arrange(stay_id, starttime_seq) %>% 
  mutate(nut_intervalnr = c(0, cumsum(lead(starttime_seq) > cummax(endtime_seq))[-n()]))

## Create groups of feeding types during an icu stay (i.e. only enteral nutrition, only parenteral nutrition, only oral intake, combination)
## NB: Nutriton - Supplements is ignored in this classification 
df_nut_all <- df_nut_all%>% group_by(stay_id) %>%
  mutate(nutritiongroup = case_when(all(category %in% c("Fluids/Intake", "Nutrition - Supplements")) ~ "Oral intake only",
                                    all(category %in% c("Nutrition - Enteral", "Nutrition - Supplements")) ~ "Enteral feeding only",
                                    all(category%in% c("Nutrition - Parenteral", "Nutrition - Supplements")) ~ "Parenteral feeding only",
                                    !(all(category %in% c("Nutrition - Enteral", "Nutrition - Supplements")) | all(category %in% c("Nutrition - Parenteral", "Nutrition - Supplements")) | all(ordercategoryname %in% c("Fluids/Intake", "Nutrition - Supplements"))) ~ "Combined"))

## Select patients that receive enteral nutrition only or that received enteral nutrition during their stay (i.e. before or after oral intake)
df_nut_enteral <- df_nut_all %>%
  filter(nutritiongroup %in% c("Combined", "Enteral feeding only")) %>%
  arrange(stay_id, starttime_seq) %>%
  group_by(stay_id, nut_intervalnr) %>%
  filter(all(category %in% c("Nutrition - Enteral", "Nutrition - Supplements")))

## Rename enteral nutrition
df_nut_enteral <- df_nut_enteral %>% rename(starttime_nut = starttime, endtime_nut=endtime, starttime_nut_seq = starttime_seq, endtime_nut_seq=endtime_seq, 
                                                   amount_nut = amount, label_nut=label, itemid_nut = itemid)

## Create arrays of included stayids
stayids_incl_nut <- unique(df_nut_enteral$stay_id) 

## Add nutrition selection to inclusion flowchart dataframe 
df_flowchart_subjects$total_incl_icustays_entfeed_frompatientselection <- length(stayids_incl_nut)
df_flowchart_subjects$total_excl_icustays_noentfeed_frompatientselection <- nrow(df_icustays_subject_selection) - df_flowchart_subjects$total_incl_icustays_entfeed_frompatientselection

# WRAP-UP: 
# df_nut_all: dataframe with all feeding episodes of included icustays (stayids_incl_subject_selection)
# df_nut_enteral: dataframe with all enteral feeding episodes of included icustays (stayids_incl_subject_selection)
# stayids_incl_nut: array with stay_ids of included subjects (stayids_incl_subject_selection) who received enteral nutrition during their ICU stay


########################################################################################################################
##### STEP 4: Extract glucose data -----
## Glucose itemids from chartevents table (both lab and fingerstick is used):
#225664, -- Glucose finger stick
#220621, -- Glucose (serum)
#226537, -- Glucose (whole blood)

## Extract all glucose events from patients in df_icustay_selection
df_gluc_all <- tbl(mimic_icu, "chartevents") %>% 
  filter(itemid %in% c(225664,220621,226537) & stay_id %in% stayids_incl_subject_selection) %>% collect() %>% 
  # Add patient information 
  left_join(df_icustays_subject_selection %>% select(stay_id, intime, outtime, los), by=c("stay_id")) %>% 
  # Only include values taken during icu stay
  filter(charttime >= intime & charttime <= outtime) %>%
  # From Robles Arevalo et al (2021): "Sometimes the STORETIME (time listed by nurses for checking glucose) was recorded
  # before the CHARTTIME (the time when the actual data entry occurred). In that case, the STORETIME timestamp 
  # was considered to be the time when the glycemic check occurred. Otherwise, the CHARTTIME timestamp was 
  # maintained as the time of glycemic check." The same strategy is followed here. 
  mutate(charttime = case_when(storetime < charttime ~ storetime,
                               storetime >= charttime ~ charttime))

## Add number of ICU stays without glucose measurment to df_flowchart_subjects
df_flowchart_subjects$total_incl_icustays_glucm_frompatientselection <- length(unique(df_gluc_all$stay_id))
df_flowchart_subjects$total_excl_icustays_noglucm_frompatientselection <- nrow(df_icustays_subject_selection) - df_flowchart_subjects$total_incl_icustays_glucm_frompatientselection

## Select glucose measurements from final included ICU stays (selected subjects that receive EN)
df_gluc_included_icustays <- df_gluc_all %>% filter(stay_id %in% stayids_incl_nut)

## Add to subject flowchart
df_flowchart_subjects$total_final_included_icustays <- length(unique(df_gluc_included_icustays$stay_id))

## Create df_flowchart_glucose with included glucose measurements from included patients
df_flowchart_glucose_m <- data.frame(total_glucose_m = nrow(df_gluc_included_icustays),
                                     total_glucose_m_patients = nrow(df_gluc_included_icustays %>% distinct(stay_id)))

## Remove in duplicate and valid duplicate glucose measurements
df_gluc_valid <- df_gluc_included_icustays %>% 
  # Remove duplicate rows; keeping serum > whole blood > fingerstick if there's a duplicated chart and lab value (i.e. max. itemid)
  group_by(subject_id, stay_id, charttime) %>% 
  mutate(item_preferred = case_when(itemid == 225664 ~ 3, #fingerstick is least preferred (more variable than other methods)
                                    itemid == 226537 ~ 2,
                                    itemid == 220621 ~ 1)) %>% 
  filter(item_preferred == min(item_preferred)) %>% 
  distinct(charttime, .keep_all=T) %>%
  select(-item_preferred) %>% ungroup() %>%
  # Remove invalid measurements; >1000 mg/dL for serum&blood, <500 mg/dL for fingerstick
  filter((valuenum <= 1000 & itemid %in% c(220621,226537)) | (valuenum <= 500 & itemid == 225664))

## Count duplicates for flowchart
df_non_dupl <- df_gluc_included_icustays %>% 
  # Remove duplicate rows; keeping serum > whole blood > fingerstick if there's a duplicated chart and lab value (i.e. max. itemid)
  group_by(subject_id, stay_id, charttime) %>% 
  mutate(item_preferred = case_when(itemid == 225664 ~ 3, #fingerstick is least preferred (more variable than other methods)
                                    itemid == 226537 ~ 2,
                                    itemid == 220621 ~ 1)) %>% 
  filter(item_preferred == min(item_preferred)) %>% 
  distinct(charttime, .keep_all=T) %>%
  select(-item_preferred) %>% ungroup()

df_flowchart_glucose_m$total_duplicates_glucose_m <- nrow(df_gluc_included_icustays) - nrow(df_non_dupl)

## Count invalids for flowchart 
df_flowchart_glucose_m$total_invalid_glucose_m <- nrow(df_gluc_included_icustays %>% filter(valuenum > 1000 & itemid %in% c(220621,226537))) + nrow(df_gluc_included_icustays %>% filter(valuenum > 500 & itemid == 225664))

## Total non-duplicate & valid glucose measurements
df_flowchart_glucose_m$total_valid_glucose_m <- nrow(df_gluc_valid)
df_flowchart_glucose_m$total_valid_glucose_m_icustays <- length(unique(df_gluc_valid$stay_id))

## Add new variables 
df_gluc_valid <- df_gluc_valid %>% group_by(stay_id) %>% 
  ## Add time variables 
  mutate(time_dec = as.numeric(local_time(charttime, tz="UTC", units="hours")),
         time_hbin = floor(time_dec),
         time_seq = as.numeric(difftime(charttime, as.Date(min(intime)), units="hours")),
         time_icudays = time_seq %/% 24) %>% 
  ungroup()

df_gluc_valid <- df_gluc_valid %>% group_by(subject_id, stay_id) %>% 
  ## Sort by time_seq to ensure correct calculation of time to next measurement 
  arrange((time_seq), .by_group = TRUE) %>%
  ## Add time since last measurement variable 
  mutate(time_last_measurement = time_seq - lag(time_seq, 1),
         time_next_measurement = lead(time_seq, 1) - time_seq) %>%
  ungroup() %>%

  ## Add sample type names 
  group_by(stay_id, charttime) %>%
  mutate(sample_type = case_when(itemid == 225664 ~ 'Fingerstick', 
                                 itemid == 226537 ~ 'Lab (serum)',
                                 itemid == 220621 ~ 'Lab (whole blood)')) %>%
  ungroup()

## WRAP-UP: 
# df_gluc_all: glucose values of all included icustays in df_icustays_subject_selection (with and without EN)
# df_gluc_included_icustays: glucose values of all included icustays in df_icustays_subject_selection with EN
# df_gluc_valid: Valid and non-duplicated glucose measurements of included icustays with EN with time variables)


########################################################################################################################
##### STEP 5: Link glucose & nutrition dataframe and extract glucose measurements taken DURING enteral feeding -----

## Link nutrition data to glucose data 
df_gluc_nutinfo <- df_gluc_valid %>% 
  left_join(df_nut_enteral %>% select(stay_id, starttime_nut_seq, endtime_nut_seq, amount_nut, label_nut, nut_intervalnr, itemid_nut), 
            by=c("stay_id"), relationship = "many-to-many") %>% 
  filter(time_seq > starttime_nut_seq & time_seq <= endtime_nut_seq) %>% 
  select(stay_id, time_seq, starttime_nut_seq, endtime_nut_seq, amount_nut, label_nut, nut_intervalnr, itemid_nut)

## Add enteral nutrition info to glucose data frame and only include glucose values taken during enteral feeding
df_gluc_nut <- df_gluc_valid %>% left_join(df_gluc_nutinfo, by=c("stay_id","time_seq")) %>% 
  arrange(stay_id, time_seq) %>%
  # Add enteral vs no enteral feeding label and rate of nutrition
  mutate(entfeeding_yn = ifelse(is.na(starttime_nut_seq), "n", "y"),
         rate_nut = amount_nut/(endtime_nut_seq - starttime_nut_seq)) %>% 
  # Only include glucose values during enteral feeding
  filter(entfeeding_yn == "y")

## Calculate carbohydrate rate per hour (gram/hour) (See Supplementary Materials for sources)
## First, check which itemids for enteral nutrition events appear in the dataset. For these itemids (nutrition types), CHO content is determined. 
unique_values_nutrition_events <- df_gluc_nut %>%
  distinct(label_nut, itemid_nut) %>%
  arrange(label_nut, itemid_nut)  # Sort alphabetically by 'label_nut' and 'itemid_nut'

## Now calculate CHO rate per hour (gram/hour) using itemid, rate (ml/hour) and the CHO content (gram/ml)
## NB: CHO content is calculated for itemids present in our dataset (df_gluc_nut)
df_gluc_nut <- df_gluc_nut %>% 
  mutate(rate_nut_cho = case_when(
    itemid_nut %in% c(225970, 229583) ~ rate_nut * 0, # Beneprotein
    itemid_nut %in% c(227979) ~ rate_nut * 0.067510549, # Boost Glucose Control
    itemid_nut %in% c(228355) ~ rate_nut * 0.189873418, # Enlive
    itemid_nut %in% c(226875, 225937) ~ rate_nut * 0.139240506, # Ensure
    itemid_nut %in% c(226877) ~ rate_nut * 0.202531646, # Ensure Plus
    itemid_nut %in% c(229574) ~ rate_nut * NA, # Fiber Supplement (i.e. Banana Flakes)
    itemid_nut %in% c(227698, 227699, 227696, 227695) ~ rate_nut * 0.156118143, # Fibersource HN
    itemid_nut %in% c(228356, 228359) ~ rate_nut * 0.11, # Glucerna
    itemid_nut %in% c(229013) ~ rate_nut * 0.113924050632911, # Glucerna 1.2
    itemid_nut %in% c(229295) ~ rate_nut * 0.132911392405063, # Glucerna 1.5
    itemid_nut %in% c(226023, 226020, 226022, 221207) ~ rate_nut * 0.132, # Impact
    itemid_nut %in% c(226027, 225928) ~ rate_nut * 0.132, # Impact with Fiber
    itemid_nut %in% c(228131, 228132, 228133, 228134, 228135) ~ rate_nut * 0.167088608, # Isosource 1.5
    itemid_nut %in% c(229010) ~ rate_nut * 0.169620253164557, # Jevity 1.2
    itemid_nut %in% c(229011) ~ rate_nut * 0.215611814345992, # Jevity 1.5
    itemid_nut %in% c(228348, 228351) ~ rate_nut * 0.147272727, # Nepro
    itemid_nut %in% c(227973, 227974, 227975) ~ rate_nut * 0.184810127, # NovaSource Renal 
    itemid_nut %in% c(226019, 226016, 226017, 227518, 225931) ~ rate_nut * 0.196, # Nutren 2.0
    itemid_nut %in% c(226882, 226881, 226880) ~ rate_nut * 0.10, # Nutren Pulmonary
    itemid_nut %in% c(226031, 226028, 226030, 221036) ~ rate_nut * 0.184810127, # Nutren Renal
    itemid_nut %in% c(229297) ~ rate_nut * 0.20337552742616, # Osmolite 1.5 
    itemid_nut %in% c(226039, 226036, 226038, 225930) ~ rate_nut * 0.188, # Peptamen 1.5
    itemid_nut %in% c(228383) ~ rate_nut * 0.078, # Peptamen Bariatric 
    itemid_nut %in% c(225929) ~ rate_nut * NA, # Probalance
    itemid_nut %in% c(229009) ~ rate_nut * 0.131223628691983, # Promote
    itemid_nut %in% c(229014) ~ rate_nut * 0.139240506329114, # Promote with Fiber
    itemid_nut %in% c(228360, 228361, 228363) ~ rate_nut * 0.105485232, # Pulmocare
    itemid_nut %in% c(226047, 226044, 226046, 225935) ~ rate_nut * 0.112, # Replete
    itemid_nut %in% c(226051, 226048, 226049, 226050, 225936) ~ rate_nut * 0.124, # Replete with Fiber
    itemid_nut %in% c(228364, 228367) ~ rate_nut * 0.219, # Two Cal HN
    itemid_nut %in% c(229012) ~ rate_nut * 0.186497890295359, # Vital 1.5
    itemid_nut %in% c(229296) ~ rate_nut * 0.112658227848101, # Vital High Protein
    itemid_nut %in% c(226059, 225934) ~ rate_nut * 0.00256, # Vivonex
    TRUE ~ NA_real_  # If none of the conditions match, specify the default value
  ))


################################################################################################## 
## Remove duplicated measurements because of overlapping feeding events and sum cho rate
df_gluc_nut <- df_gluc_nut %>% group_by(subject_id, hadm_id, stay_id, valuenum, valueuom, charttime, storetime, time_seq, time_dec, time_hbin, 
                                         itemid, sample_type, intime, outtime, los, time_icudays, time_last_measurement, time_next_measurement) %>% 
  #distinct(valuenum, .keep_all = T)
  summarize(rate_nut_cho = sum(rate_nut_cho),
            rate_nut = sum(rate_nut),
            label_nut = paste(label_nut, collapse="_"),
            itemid_nut = paste(itemid_nut, collapse="_"), 
            starttime_nut_seq = min(starttime_nut_seq), 
            endtime_nut_seq = max(endtime_nut_seq),
            .groups = "drop") %>%
  select(everything())

## Update glucose measurement flowchart
df_flowchart_glucose_m$total_noENduring_glucose_m <- nrow(df_gluc_valid) - nrow(df_gluc_nut) 

## Remove NA values (missing CHO content)
n_gluc_measurements_missing_cho <- sum(is.na(df_gluc_nut$rate_nut_cho)) # Save for flowchart
n_gluc_measurements_zero_cho <- df_gluc_nut %>% filter(rate_nut_cho == 0) %>% nrow() # Save for flowchart
df_gluc_nut <- df_gluc_nut %>%
  filter(!is.na(rate_nut_cho) & rate_nut_cho != 0)

## Extract ids of current selection
stayids_incl_gluc_selection <- unique(df_gluc_nut$stay_id)
subjectids_incl_gluc_selection <- unique(df_gluc_nut$subject_id)
hadmids_incl_gluc_selection <- unique(df_gluc_nut$hadm_id)

## Create dataframe including all data needed for analysis and rename columns
df_gluc_nut_variables <- df_gluc_nut %>% select(subject_id, hadm_id, stay_id, glucose_value=valuenum, glucose_unit=valueuom, charttime, storetime, time_seq, time_dec, time_hbin, 
                                                glucose_itemid=itemid, sample_type, intime, outtime, los, time_icudays, time_last_measurement, time_next_measurement,
                                                label_nut, rate_nut, rate_nut_cho, itemid_nut, starttime_nut_seq, endtime_nut_seq) %>% 
  ## Convert glucose value to mmol/L
  mutate(glucose_value_mmolL = glucose_value/18.0182)
  
## Update glucose measurement flowchart
df_flowchart_glucose_m$total_gluc_m_during_gluc_m_missing_cho<- n_gluc_measurements_missing_cho
df_flowchart_glucose_m$total_gluc_m_during_gluc_m_zero_cho<- n_gluc_measurements_zero_cho
df_flowchart_glucose_m$total_noENduring_glucose_m_patients <- length(stayids_incl_nut) - length(stayids_incl_gluc_selection)
df_flowchart_glucose_m$total_final_glucose_m <- nrow(df_gluc_nut_variables)
df_flowchart_glucose_m$total_final_glucose_m_patients <- length(stayids_incl_gluc_selection)

## WRAP-UP: 
# df_gluc_nut_variables: time-stamped dataframe with glucose & nutrition data with only glucose measurements during EN
# stayids_incl_gluc_selection: array with stay_ids of final included subject (who have glucose measurements during enteral nutrition)
# subjectids_incl_gluc_selection: array with subject_ids of final included subject (who have glucose measurements during enteral nutrition)                                                                                               
# hadmids_incl_gluc_selection: array with hamd_ids of final included subject (who have glucose measurements during enteral nutrition)                                                                                                 


########################################################################################################################
##### STEP 6: Export flowchart dataframes to csv and array of included subject_ids/hadm_ids/stay_ids -----

write_csv(df_flowchart_subjects, file=paste0("2-output/mimic-iv_24h-glucose-variation_01_flowchart-data-subject-inclusion.csv"))
write_csv(df_flowchart_glucose_m, file=paste0("2-output/mimic-iv_24h-glucose-variation_01_flowchart-data-glucose-measurements-inclusion.csv"))

df_included_ids <- data.frame(stayids_incl_gluc_selection, subjectids_incl_gluc_selection, hadmids_incl_gluc_selection)
write_csv(df_flowchart_subjects, file=paste0("2-output/mimic-iv_24h-glucose-variation_01_included-subject-hadm-stay-ids.csv"))


########################################################################################################################
##### STEP 7: Add patient characteristics to dataframe -----
## (age, sex, race, in-hospital mortality, diabetes diagnosis, comorbidity index, admission type, SOFA and OASIS scores )

### Extract and prepare relevant tables -----
## Extract patient info table
df_pat_info_gluc_selection <- tbl(mimic_hosp, "patients") %>% filter(subject_id %in% subjectids_incl_gluc_selection) %>% collect()

## Extract ICU stay info table
df_icustays_gluc_selection <- df_icustays_subject_selection  %>% filter(stay_id %in% stayids_incl_gluc_selection)

## Extract admission info table
df_adm_gluc_selection <- tbl(mimic_hosp, "admissions") %>% filter(hadm_id  %in% hadmids_incl_gluc_selection) %>% collect() %>% 
  select(subject_id, hadm_id, hosp_admittime=admittime, hosp_dischtime=dischtime, hosp_admission_type = admission_type, hosp_discharge_location = discharge_location)

## Combine patient, ICU stay and admission info tables 
df_general_patient_info <- df_pat_info_gluc_selection %>% 
  left_join(df_icustays_gluc_selection %>% select(subject_id, stay_id, intime, outtime, los), by="subject_id") %>% 
  left_join(df_adm_gluc_selection, by="subject_id")  

## Intime & outtime in hours since midnight of first day of ICU admission (faster downstream) 
df_general_patient_info <- df_general_patient_info %>% 
  mutate(intime_seq = as.numeric(difftime(intime, as.Date(intime), units="hours")),
         outtime_seq = as.numeric(difftime(outtime, as.Date(intime), units="hours")))


### Age and sex -----
df_general_patient_info <- df_general_patient_info %>% 
  ## Recalculate age because age is dependent on ICU stay time 
  mutate(age = as.numeric(format(intime, "%Y")) - anchor_year + anchor_age) %>% 
  select(-c(anchor_age, anchor_year))

## Add to df_gluc_nut_variables
df_gluc_nut_variables <- df_gluc_nut_variables %>% 
  left_join(df_general_patient_info %>% select(subject_id, gender, age), by="subject_id") %>%
  rename(sex=gender)

## Center age by the average age of all subjects
#age_per_subject <- df_gluc_nut_variables %>%
#  group_by(stay_id) %>%
#  summarize(age = mean(age, na.rm = TRUE))
#mean_age <- mean(age_per_subject$age)
#df_gluc_nut_variables <- df_gluc_nut_variables %>% 
#  mutate(age_centered = age - mean_age)

### Race -----
## Extract race data from admissions table 
admissions <- tbl(mimic_hosp, "admissions")  %>% 
  collect() %>% 
  #race is often coded in two parts, e.g. "ASIAN - KOREAN" --> only keep first part before hyphen using separate function
  separate(race, c("race_abbrv", NA), sep=" -", fill="right", remove=F) %>% 
  #categorize into Asian,  Black, Hispanic, White, Other/Unknown as in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10524813/; published code: https://github.com/joamats/mit-sepsis-tx/blob/master/src/1_queries/main.SQL
  mutate(race_cat = case_when(race_abbrv %in% c("HISPANIC/LATINO","SOUTH AMERICAN","HISPANIC OR LATINO") ~ "Hispanic", 
                              race_abbrv %in% c("BLACK/AFRICAN AMERICAN", "BLACK/CAPE VERDEAN", "BLACK/AFRICAN", "BLACK/CARIBBEAN ISLAND") ~ "Black",
                              race_abbrv %in% c("PORTUGUESE", "WHITE") ~ "White",
                              race_abbrv %in% c("ASIAN") ~ "Asian",
                              race_abbrv %in% c("UNABLE TO OBTAIN", "UNKNOWN", "PATIENT DECLINED TO ANSWER") ~ 'Unknown', 
                              race_abbrv %in% c("OTHER","AMERICAN INDIAN/ALASKA NATIVE","AMERICAN INDIAN OR ALASKA NATIVE", "NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER", "MULTIPLE RACE/ETHNICITY") ~ "Other",
                              TRUE ~ race_abbrv))

## Add to df_gluc_nut_variables
df_gluc_nut_variables <- df_gluc_nut_variables %>% 
  left_join(admissions %>% select(hadm_id, race_cat), by="hadm_id") 


### In-hospital mortality -----
df_general_patient_info <- df_general_patient_info %>% 
  ## Use hospital discharge info to determine whether patient died during hospital stay or sent to hospice
  mutate(diff_death = difftime(dod, as.Date(hosp_dischtime), units="days"),
         died_hospstay_yn = case_when(as.numeric(diff_death) <= 0 | hosp_discharge_location %in% c("HOSPICE","DIED") ~ "y",
                                      TRUE ~ "n")) %>% 
  select(-diff_death)

## Add to df_gluc_nut_variables
df_gluc_nut_variables <- df_gluc_nut_variables %>% 
  left_join(df_general_patient_info %>% select(subject_id, dod, died_hospstay_yn), by="subject_id")


### Diabetes -----
# Using icd codes from https://github.com/MIT-LCP/mimic-code/blob/main/mimic-iv/concepts/comorbidity/charlson.sql
dfdiag_diab <- tbl(mimic_hosp, "diagnoses_icd") %>% filter(hadm_id %in% hadmids_incl_gluc_selection) %>% 
  #icd codes from script linked above; ^ indicates that string should start with character string; | indicates OR
  filter(str_detect(icd_code, "^2500|^2501|^2502|^2503|^2508|^2509|^E100|^E10l|^E106|^E108|^E109|^E110|^E111|^E116|^E118|^E119|^E120|^E121|^E126|^E128|^E129|^E130|^E131|^E136|^E138|^E139|^E140|^E141|^E146|^E148|^E149|^2504|^2505|^2506|^2507|^E102|^E103|^E104|^E105|^E107|^E112|^E113|^E114|^E115|^E117|^E122|^E123|^E124|^E125|^E127|^E132|^E133|^E134|^E135|^E137|^E142|^E143|^E144|^E145|^E147")) %>% 
  collect()

## Add to df_gluc_nut_variables
df_gluc_nut_variables <- df_gluc_nut_variables %>% mutate(diabetes_status = ifelse(hadm_id %in% dfdiag_diab$hadm_id, "y", "n"))


### Charlson comorbidity index -----
# NB: CCI is available per hadm, not per icu stay
charlson <- tbl(mimic_derived, "charlson")  %>% 
  collect() %>% filter(hadm_id  %in% hadmids_incl_gluc_selection) %>% 
  #combine diabetes with & without chronic complications
  mutate(diabetes_charlson = ifelse(diabetes_without_cc == 1 | diabetes_with_cc == 1, 1, 0)) %>% 
  select(-c("diabetes_without_cc","diabetes_with_cc")) #note: group of patients with diabetes using charlson table is identical to what we already have in preprocessing script as the ICD codes were based on what is used to determine diabetes diagnosis for Charlson table

## Extract top 5 comorbidities
charlson_top5 <- charlson %>% select(-c("subject_id","hadm_id","age_score", "charlson_comorbidity_index")) %>% 
  summarize_all(sum) %>% pivot_longer(everything()) %>% arrange(desc(value)) %>% head(5) %>% 
  pull(name)

# Add to df_gluc_nut_variables
df_gluc_nut_variables <- df_gluc_nut_variables %>% 
  left_join(charlson %>% select(hadm_id, all_of(charlson_top5), charlson_comorbidity_index), by="hadm_id") 


### SOFA at admission (first day)-----
#SOFA is available per ICU stay
sofa <- tbl(mimic_derived, "first_day_sofa")  %>% 
  collect() %>% filter(stay_id  %in% stayids_incl_gluc_selection)

## Add to df_gluc_nut_variables
df_gluc_nut_variables <- df_gluc_nut_variables %>% 
  left_join(sofa %>% select(stay_id, sofa), by="stay_id") 

### Admission type (urgent vs elective) ----
# Mapping as in https://github.com/joamats/mit-sepsis-tx/blob/master/src/1_queries/main.SQL from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10524813/: 
# Emergency: "AMBULATORY OBSERVATION’, ‘DIRECT EMER.’, ‘URGENT’, ‘EW EMER.’, ‘DIRECT OBSERVATION’, ‘EU OBSERVATION’, ‘OBSERVATION ADMIT’
# Elective: ‘ELECTIVE’, ‘SURGICAL SAME DAY ADMISSION’

admissiontype <- tbl(mimic_hosp, "admissions")  %>% select(hadm_id, admission_type) %>% 
  collect() %>%  filter(hadm_id  %in% hadmids_incl_gluc_selection)

admissiontype %>% select(-hadm_id) %>% tbl_summary(sort = list(everything() ~ "frequency"))

## Map to elective or emergency
admissiontype <- admissiontype %>% mutate(admtype_cat = case_when(admission_type %in% c("AMBULATORY OBSERVATION", "DIRECT EMER.", "URGENT", "EW EMER.", "DIRECT OBSERVATION", "EU OBSERVATION", "OBSERVATION ADMIT") ~ "Emergency",
                                                                  admission_type %in% c("ELECTIVE", "SURGICAL SAME DAY ADMISSION") ~ "Elective",
                                                                  TRUE ~ NA_character_))
admissiontype %>% select(admtype_cat) %>% tbl_summary(sort = list(everything() ~ "frequency")) #no missing, all correctly mapped

## Add to df_gluc_nut_variables
df_gluc_nut_variables <- df_gluc_nut_variables %>% 
  left_join(admissiontype %>% select(hadm_id, admtype_cat), by="hadm_id") 


### OASIS ----
# Oxford Acute Severity of Illness Score (OASIS) is a novel severity of illness score that is determined upon admission to ICU, based on fewer parameters compared to APACHE IV but comparable performance (APACHE is not readily available in MIMIC-IV)
# ref: https://journals.lww.com/ccmjournal/fulltext/2013/07000/A_New_Severity_of_Illness_Scale_Using_a_Subset_of.15.aspx
oasis <- tbl(mimic_derived, "oasis")  %>% 
  collect() %>% filter(stay_id  %in% stayids_incl_subject_selection)

## Add to df_gluc_nut_variables
df_gluc_nut_variables <- df_gluc_nut_variables %>% 
  left_join(oasis %>% select(stay_id, oasis), by="stay_id") 


########################################################################################################################
##### STEP 8: Add dextrose administration to dataframe -----
items_dext <- dfItems %>% filter(grepl("dextrose", tolower(label)))  
itemid_dext <- items_dext$itemid

## Extract dextrose administration data of included subjects
df_dext <- tbl(mimic_icu, "inputevents") %>% filter(itemid %in% itemid_dext) %>% 
  # Only select events from included patients
  filter(stay_id %in% stayids_incl_gluc_selection) %>% filter(amount > 0) %>% collect() %>% 
  # Link selection from inputevents table to df_general_patient_info
  left_join(items_dext, by="itemid") %>% 
  left_join(df_general_patient_info %>% select(stay_id, intime), by=c("stay_id")) %>% 
  group_by(stay_id) %>% 
  # Change start- & endtimes to times 
  mutate(starttime_dext_seq = as.numeric(difftime(starttime, as.Date(min(intime)), units="hours")),
         endtime_dext_seq = as.numeric(difftime(endtime, as.Date(min(intime)), units="hours"))) %>% # # 
  arrange(subject_id, starttime) %>% 
  select(stay_id, starttime_dext_seq, endtime_dext_seq, rate, rateuom, itemid_dext = itemid, label_dext=label, ordercat_dext = ordercategorydescription, amount_total = amount) 

# Quantify dextrose administration for each concentration 
df_dext <- df_dext %>%
  mutate(rate_dext = case_when(
    label_dext == 'Dextrose 5%' & rateuom == 'mL/hour' & (ordercat_dext == 'Continuous Med' | ordercat_dext == 'Continuous IV')  ~ 0.05 * rate, # Rate is always in mL/hour for continuous administrations
    label_dext == 'Dextrose 10%' & rateuom == 'mL/hour' & (ordercat_dext == 'Continuous Med' | ordercat_dext == 'Continuous IV') ~ 0.1 * rate,
    label_dext == 'Dextrose 20%' & rateuom == 'mL/hour' & (ordercat_dext == 'Continuous Med' | ordercat_dext == 'Continuous IV') ~ 0.2 * rate,
    label_dext == 'Dextrose 40%' & rateuom == 'mL/hour' & (ordercat_dext == 'Continuous Med' | ordercat_dext == 'Continuous IV') ~ 0.4 * rate,
    label_dext == 'Dextrose 50%' & rateuom == 'mL/hour' & (ordercat_dext == 'Continuous Med' | ordercat_dext == 'Continuous IV') ~ 0.5 * rate,
    ## For drug push and bolus administrations, the period is extended to 10 minutes 
    label_dext == 'Dextrose 5%' & (ordercat_dext == 'Drug Push' | ordercat_dext == 'Bolus') ~ 0.05 * amount_total/0.1667,
    label_dext == 'Dextrose 10%' & (ordercat_dext == 'Drug Push' | ordercat_dext == 'Bolus') ~ 0.10 * amount_total/0.1667,
    label_dext == 'Dextrose 20%' & (ordercat_dext == 'Drug Push' | ordercat_dext == 'Bolus') ~ 0.20 * amount_total/0.1667,
    label_dext == 'Dextrose 40%' & (ordercat_dext == 'Drug Push' | ordercat_dext == 'Bolus') ~ 0.40 * amount_total/0.1667,
    label_dext == 'Dextrose 50%' & (ordercat_dext == 'Drug Push' | ordercat_dext == 'Bolus') ~ 0.50 * amount_total/0.1667,
    label_dext == 'Dextrose PN' ~ amount_total/(endtime_dext_seq - starttime_dext_seq ),   # Amount of dextrose in parenteral nutrition,
                          # (adapted from: https://github.com/oizin/glucose-data-driven-prediction/blob/main/sql/1-feeding.sql --> not present in our dataset)
    TRUE ~ NA_real_ # Fill with NAs if none of the above conditions are true --> All conditions are met, no Nas. 
    # NB: Unit of rate_dext is grams dextrose/hour 
  ))

## Link dextrose administration data to df_gluc_nut_variables dataframe
dfgluc_dextinfo <- df_gluc_nut_variables %>% ungroup() %>% left_join(df_dext, by=c("stay_id"), relationship="many-to-many") %>% 
  filter(time_seq > starttime_dext_seq & time_seq <= endtime_dext_seq) %>% ungroup() %>% 
  group_by(stay_id, time_seq) %>% 
  reframe(starttime_dext_seq  = min(starttime_dext_seq),
            endtime_dext_seq = max(endtime_dext_seq),
            rate_dext = sum(rate_dext),# Account for overlapping events
            label_dext = paste(label_dext, collapse="_"),
            ordercat_dext = paste(ordercat_dext, collapse="_"))

## Add dextrose administration data to df_gluc_nut_variables dataframe
df_gluc_nut_variables <- df_gluc_nut_variables %>% 
  left_join(dfgluc_dextinfo %>% ungroup(), by=c("stay_id", "time_seq")) %>%
  mutate(dext_yn = ifelse(!is.na(starttime_dext_seq), "y","n"),
         rate_dext = ifelse(is.na(rate_dext), 0, rate_dext))

########################################################################################################################
##### STEP 9: Add insulin administration to dataframe -----
items_insulin <- dfItems %>% filter(grepl("insulin", tolower(label)) & category=="Medications") %>%
  select(itemid, label)
itemid_insulin <- items_insulin$itemid

## Extract insulin administration data
df_insulin <- tbl(mimic_icu, "inputevents") %>% filter(itemid %in% itemid_insulin) %>% 
  filter(stay_id %in% stayids_incl_gluc_selection) %>%   filter(amount > 0) %>% collect() %>% 
  left_join(items_insulin, by="itemid") %>% 
  left_join(df_general_patient_info %>% select(stay_id, intime), by=c("stay_id")) %>% 
  group_by(stay_id) %>% 
  mutate(starttime_insu_seq = as.numeric(difftime(starttime, as.Date(min(intime)), units="hours")),
         endtime_insu_seq = as.numeric(difftime(endtime, as.Date(min(intime)), units="hours"))) %>% 
  arrange(stay_id, starttime_insu_seq) %>% 
  select(stay_id, starttime_insu_seq, endtime_insu_seq, itemid_insu = itemid, label_insu=label, ordercat_insu = ordercategorydescription, units_insu = amount) %>% 
  # Correction of start times of insulin (+0.5h) and end times of insulin depending on their mode of action; see e.g. ref. https://www.ncbi.nlm.nih.gov/pubmed/33692359) 
  mutate(
    endtime_insu_seq_corr = case_when(
      label_insu == "Insulin - Regular" ~ endtime_insu_seq + 4, #short
      label_insu == "Insulin - Humalog" ~ endtime_insu_seq + 4, #short
      label_insu == "Insulin - Novolog" ~ endtime_insu_seq + 2, #rapid
      label_insu == "Insulin - NPH" ~ endtime_insu_seq + 10, #intermediate
      label_insu == "Insulin - 70/30" ~ endtime_insu_seq + 10, #intermediate
      label_insu == "Insulin - Humalog 75/25"~ endtime_insu_seq + 10, #intermediate
      label_insu == "Insulin - U500" ~ endtime_insu_seq + 10, #intermediate
      label_insu == "Insulin - Glargine" ~ endtime_insu_seq + 12, #long
      TRUE ~ endtime_insu_seq),
    starttime_insu_seq_corr = case_when(label_insu == "Insulin - Regular"~ starttime_insu_seq + 0.5,
                                        label_insu == "Insulin - Humalog"~ starttime_insu_seq + 0.5,
                                        label_insu == "Insulin - Novolog"~ starttime_insu_seq + 0.5,
                                        label_insu == "Insulin - NPH"~ starttime_insu_seq  + 0.5,
                                        label_insu == "Insulin - 70/30"~ starttime_insu_seq + 0.5,
                                        label_insu == "Insulin - Humalog 75/25"~ starttime_insu_seq + 0.5,
                                        label_insu == "Insulin - U500"~ starttime_insu_seq + 0.5,
                                        label_insu == "Insulin - Glargine" ~ starttime_insu_seq + 0.5,
                                        TRUE ~ starttime_insu_seq)) %>% 
  relocate(endtime_insu_seq_corr, .after=endtime_insu_seq) %>% 
  relocate(starttime_insu_seq_corr, .after=starttime_insu_seq) %>% 
  # Add variable to indicate whether insulin administration is bolus or infusion
  mutate(insu_cat = case_when(ordercat_insu %in% c("Continuous Med","Continuous IV") ~ "infusion",
                              ordercat_insu %in% c("Drug Push") ~ "bolus"),
         # Calculate rate of insulin administration (units/hour)
         rate_insu = units_insu /(endtime_insu_seq - starttime_insu_seq))

## Link insulin administration data to df_gluc_nut_variables dataframe
dfgluc_insuinfo <- df_gluc_nut_variables %>% ungroup() %>% left_join(df_insulin, by=c("stay_id"), relationship="many-to-many") %>% 
  filter(time_seq > starttime_insu_seq_corr & time_seq <= endtime_insu_seq_corr) %>% ungroup() %>% 
  group_by(stay_id, time_seq) %>% 
  summarize(starttime_insu_seq_corr = min(starttime_insu_seq_corr),
            endtime_insu_seq_corr = max(endtime_insu_seq_corr)) %>%  
  ungroup() %>% 
  left_join(df_insulin, by=c("stay_id", "starttime_insu_seq_corr", "endtime_insu_seq_corr"), relationship="many-to-many") %>% 
  group_by(stay_id, time_seq, starttime_insu_seq_corr, endtime_insu_seq_corr) %>%
  reframe(insu_cat = paste(insu_cat, collapse="_"),
          rate_insu = sum(rate_insu)) # Account for overlapping events

## Add insulin administration data to df_gluc_nut_variables dataframe
df_gluc_nut_variables <- df_gluc_nut_variables %>% left_join(dfgluc_insuinfo %>% ungroup(), by=c("stay_id", "time_seq"),  relationship="many-to-many") %>%
  mutate(insu_yn = ifelse(!is.na(starttime_insu_seq_corr), "y","n"),
         rate_insu = ifelse(is.na(rate_insu), 0, rate_insu))  

## Add insulin requirement per day (units per day in ICU; measure of insulin resistance) -----
df_insu_req <- df_insulin %>% mutate(time_icudays = starttime_insu_seq %/% 24)  %>% 
  group_by(stay_id, time_icudays) %>% 
  summarize(unitsperday_insu = sum(units_insu)) %>% ungroup()

df_gluc_nut_variables <- df_gluc_nut_variables %>% left_join(df_insu_req, by=c("stay_id", "time_icudays")) %>% 
  mutate(unitsperday_insu = ifelse(is.na(unitsperday_insu), 0, unitsperday_insu))

########################################################################################################################
##### STEP 10: Add covariate RASS scores to dataframe -----
## Extract RASS scores data
df_ras <- tbl(mimic_icu, "chartevents") %>% 
  filter(itemid %in% c(228096) & stay_id %in% stayids_incl_gluc_selection) %>% collect() %>% 
  # Add intime_seq and outtime_seq to determine if values are registered during ICU stay 
  left_join(df_general_patient_info %>% select(subject_id, stay_id, intime, intime_seq, outtime_seq), by=c("stay_id", "subject_id")) %>% 
  #group_by(stay_id) %>% 
  mutate(charttime_seq = as.numeric(difftime(charttime, as.Date(intime), units="hours"))) %>% 
  # Only include values taken during icu stay
  filter(charttime_seq >= intime_seq & charttime_seq <= outtime_seq) %>% 
  arrange(stay_id, charttime_seq) %>% 
  group_by(stay_id) %>% 
  mutate(starttime_ras_seq = charttime_seq,
         endtime_ras_seq = lead(charttime_seq)) %>% 
  mutate(endtime_ras_seq = ifelse(is.na(endtime_ras_seq), outtime_seq, endtime_ras_seq))

## Combine periods with similar RAS scores to reduce size of table
dfras_int <- df_ras %>% group_by(stay_id) %>% arrange(stay_id, charttime_seq) %>% 
  mutate(intnr = cumsum(valuenum != lag(valuenum, default = 99))) %>% 
  group_by(stay_id, intnr) %>% 
  summarize(starttime_ras_seq = min(starttime_ras_seq),
            endtime_ras_seq = max(endtime_ras_seq),
            ras_score = mean(valuenum),
            ras_label = unique(value)) 

## Link RASS score data to df_gluc_nut_variables dataframe
dfgluc_rasinfo <- df_gluc_nut_variables %>% ungroup() %>% 
  left_join(dfras_int %>% select(stay_id, starttime_ras_seq, endtime_ras_seq, ras_score, ras_label), by=c("stay_id"), relationship="many-to-many") %>% 
  filter(time_seq > starttime_ras_seq & time_seq <= endtime_ras_seq) 

## Add RASS score data to df_gluc_nut_variables dataframe
df_gluc_nut_variables <- df_gluc_nut_variables %>% ungroup() %>% 
  left_join(dfgluc_rasinfo %>% select(stay_id, time_seq, ras_score, ras_label), 
            by=c("stay_id", "time_seq"), relationship="many-to-many") %>% 
  arrange(stay_id, time_seq) %>% 
  distinct(stay_id, time_seq, .keep_all = TRUE)


########################################################################################################################
##### STEP 12: Add ventilation mode to dataframe -----
## Extract ventilation data
dfvent <- tbl(mimic_derived, "ventilation") %>% collect() %>% 
  #add intime to compute sequential time of day
  left_join(df_icustays_subject_selection %>% select(stay_id, intime), by=c("stay_id")) %>% 
  group_by(stay_id) %>% 
  mutate(starttime_vent_seq = as.numeric(difftime(starttime, as.Date(min(intime)), units="hours")),
         endtime_vent_seq = as.numeric(difftime(endtime, as.Date(min(intime)), units="hours")))

## Link ventilation data to df_gluc_nut_variables dataframe
dfgluc_ventinfo <- df_gluc_nut_variables %>% 
  left_join(dfvent, by=c("stay_id"), relationship = "many-to-many") %>% 
  filter(time_seq > starttime_vent_seq & time_seq <= endtime_vent_seq) %>% ungroup()

## Add ventilation data to df_gluc_nut_variables dataframe
df_gluc_nut_variables <- df_gluc_nut_variables %>% ungroup() %>% 
  left_join(dfgluc_ventinfo %>% select(stay_id, time_seq, ventilation_status), by=c("stay_id", "time_seq"), relationship="many-to-many") %>% 
  arrange(stay_id, time_seq) %>% 
  # Reduce categories (from MIMIC code ventilation.sql); categories changed since MIMIC-IV 2.0
  mutate(ventilation_status_cat = case_when(ventilation_status %in% c("InvasiveVent", "Tracheostomy") ~ "invasive",
                                            ventilation_status %in% c("NonInvasiveVent", "HFNC", "SupplementalOxygen") ~ "noninvasive",
                                            ventilation_status %in% c("None") ~ "none", 
                                            is.na(ventilation_status) ~ "NaN")) 


########################################################################################################################
##### STEP 13: Add covariate corticosteroid administration to dataframe -----
# Corticosteroids include: dexamethasone, *cortisone, *prednisone, and *prednisolone (from: https://www.sciencedirect.com/science/article/pii/S0735675720300759)

# Use pharmacy table to determine when steroids were released by pharmacy (emar table that contains actual administrations is not available for all patients)
cort_pharm_tot <- tbl(mimic_hosp, "pharmacy") %>% filter(hadm_id %in% hadmids_incl_gluc_selection) %>%  
  filter(grepl("dexamethasone", tolower(medication)) | grepl("cortisone", tolower(medication)) | grepl("prednisolone", tolower(medication)) | grepl("prednisone", tolower(medication))) %>% 
  collect() %>% 
  # Add icustay information
  left_join(df_general_patient_info %>% select(stay_id, hadm_id, intime, intime_seq, outtime_seq), by=c("hadm_id")) %>% 
  group_by(hadm_id) %>% 
  mutate(starttime_seq = as.numeric(difftime(starttime, min(as.Date(intime)), units="hours")),
         stoptime_seq = as.numeric(difftime(stoptime, min(as.Date(intime)), units="hours"))) 

# Filter cort table
cort_pharm <- cort_pharm_tot %>% 
  #exclude glucocorticoids administred as study medication or when pharmacy_id is missing
  filter(!grepl("study", tolower(medication)) & !grepl("placebo", tolower(medication)) & !is.na(pharmacy_id)) %>% 
  #only include oral or IV administrations
  filter(route %in% c("IV","PO","ORAL","PO/NG")) %>% 
  #restrict data to administrations that occurred during ICU stay (stoptime should be not be before intime OR starttime should not be before outtime) OR Starttime is before intime and stoptime after outtime
  filter(!(stoptime_seq <= intime_seq | starttime_seq >= outtime_seq) | (starttime_seq <= intime_seq & stoptime_seq >= outtime_seq)) %>% 
  arrange(stay_id, starttime_seq) %>% 
  ungroup()


# Determine intervals during which glucocorticoids are active
cort_pharm <- cort_pharm %>% 
  #sometimes stoptime is after starttime; swap these so order of starttimes can be determined
  mutate(starttime_seq_corr = ifelse(starttime_seq > stoptime_seq, stoptime_seq, starttime_seq),
         stoptime_seq_corr = ifelse(starttime_seq > stoptime_seq, starttime_seq, stoptime_seq),
         #add 24h to stoptime to take into account (approximate) duration of action of glucocorticoids
         stoptime_seq_corr = stoptime_seq_corr + 24,
         #if starttime is negative (ie happens before start of icu stay), adjust to 0
         starttime_seq_corr = ifelse(starttime_seq_corr < 0, 0, starttime_seq_corr)) %>% 
  #determine intervals using ivs package
  mutate(iv = iv(starttime_seq_corr, stoptime_seq_corr)) %>%  
  group_by(stay_id)  %>%  
  mutate(iv = iv_identify_group(iv),
         iv_start = iv_start(iv),
         iv_end = iv_end(iv),
         cort_intervalno = cumsum(iv_start(iv) != lag(iv_start(iv), default=-Inf))) %>% 
  arrange(stay_id, starttime_seq_corr)  

cort_pharm_interval <- cort_pharm %>% group_by(stay_id, cort_intervalno) %>% 
  summarize(cort_interval_starttime_seq = unique(iv_start(iv)),
            cort_interval_endtime_seq = unique(iv_end(iv)))  

# Join with glucose dataframe
dfgluc_cortinfo <- df_gluc_nut_variables %>% ungroup() %>%  left_join(cort_pharm_interval, by=c("stay_id"), relationship="many-to-many") %>% 
  filter(time_seq > cort_interval_starttime_seq  & time_seq <= cort_interval_endtime_seq ) %>% ungroup() %>% 
  group_by(stay_id, time_seq) %>% 
  summarize(cort_interval_starttime_seq  = min(cort_interval_starttime_seq),
            cort_interval_endtime_seq = max(cort_interval_endtime_seq))

df_gluc_nut_variables <- df_gluc_nut_variables %>% left_join(dfgluc_cortinfo %>% ungroup(), by=c("stay_id", "time_seq")) %>%
  mutate(glucocorticoid_yn = ifelse(!is.na(cort_interval_starttime_seq), "y","n")) 


########################################################################################################################
##### STEP 14: Categorize continuous variables age, dextrose rate, insulin rate and cho rate -----

## Categorize age variable
# Check quantiles of continuous age variable (boundaries are based on data distribution)
quantiles <- quantile(df_gluc_nut_variables$age, probs = c(0, 0.25, 0.5, 0.75, 1))
quantiles

# Add categorized age variable
df_gluc_nut_variables$age_cat <- cut(df_gluc_nut_variables$age, 
                                     breaks = c(-Inf, 55, 65, 75, Inf), 
                                     include.lowest = TRUE)
## Check distribution of data points per category
table(df_gluc_nut_variables$age_cat)


## Categorize dextrose rate variable
# Check quantiles of dextrose rate (when rate_dext > 0) (boundaries are based on data distribution)
subset <- df_gluc_nut_variables[df_gluc_nut_variables$rate_dext > 0, ]
quantiles <- quantile(subset$rate_dext, probs = c(0, 0.33, 0.66, 1))
quantiles

# Add categorized dextrose rate variable
df_gluc_nut_variables$rate_dext_cat <- cut(df_gluc_nut_variables$rate_dext, 
                                           breaks = c(-Inf, 0, 0.5, 2, Inf), 
                                           include.lowest = TRUE)
## Check distribution of data points per category
table(df_gluc_nut_variables$rate_dext_cat)


## Categorize insulin rate variable
# Check quantiles of insulin rate (when rate_insu > 0) (boundaries are based on data distribution)
subset <- df_gluc_nut_variables[df_gluc_nut_variables$rate_insu > 0, ]
quantiles <- quantile(subset$rate_insu, probs = c(0, 0.33, 0.66, 1))
quantiles

# Add categorized insulin rate variable
df_gluc_nut_variables$rate_insu_cat <- cut(df_gluc_nut_variables$rate_insu, 
                                           breaks = c(-Inf, 0, 200, 750, Inf), 
                                           include.lowest = TRUE)
## Check distribution of data points per category
table(df_gluc_nut_variables$rate_insu_cat)


## Categorize nutrition carbohydrate (cho) rate variable
# Check quantiles of nutrition cho rate (boundaries are based on data distribution)
quantiles <- quantile(df_gluc_nut_variables$rate_nut_cho, probs = c(0, 0.25, 0.5, 0.75, 1))
quantiles

# Add categorized nutrition cho rate variable
df_gluc_nut_variables$rate_nut_cho_cat <- cut(df_gluc_nut_variables$rate_nut_cho, 
                                              breaks = c(-Inf, 4.5, 6.5, 8.5, Inf), 
                                              include.lowest = TRUE)
## Check distribution of data points per category
table(df_gluc_nut_variables$rate_nut_cho_cat)


########################################################################################################################
##### STEP 15: Write output to csv file -----
write_csv(df_gluc_nut_variables, file=paste0("2-output/mimic-iv_24h-glucose-variation_01_glucose-incl-covariates.csv"))
write_csv(df_insulin, file=paste0("2-output/mimic-iv_24h-glucose-variation_01_insulinadministrations.csv"))
write_csv(df_dext, file=paste0("2-output/mimic-iv_24h-glucose-variation_01_dextroseadministrations.csv"))
write_csv(df_nut_enteral, file=paste0("2-output/mimic-iv_24h-glucose-variation_01_enteralfeeding.csv"))
write_csv(df_gluc_all, file=paste0("2-output/mimic-iv_24h-glucose-variation_01_all_gluc_measurements.csv"))
write_csv(dfgluc_cortinfo, file=paste0("2-output/mimic-iv_24h-glucose-variation_01_cortadministrations.csv"))


