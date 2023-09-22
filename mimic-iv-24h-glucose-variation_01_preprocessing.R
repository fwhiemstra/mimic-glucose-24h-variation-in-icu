# Run-Time Environment: 	R version 4.0.3 & R Studio 1.4.1106
# Authors:				        Floor Hiemstra & Laura Kervezee
# Date created:           August 2023
# Project:      			    24-hour glucose variation in ICU patients
# Filename:               mimic-iv-24h-glucose-variation_01_preprocessing
# Purpose:  			        Extraction and preprocessing for the time-stamped dataframe including:
#                         - Glucose data and info
#                         - Nutrition data
#                         - Patient & ICU stay info 
#                         - Dextrose, insulin, glucocorticoids administration
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

df_flowchart_subjects <- data.frame(total_icustays_all = nrow(df_icustays_all), # Store inclusion flowchart data 
                           total_patients_all = length(unique(df_icustays_all$subject_id)))

## Exclude stays with length of stay (LOS) < 4 days
df_icustays_subject_selection <- df_icustays_all %>% filter(los >= 4) %>% collect()

df_flowchart_subjects$total_incl_icustays_los4 <- nrow(df_icustays_subject_selection)
df_flowchart_subjects$total_incl_patients_los4 <- length(unique(df_icustays_subject_selection$subject_id))
df_flowchart_subjects$total_excl_icustays_los4_from_all <- df_flowchart_subjects$total_icustays_all - df_flowchart_subjects$total_incl_icustays_los4
df_flowchart_subjects$total_excl_patients_los4_from_all <- df_flowchart_subjects$total_patients_all - df_flowchart_subjects$total_incl_patients_los4

## Exclude readmission stays 
df_icustays_subject_selection <- df_icustays_subject_selection %>%  group_by(subject_id) %>% filter(intime == min(intime)) %>% ungroup()

df_flowchart_subjects$total_incl_icustays_los4_readm <- nrow(df_icustays_subject_selection)
df_flowchart_subjects$total_incl_patients_los4_readm <- length(unique(df_icustays_subject_selection$subject_id))
df_flowchart_subjects$total_excl_icustays_los4_readm <- df_flowchart_subjects$total_incl_icustays_los4 - df_flowchart_subjects$total_incl_icustays_los4_readm
df_flowchart_subjects$total_excl_patients_los4_readm <- df_flowchart_subjects$total_incl_patients_los4 - df_flowchart_subjects$total_incl_patients_los4_readm

## Create arrays of included stayids en patientids
incl_stayids_subject_selection <- df_icustays_subject_selection$stay_id
incl_subjectids_subject_selection <- df_icustays_subject_selection$subject_id

## WRAP-UP:
# df_icustays_all: dataframe with all ICU stays
# df_icustays_subject_selection: dataframe with selected ICU stays
# incl_stayids_subject_selection: array with stay_ids of included ICU stays
# incl_subjectids_subject_selection: array with subject_ids of included ICU stays


########################################################################################################################
##### STEP 3: Extract nutrition data -----
## Collect itemids of nutrition events 
dfItems <- tbl(mimic_icu, "d_items") %>% collect() #table with definition of itemids (e.g. types of nutrition)
dfitems_nut <-  dfItems %>% filter(grepl("nutrition", tolower(category)) | label == "PO Intake") %>% select(itemid, label, category)
itemid_nut <- dfitems_nut %>% pull(itemid)

## Obtain all feeding events in selected ICU stays from 'inputevents' table 
df_nut_all <- tbl(mimic_icu, "inputevents") %>% 
  filter(stay_id %in% incl_stayids_subject_selection & itemid %in% itemid_nut & amount > 0) %>%  collect() %>% 
  # Add feeding label
  left_join(dfitems_nut %>% select(itemid, category, label), by="itemid") %>% 
  # Combine information from same order (= e.g. enteral feeding + added beneprotein; sum their amounts)
  group_by(stay_id, starttime, endtime, orderid, ordercategoryname, ordercategorydescription, patientweight) %>% 
  summarize(amount_combined = sum(amount), 
            label_combined = paste(label, collapse="_"),
            itemid_combined = paste(itemid, collapse="_"),
            category_combined = paste(unique(category), collapse="_")) %>% 
  # Add icustay information
  left_join(df_icustays_subject_selection %>% select(stay_id, intime, outtime, first_careunit, last_careunit, los), by=c("stay_id")) %>% 
  # Starttime and endtime as seq time (time since midnight of first day), so time of day info is preserved
  mutate(starttime_seq = as.numeric(difftime(starttime, min(as.Date(intime)), units="hours")),
         endtime_seq = as.numeric(difftime(endtime, min(as.Date(intime)), units="hours")),
         outtime_seq = as.numeric(difftime(outtime, min(as.Date(intime)), units="hours")),
         intime_seq = as.numeric(difftime(intime, min(as.Date(intime)), units="hours"))) %>% 
  # Add index to indicate if feeding periods are adjacent
  ungroup() %>% group_by(stay_id) %>% arrange(stay_id, starttime_seq) %>% 
  mutate(nut_intervalnr = c(0, cumsum(lead(starttime_seq) > cummax(endtime_seq))[-n()]))  

## Create groups of feeding types 
df_nut_all <- df_nut_all%>% group_by(stay_id) %>%
  mutate(nutritiongroup = case_when(all(category_combined  %in% c("Fluids/Intake", "Nutrition - Supplements")) ~ "Oral intake only",
                                    all(category_combined  %in% c("Nutrition - Enteral", "Nutrition - Supplements")) ~ "Enteral feeding only",
                                    all(category_combined  %in% c("Nutrition - Parenteral", "Nutrition - Supplements")) ~ "Parenteral feeding only",
                                    !(all(category_combined  %in% c("Nutrition - Enteral", "Nutrition - Supplements")) | all(category_combined  %in% c("Nutrition - Parenteral", "Nutrition - Supplements")) | all(ordercategoryname %in% c("Fluids/Intake", "Nutrition - Supplements"))) ~ "Combined"))

## Select patients that receive enteral nutrition only or that received enteral nutrition during their stay (i.e. before or after oral intake)
df_nut_enteral <- df_nut_all %>% filter(nutritiongroup %in% c("Combined", "Enteral feeding only")) %>% 
  arrange(stay_id, starttime_seq) %>% 
  group_by(stay_id, nut_intervalnr) %>% 
  filter(all(category_combined %in% c("Nutrition - Enteral","Nutrition - Supplements")))  

## Rename enteral nutrition
df_nut_enteral <- df_nut_enteral %>% rename(starttime_nut = starttime, endtime_nut=endtime, starttime_nut_seq = starttime_seq, endtime_nut_seq=endtime_seq, 
                                                   amount_nut = amount_combined, label_nut=label_combined)


## Create arrays of included stayids
incl_stayids_nut <- unique(df_nut_enteral$stay_id) 

## Add nutrition selection to inclusion flowchart dataframe 
df_flowchart_subjects$total_incl_icustays_entfeed_frompatientselection <- length(incl_stayids_nut)
df_flowchart_subjects$total_excl_icustays_noentfeed_frompatientselection <- nrow(df_icustays_subject_selection) - df_flowchart_subjects$total_incl_icustays_entfeed_frompatientselection

# WRAP-UP: 
# df_nut_all: dataframe with all feeding episodes of included icustays (incl_stayids_subject_selection)
# df_nut_enteral: dataframe with all enteral feeding episodes of included icustays (incl_stayids_subject_selection)
# incl_stayids_nut: array with stay_ids of included subjects (incl_stayids_subject_selection) who received enteral nutrition during their ICU stay


########################################################################################################################
##### STEP 4: Extract glucose data -----
## Glucose itemids from chartevents table (both lab and fingerstick is used:
#225664, -- Glucose finger stick
#220621, -- Glucose (serum)
#226537, -- Glucose (whole blood)

## Extract all glucose events from patients in df_icustay_selection
df_gluc_all <- tbl(mimic_icu, "chartevents") %>% 
  filter(itemid %in% c(225664,220621,226537) & stay_id %in% incl_stayids_subject_selection) %>% collect() %>% 
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
df_gluc_included_icustays <- df_gluc_all %>% filter(stay_id %in% incl_stayids_nut)

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
  left_join(df_nut_enteral %>% select(stay_id, starttime_nut_seq, endtime_nut_seq, amount_nut, label_nut,nut_intervalnr), 
            by=c("stay_id"), type = "many-to-many") %>% 
  filter(time_seq > starttime_nut_seq & time_seq <= endtime_nut_seq) %>% 
  select(stay_id, time_seq, starttime_nut_seq, endtime_nut_seq, amount_nut, label_nut, nut_intervalnr)

## Add enteral nutrition info to glucose data frame and only include glucose values taken during enteral feeding
df_gluc_nut <- df_gluc_valid %>% left_join(df_gluc_nutinfo, by=c("stay_id","time_seq")) %>% 
  arrange(stay_id, time_seq) %>% 
  # Remove duplicated measurements generated because feeding intervals occasionally overlap
  group_by(hadm_id, time_seq) %>% 
  distinct(valuenum, .keep_all = T) %>% 
  # Add enteral vs no enteral feeding label and rate of nutrition
  mutate(entfeeding_yn = ifelse(is.na(starttime_nut_seq), "n", "y"),
         rate_nut = amount_nut/(endtime_nut_seq - starttime_nut_seq)) %>% 
  # Only include glucose values during enteral feeding
  filter(entfeeding_yn == "y") %>% 
  # Add nutrition rate category
  mutate(rate_nut_cat = case_when(rate_nut > 50 ~ '>50 mL/h',
                                  rate_nut <= 50 ~ '=< 50 mL/h'))

## Extract ids of current selection
incl_stayid_gluc_selection <- unique(df_gluc_nut$stay_id)
incl_subjectid_gluc_selection <- unique(df_gluc_nut$subject_id)
incl_hadmid_gluc_selection <- unique(df_gluc_nut$hadm_id)

## Create dataframe including all data needed for analysis and rename columns
df_gluc_nut_variables <- df_gluc_nut %>% select(subject_id, hadm_id, stay_id, glucose_value=valuenum, glucose_unit=valueuom, charttime, storetime, time_seq, time_dec, time_hbin, 
                                                glucose_itemid=itemid, sample_type, intime, outtime, los, time_icudays, time_last_measurement, time_next_measurement,
                                                starttime_nut_seq, endtime_nut_seq, amount_nut, label_nut, nut_intervalnr, rate_nut, rate_nut_cat)

## Convert glucose value to mmol/L
df_gluc_nut_variables <- df_gluc_nut_variables %>% mutate(glucose_value_mmolL = glucose_value/18.0182)

## Update glucose measurement flowchart
df_flowchart_glucose_m$len_noENduring_glucose_m <- nrow(df_gluc_valid) - nrow(df_gluc_nut_variables) 
df_flowchart_glucose_m$len_yesENduring_glucose_m <- nrow(df_gluc_nut_variables)
df_flowchart_glucose_m$len_noENduring_glucose_m_patients <- length(incl_stayids_nut) - length(incl_stayid_gluc_selection)
df_flowchart_glucose_m$len_yesENduring_glucose_m_patients <- length(incl_stayid_gluc_selection)


## WRAP-UP: 
# df_gluc_nut_variables: time-stamped dataframe with glucose & nutrition data with only glucose measurements during EN
# incl_stayid_gluc_selection: array with stay_ids of final included subject (who have glucose measurements during enteral nutrition)
# incl_subjectid_gluc_selection: array with subject_ids of final included subject (who have glucose measurements during enteral nutrition)                                                                                               
# incl_hadmid_gluc_selection: array with hamd_ids of final included subject (who have glucose measurements during enteral nutrition)                                                                                                 


########################################################################################################################
##### STEP 6: Export flowchart dataframes to csv and array of included subject_ids/hadm_ids/stay_ids -----

write_csv(df_flowchart_subjects, file=paste0("2-output/mimic-iv_24h-glucose-variation_01_flowchart-data-subject-inclusion.csv"))
write_csv(df_flowchart_glucose_m, file=paste0("2-output/mimic-iv_24h-glucose-variation_01_flowchart-data-glucose-measurements-inclusion.csv"))

df_included_ids <- data.frame(incl_stayid_gluc_selection, incl_subjectid_gluc_selection, incl_hadmid_gluc_selection)
write_csv(df_flowchart_subjects, file=paste0("2-output/mimic-iv_24h-glucose-variation_01_included-subject-hadm-stay-ids.csv"))

## Count how many ICU days
icu_days_per_stay_id <- df_gluc_nut_variables %>%
  group_by(stay_id) %>%
  summarise(unique_count = n_distinct(time_icudays))
total_icu_days <- sum(icu_days_per_stay_id$unique_count)

########################################################################################################################
##### STEP 7: Extract and prepare relevant tables to add covariates to dataframe -----

## Extract patient info table
df_pat_info_gluc_selection <- tbl(mimic_hosp, "patients") %>% filter(subject_id %in% incl_subjectid_gluc_selection) %>% collect()

## Extract ICU stay info table
df_icustays_gluc_selection <- df_icustays_subject_selection  %>% filter(stay_id %in% incl_stayid_gluc_selection)

## Extract admission info table
df_adm_gluc_selection <- tbl(mimic_hosp, "admissions") %>% filter(hadm_id  %in% incl_hadmid_gluc_selection) %>% collect() %>% 
  select(subject_id, hadm_id, hosp_admittime=admittime, hosp_dischtime=dischtime, hosp_admission_type = admission_type, hosp_discharge_location = discharge_location)

## Combine patient, ICU stay and admission info tables 
df_general_patient_info <- df_pat_info_gluc_selection %>% 
  left_join(df_icustays_gluc_selection %>% select(subject_id, stay_id, intime, outtime, los), by="subject_id") %>% 
  left_join(df_adm_gluc_selection, by="subject_id")  

## Intime & outtime in hours since midnight of first day of ICU admission (faster downstream) 
df_general_patient_info <- df_general_patient_info %>% 
  mutate(intime_seq = as.numeric(difftime(intime, as.Date(intime), units="hours")),
         outtime_seq = as.numeric(difftime(outtime, as.Date(intime), units="hours")))


########################################################################################################################
##### STEP 8: Add covariates age & gender to dataframe -----
df_general_patient_info <- df_general_patient_info %>% 
  ## Recalculate age because age is dependent on ICU stay time 
  mutate(age = as.numeric(format(intime, "%Y")) - anchor_year + anchor_age) %>% 
  select(-c(anchor_age, anchor_year)) %>%
  ## Categorize age variable
  mutate(age_cat = cut(age, breaks=c(-Inf, 50, 70, Inf))) 

## Add age and gender to df_gluc_nut_variables
df_gluc_nut_variables <- df_gluc_nut_variables %>% 
  left_join(df_general_patient_info %>% select(subject_id, gender, age, age_cat), by="subject_id")


########################################################################################################################
##### STEP 9: Add covariate in-hospital mortality to dataframe -----
df_general_patient_info <- df_general_patient_info %>% 
  ## Use hospital discharge info to determine whether patient died during hospital stay or sent to hospice
  mutate(diff_death = difftime(dod, as.Date(hosp_dischtime), units="days"),
         died_hospstay_yn = case_when(as.numeric(diff_death) <= 0 | hosp_discharge_location %in% c("HOSPICE","DIED") ~ "y",
                                      TRUE ~ "n")) %>% 
  select(-diff_death)

## Add death-of-death and in-hospital mortality to df_gluc_nut_variables
df_gluc_nut_variables <- df_gluc_nut_variables %>% 
  left_join(df_general_patient_info %>% select(subject_id, dod, died_hospstay_yn), by="subject_id")


########################################################################################################################
##### STEP 10: Add covariate diabetes to dataframe -----
# Using icd codes from https://github.com/MIT-LCP/mimic-code/blob/main/mimic-iv/concepts/comorbidity/charlson.sql
dfdiag_diab <- tbl(mimic_hosp, "diagnoses_icd") %>% filter(hadm_id %in% incl_hadmid_gluc_selection) %>% 
  #icd codes from script linked above; ^ indicates that string should start with character string; | indicates OR
  filter(str_detect(icd_code, "^2500|^2501|^2502|^2503|^2508|^2509|^E100|^E10l|^E106|^E108|^E109|^E110|^E111|^E116|^E118|^E119|^E120|^E121|^E126|^E128|^E129|^E130|^E131|^E136|^E138|^E139|^E140|^E141|^E146|^E148|^E149|^2504|^2505|^2506|^2507|^E102|^E103|^E104|^E105|^E107|^E112|^E113|^E114|^E115|^E117|^E122|^E123|^E124|^E125|^E127|^E132|^E133|^E134|^E135|^E137|^E142|^E143|^E144|^E145|^E147")) %>% 
  collect()

df_gluc_nut_variables <- df_gluc_nut_variables %>% mutate(diabetes_status = ifelse(hadm_id %in% dfdiag_diab$hadm_id, "y", "n"))


########################################################################################################################
##### STEP 11: Add covariates dextrose administration to dataframe -----
items_dext <- dfItems %>% filter(grepl("dextrose", tolower(label)))  
itemid_dext <- items_dext$itemid

## Extract dextrose administration data
df_dext <- tbl(mimic_icu, "inputevents") %>% filter(itemid %in% itemid_dext) %>% 
  filter(stay_id %in% incl_stayid_gluc_selection) %>% filter(amount > 0) %>% collect() %>% 
  left_join(items_dext, by="itemid") %>% 
  left_join(df_general_patient_info %>% select(stay_id, intime), by=c("stay_id")) %>% 
  group_by(stay_id) %>% 
  mutate(starttime_dext_seq = as.numeric(difftime(starttime, as.Date(min(intime)), units="hours")),
         endtime_dext_seq = as.numeric(difftime(endtime, as.Date(min(intime)), units="hours"))) %>% 
  arrange(subject_id, starttime) %>% 
  select(stay_id, starttime_dext_seq, endtime_dext_seq, itemid_dext = itemid, label_dext=label, ordercat_dext = ordercategorydescription) 

## Link dextrose administration data to df_gluc_nut_variables dataframe
dfgluc_dextinfo <- df_gluc_nut_variables %>% ungroup() %>%  left_join(df_dext, by=c("stay_id"), type="many-to-many") %>% 
  filter(time_seq > starttime_dext_seq & time_seq <= endtime_dext_seq) %>% ungroup() %>% 
  group_by(stay_id, time_seq) %>% 
  summarize(starttime_dext_seq  = min(starttime_dext_seq),
            endtime_dext_seq = max(endtime_dext_seq))

## Add dextrose administration data to df_gluc_nut_variables dataframe
df_gluc_nut_variables <- df_gluc_nut_variables %>% 
  left_join(dfgluc_dextinfo %>% ungroup(), by=c("stay_id", "time_seq")) %>%
  mutate(dext_yn = ifelse(!is.na(starttime_dext_seq), "y","n"))  


########################################################################################################################
##### STEP 12: Add covariate insulin administration to dataframe -----
items_insulin <- dfItems %>% filter(grepl("insulin", tolower(label)) & category=="Medications") %>%
  select(itemid, label)
itemid_insulin <- items_insulin$itemid

## Extract insulin administration data
df_insulin <- tbl(mimic_icu, "inputevents") %>% filter(itemid %in% itemid_insulin) %>% 
  filter(stay_id %in% incl_stayid_gluc_selection) %>%   filter(amount > 0) %>% collect() %>% 
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
                              ordercat_insu %in% c("Drug Push") ~ "bolus"))

## Link insulin administration data to df_gluc_nut_variables dataframe
dfgluc_insuinfo <- df_gluc_nut_variables %>% ungroup() %>% left_join(df_insulin, by=c("stay_id"), type = "many-to-many") %>% 
  filter(time_seq > starttime_insu_seq_corr & time_seq <= endtime_insu_seq_corr) %>% ungroup() %>% 
  group_by(stay_id, time_seq) %>% 
  summarize(starttime_insu_seq_corr = min(starttime_insu_seq_corr),
            endtime_insu_seq_corr = max(endtime_insu_seq_corr)) %>%  
  ungroup() %>% 
  left_join(df_insulin, by=c("stay_id", "starttime_insu_seq_corr", "endtime_insu_seq_corr"), type="many-to-many") %>% 
  group_by(stay_id, time_seq, starttime_insu_seq_corr, endtime_insu_seq_corr) %>%
  summarize(insu_cat = unique(insu_cat))

## Add insulin administration data to df_gluc_nut_variables dataframe
df_gluc_nut_variables <- df_gluc_nut_variables %>% left_join(dfgluc_insuinfo %>% ungroup(), by=c("stay_id", "time_seq")) %>%
  mutate(insu_yn = ifelse(!is.na(starttime_insu_seq_corr), "y","n"))  

## Add insulin requirement per day (units per day in ICU; measure of insulin resistance) -----
df_insu_req <- df_insulin %>% mutate(time_icudays = starttime_insu_seq %/% 24)  %>% 
  group_by(stay_id, time_icudays) %>% 
  summarize(unitsperday_insu = sum(units_insu)) %>% ungroup()

df_gluc_nut_variables <- df_gluc_nut_variables %>% left_join(df_insu_req, by=c("stay_id", "time_icudays")) %>% 
  mutate(unitsperday_insu = ifelse(is.na(unitsperday_insu), 0, unitsperday_insu))


########################################################################################################################
##### STEP 13: Add covariate RASS scores to dataframe -----
## Extract RASS scores data
df_ras <- tbl(mimic_icu, "chartevents") %>% 
  filter(itemid %in% c(228096) & stay_id %in% incl_stayid_gluc_selection) %>% collect() %>% 
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
  left_join(dfras_int %>% select(stay_id, starttime_ras_seq, endtime_ras_seq, ras_score, ras_label), by=c("stay_id"), type="many-to-many") %>% 
  filter(time_seq > starttime_ras_seq & time_seq <= endtime_ras_seq) 

## Add RASS score data to df_gluc_nut_variables dataframe
df_gluc_nut_variables <- df_gluc_nut_variables %>% ungroup() %>% 
  left_join(dfgluc_rasinfo %>% select(stay_id, time_seq, ras_score, ras_label), 
            by=c("stay_id", "time_seq")) %>% 
  arrange(stay_id, time_seq)


########################################################################################################################
##### STEP 14: Add ventilation mode to dataframe -----
## Extract ventilation data
dfvent <- tbl(mimic_derived, "ventilation") %>% collect() %>% 
  #add intime to compute sequential time of day
  left_join(df_icustays_subject_selection %>% select(stay_id, intime), by=c("stay_id")) %>% 
  group_by(stay_id) %>% 
  mutate(starttime_vent_seq = as.numeric(difftime(starttime, as.Date(min(intime)), units="hours")),
         endtime_vent_seq = as.numeric(difftime(endtime, as.Date(min(intime)), units="hours")))

## Link ventilation data to df_gluc_nut_variables dataframe
dfgluc_ventinfo <- df_gluc_nut_variables %>% 
  left_join(dfvent, by=c("stay_id")) %>% 
  filter(time_seq > starttime_vent_seq & time_seq <= endtime_vent_seq) %>% ungroup()

## Add ventilation data to df_gluc_nut_variables dataframe
df_gluc_nut_variables <- df_gluc_nut_variables %>% ungroup() %>% 
  left_join(dfgluc_ventinfo %>% select(stay_id, time_seq, ventilation_status), by=c("stay_id", "time_seq"), type="many-to-many") %>% 
  arrange(stay_id, time_seq) %>% 
  # Reduce categories (from MIMIC code ventilation.sql); categories changed since MIMIC-IV 2.0
  mutate(ventilation_status_cat = case_when(ventilation_status %in% c("InvasiveVent", "Tracheostomy") ~ "invasive",
                                            ventilation_status %in% c("NonInvasiveVent", "HFNC", "SupplementalOxygen") ~ "noninvasive",
                                            ventilation_status %in% c("None") ~ "none", 
                                            is.na(ventilation_status) ~ "NaN")) 


########################################################################################################################
##### STEP 15: Add covariate corticosteroid administration to dataframe -----
# Corticosteroids include: dexamethasone, *cortisone, *prednisone, and *prednisolone (from: https://www.sciencedirect.com/science/article/pii/S0735675720300759)

cort_emar_tot <- tbl(mimic_hosp, "emar") %>% filter(hadm_id %in% incl_hadmid_gluc_selection) %>% 
  filter(grepl("dexamethasone", tolower(medication)) | grepl("cortisone", tolower(medication)) | grepl("prednisolone", tolower(medication)) | grepl("prednisone", tolower(medication))) %>% 
  collect() %>% 
  # Add icustay information
  left_join(df_general_patient_info %>% select(stay_id, hadm_id, intime, intime_seq, outtime_seq), by=c("hadm_id")) %>% 
  # Chartime to numeric time in hours since ICU admission
  group_by(hadm_id) %>% 
  # Use charttime if schedule time is missing 
  mutate(scheduletime = case_when(is.na(scheduletime) ~ charttime, 
                                  !is.na(scheduletime) ~ scheduletime)) %>% 
  #time as hours from midnight since start of icu stay (computationally faster than datetimes)
  mutate(charttime_seq = as.numeric(difftime(charttime, min(as.Date(intime)), units="hours")),
         scheduletime_seq = as.numeric(difftime(scheduletime, min(as.Date(intime)), units="hours"))) 

# Use pharmacy table to filter glucocorticoid creams and other non-oral non-iv routes of administration
cort_pharmid_tot <- unique(cort_emar_tot$pharmacy_id)

cort_pharm_tot <- tbl(mimic_hosp, "pharmacy") %>% filter(hadm_id %in% incl_hadmid_gluc_selection) %>%  
  filter(grepl("dexamethasone", tolower(medication)) | grepl("cortisone", tolower(medication)) | grepl("prednisolone", tolower(medication)) | grepl("prednisone", tolower(medication))) %>% 
  filter(pharmacy_id %in% cort_pharmid_tot) %>% collect() %>% 
  # Add icustay information
  left_join(df_general_patient_info %>% select(stay_id, hadm_id, intime), by=c("hadm_id")) %>% 
  group_by(hadm_id) %>% 
  mutate(starttime_seq = as.numeric(difftime(starttime, min(as.Date(intime)), units="hours")),
         stoptime_seq = as.numeric(difftime(stoptime, min(as.Date(intime)), units="hours")))

# Join emar and pharmacy tables
cort_emar_tot <- cort_emar_tot %>% 
  left_join(cort_pharm_tot %>% select(hadm_id, pharmacy_id, medication_pharm=medication, route, frequency), by=c("hadm_id", "pharmacy_id"))

# Filter final emar table
cort_emar <- cort_emar_tot %>% 
  #restrict data to administered doses (occurrence of other categories is very limited)
  filter(event_txt == "Administered" & 
           #restrict data to administrations that occurred during ICU stay
           charttime_seq >= intime_seq & charttime_seq <= outtime_seq) %>% 
  #exclude glucocorticoids administred as study medication or when pharmacy_id is missing
  filter(!grepl("study", tolower(medication)) & !grepl("placebo", tolower(medication)) & !is.na(pharmacy_id)) %>% 
  filter(route %in% c("IV","PO","ORAL","PO/NG")) %>% 
  arrange(stay_id, charttime_seq)

# Use scheduled administration times in EMAR to determine intervals during which glucocorticoids are active
dfemar_cort_interval <- cort_emar %>% mutate(scheduletime_seq = as.numeric(difftime(scheduletime, min(as.Date(intime)), units="hours")),
                                             scheduletime_seq = ifelse(is.na(scheduletime_seq), charttime_seq, scheduletime_seq),
                                             prev_dosetime_seq = scheduletime_seq - lag(scheduletime_seq, default=-Inf),
                                             cort_intervalno = cumsum(scheduletime_seq - lag(scheduletime_seq, default=-Inf) > 24)) %>% 
  # Use cort interval numbers to define start and end of periods of glucocorticoid administration; from first administration to 24h after last
  group_by(stay_id, cort_intervalno) %>% 
  summarize(cort_interval_starttime_seq = min(scheduletime_seq),
            cort_interval_endtime_seq = max(scheduletime_seq)+24,
            ndoses = n())  

dfgluc_cortinfo <- df_gluc_nut_variables %>% ungroup() %>%  left_join(dfemar_cort_interval, by=c("stay_id"), type="many-to-many") %>% 
  filter(time_seq > cort_interval_starttime_seq  & time_seq <= cort_interval_endtime_seq ) %>% ungroup() %>% 
  group_by(stay_id, time_seq) %>% 
  summarize(cort_interval_starttime_seq  = min(cort_interval_starttime_seq),
            cort_interval_endtime_seq = max(cort_interval_endtime_seq))

df_gluc_nut_variables <- df_gluc_nut_variables %>% left_join(dfgluc_cortinfo %>% ungroup(), by=c("stay_id", "time_seq")) %>%
  mutate(glucocorticoid_yn = ifelse(!is.na(cort_interval_starttime_seq), "y","n"))  

########################################################################################################################
##### STEP 16: Write output to csv file -----
write_csv(df_gluc_nut_variables, file=paste0("2-output/mimic-iv_24h-glucose-variation_01_glucose-incl-covariates.csv"))
write_csv(df_insulin, file=paste0("2-output/mimic-iv_24h-glucose-variation_01_insulinadministrations.csv"))
write_csv(df_dext, file=paste0("2-output/mimic-iv_24h-glucose-variation_01_dextroseadministrations.csv"))
write_csv(df_nut_enteral, file=paste0("2-output/mimic-iv_24h-glucose-variation_01_enteralfeeding.csv"))
write_csv(df_gluc_all, file=paste0("2-output/mimic-iv_24h-glucose-variation_01_all_gluc_measurements.csv"))
write_csv(dfgluc_cortinfo, file=paste0("2-output/mimic-iv_24h-glucose-variation_01_cortadministrations.csv"))

