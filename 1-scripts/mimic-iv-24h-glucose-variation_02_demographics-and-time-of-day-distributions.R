# Run-Time Environment: 	R version 4.3.1 
# Authors:				        Floor Hiemstra & Laura Kervezee
# Date created:           January 2024
# Project:      			    24-hour glucose variation in ICU patients
# Filename:               mimic-iv-24h-glucose-variation_02_demographics-and-time-of-day-distributions
# Purpose:                - Calculate demographics and characteristics of patients and glucose measurements
#                         - Plot time distributions of 
# Dependencies:           mimic-iv-24h-glucose-variation_01_preprocessing			    


########################################################################################################################
##### Packages -----
library(RPostgreSQL)
library(tidyverse)
library(gtsummary)
library(cowplot)
library(RColorBrewer)
library(lubridate)

########################################################################################################################
##### Load data from csv -----
df_gluc_nut_variables <- read_csv("2-output/mimic-iv_24h-glucose-variation_01_glucose-incl-covariates.csv", lazy=F, show_col_types=FALSE)
df_gluc_all <- read_csv("2-output/mimic-iv_24h-glucose-variation_01_all_gluc_measurements.csv", lazy=F, show_col_types=FALSE)
df_nut_enteral <- read_csv("2-output/mimic-iv_24h-glucose-variation_01_enteralfeeding.csv", lazy=F, show_col_types=FALSE)

## OR ##

########################################################################################################################
##### Load data from script -----
#source("C:\Users\fwhiemstra\Documents\mimic-glucose\1-scripts\R script versie 4\mimic-iv-24h-glucose-variation_01_preprocessing.R")

########################################################################################################################
##### Patient & glucose measurement characteristics -----
stayids_incl_subjects <- unique(df_gluc_nut_variables$stay_id)
hadmids_incl_subjects <- unique(df_gluc_nut_variables$hadm_id)
subjectids_incl_subjects <- unique(df_gluc_nut_variables$subject_id)

##### Patient characteristics -----
df_gluc_nut_variables %>% 
  group_by(stay_id, los, age, sex, diabetes_status, died_hospstay_yn, admtype_cat, race_cat, oasis, sofa) %>% 
  summarize() %>%
  ungroup() %>% select(-stay_id) %>% 
  tbl_summary()

## Count how many ICU days
icu_days_per_stay_id <- df_gluc_nut_variables %>%
  group_by(stay_id) %>%
  summarise(unique_count = n_distinct(time_icudays))
total_icu_days <- sum(icu_days_per_stay_id$unique_count)

##### Patient characteristics: Ventilation -----
# Summarize ventilation status during ICU stay as follows: 
# - invasive = patient had invasive ventilation at any moment during ICU stay
# - non-invasive = patient had non-invasive ventilation at any moment during ICU stay 
# - none = patient did not receive any invasive/non-invasive ventilation during ICU stay
mimic_derived <- dbConnect(RPostgres::Postgres(),
                           host     = "localhost",
                           dbname   = "mimic",
                           user     = "postgres",
                           password = "postgres",
                           bigint   = "integer",
                           port="5432",
                           options  = "-c search_path=mimiciv_derived")

df_ventilation_incl_patients <- tbl(mimic_derived, "ventilation") %>% collect() %>% filter(stay_id %in% stayids_incl_subjects)  %>% 
    # ventilation data is not available for all patients. Absence of ventilation data is considered as no ventilation. To include ICU_stays without ventilation
    # data to the dataframe, a right join is with a dataframe consisting of all stay_ids is done. 
    right_join(data.frame(stay_id = stayids_incl_subjects, extra_row_for_merge = stayids_incl_subjects), by="stay_id") %>% 
    mutate(ventilation_status_cat = case_when(ventilation_status %in% c("InvasiveVent", "Tracheostomy") ~ "invasive",
                                            ventilation_status %in% c("NonInvasiveVent", "HFNC", "SupplementalOxygen") ~ "noninvasive",
                                            ventilation_status %in% c("None") ~ "none", 
                                            is.na(ventilation_status) ~ "NaN")) %>% 
  group_by(stay_id) %>% 
  mutate(ventilation_status_stay = case_when(any(ventilation_status_cat == "invasive") ~ "invasive",
                                                any(ventilation_status_cat == "noninvasive") & all(ventilation_status_cat != "invasive") ~ "non-invasive only",
                                                all(ventilation_status_cat == "none") ~ "none", 
                                                is.na(ventilation_status) ~ "NaN")) %>%  
  ungroup() 

df_ventilation_incl_patients %>% group_by(stay_id) %>% 
  summarize(ventilation_status_stay = case_when(any(ventilation_status_cat == "invasive") ~ "invasive",
                                             any(ventilation_status_cat == "noninvasive") & all(ventilation_status_cat != "invasive") ~ "non-invasive only",
                                             all(ventilation_status_cat == "none") ~ "none",
                                             all(ventilation_status_cat == "NaN") ~ "NaN")) %>%  
  select(-stay_id) %>% tbl_summary()


##### Patient characteristics: Insulin administration -----
mimic_icu <- dbConnect(RPostgres::Postgres(),
                       host     = "localhost",
                       dbname   = "mimic",
                       user     = "postgres",
                       password = "postgres",
                       bigint   = "integer",
                       port="5432",
                       options  = "-c search_path=mimiciv_icu")
dfItems <- tbl(mimic_icu, "d_items") %>% collect() #table with definition of itemids (e.g. types of nutrition)
items_insulin <- dfItems %>% filter(grepl("insulin", tolower(label)) & category=="Medications") %>%
  select(itemid, label)
itemid_insulin <- items_insulin$itemid
df_insulin <- tbl(mimic_icu, "inputevents") %>% filter(itemid %in% itemid_insulin) %>% 
  filter(stay_id %in% stayids_incl_subjects) %>% filter(amount > 0) %>% collect() %>% 
  left_join(items_insulin, by="itemid")
## df_insulin = dataframe from patients with any insulin administration during their ICU stay. 
length(unique(df_insulin$stay_id)) # number of patients/ICU stays with insulin at any moment during their ICU stay

##### Patient characteristics: Dextrose administration -----
items_dext <- dfItems %>% filter(grepl("dextrose", tolower(label)))  
itemid_dext <- items_dext$itemid
df_dext <- tbl(mimic_icu, "inputevents") %>% filter(itemid %in% itemid_dext) %>% 
  filter(stay_id %in% stayids_incl_subjects) %>% filter(amount > 0) %>% collect() %>% 
  left_join(items_dext, by="itemid")
## df_dext = dataframe from patients with any dextrose administration during their ICU stay. 
length(unique(df_dext$stay_id)) # number of patients/ICU stays with dextrose at any moment during their ICU stay


##### Patient characteristics: Glucocorticoids administration -----
df_icustays_selection <- tbl(mimic_icu, "icustays") %>% filter(stay_id %in% stayids_incl_subjects) %>% collect() %>% 
  select(subject_id, stay_id, hadm_id, intime, outtime, los) %>% 
  mutate(intime_seq = as.numeric(difftime(intime, as.Date(intime), units="hours")),
         outtime_seq = as.numeric(difftime(outtime, as.Date(intime), units="hours")))

mimic_hosp <- dbConnect(RPostgres::Postgres(),
                        host     = "localhost",
                        dbname   = "mimic",
                        user     = "postgres",
                        password = "postgres",
                        bigint   = "integer",
                        port="5432",
                        options  = "-c search_path=mimiciv_hosp")

cort_emar_tot <- tbl(mimic_hosp, "emar") %>% filter(hadm_id %in% hadmids_incl_subjects) %>% 
  filter(grepl("dexamethasone", tolower(medication)) | grepl("cortisone", tolower(medication)) | grepl("prednisolone", tolower(medication)) | grepl("prednisone", tolower(medication))) %>% 
  collect() %>% 
  # Add icustay information
  left_join(df_icustays_selection %>% select(stay_id, hadm_id, intime, intime_seq, outtime_seq), by=c("hadm_id")) %>% 
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
cort_pharm_tot <- tbl(mimic_hosp, "pharmacy") %>% filter(hadm_id %in% hadmids_incl_subjects) %>%  
  filter(grepl("dexamethasone", tolower(medication)) | grepl("cortisone", tolower(medication)) | grepl("prednisolone", tolower(medication)) | grepl("prednisone", tolower(medication))) %>% 
  filter(pharmacy_id %in% cort_pharmid_tot) %>% collect() %>% 
  # Add icustay information
  left_join(df_icustays_selection %>% select(stay_id, hadm_id, intime), by=c("hadm_id")) %>% 
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
## cort_emar = dataframe from patients with any glucocorticoids administration during their ICU stay. 
length(unique(cort_emar$stay_id)) # number of patients/ICU stays with glucocorticoids at any moment during their ICU stay


##### Glucose measurements characteristics - N measurements per patient per day -----
median_glucose_per_patient_per_day <- df_gluc_nut_variables %>% group_by(stay_id, time_icudays) %>%
  summarize(nglucperday = length(glucose_value)) %>% ungroup() %>% 
  summarize(nglucperday_mean = median(nglucperday),
            nglucperday_IQR_low = quantile(nglucperday, 0.25),
            nglucperday_IQR_upp = quantile(nglucperday, 0.75))
median_glucose_per_patient_per_day

##### Glucose measurements characteristics - Medication administration and sample type -----
df_gluc_nut_variables %>% select(ventilation_status_cat, insu_yn, dext_yn, glucocorticoid_yn, sample_type) %>% 
  tbl_summary()



########################################################################################################################
##### Time-of-day distributions -----

theme_plot <- theme(panel.grid = element_blank(), 
                    
                    # Legend settings 
                    legend.position="bottom", legend.justification = c("center", "top"), 
                    legend.text = element_text(colour="black", size=10, hjust=0.5), 
                    legend.key.height = unit(6, units="pt"), 
                    legend.title = element_text(colour="black", size=10, face="bold"),
                    legend.background = element_blank(),
                    
                    # Axes  
                    axis.text=element_text(colour="black", size=8),
                    axis.ticks = element_line(size = 0.25, color="black"),
                    
                    # Axis title 
                    axis.title = element_text(colour="black", size=10), 
                    plot.background=element_rect(fill="transparent", color=NA),
                    
                    # Title 
                    plot.title=element_text(colour="black", size=10, hjust=0.5, face="bold"), 
                    
                    strip.text=element_text(face="bold"), 
                    panel.border=element_rect(fill=NA, colour="black", size=0.25))

##### Glucose measurements ----- 
color_gluc_plot <- "#882255" 

plt_gluc_time_distribution <- ggplot(df_gluc_nut_variables, aes(x=time_hbin)) + 
  geom_bar(color=color_gluc_plot, width=0.8, size=0.2, fill=color_gluc_plot)+
  xlab("Time of day (h)")+
  ylab("Number of measurements") +
  scale_x_continuous(expand=c(0.01,0), breaks=seq(0, 24, by = 4))+
  scale_y_continuous(expand=expansion(mult = c(0, 0.04))) +
  theme_bw() + # theme_bw will remove the background
  theme_plot + 
  ggtitle("Glucose measurements") +
  guides(fill = guide_legend(title.position="top", title.hjust=0.5, label.hjust=0.5))+coord_cartesian(clip = "off")
plt_gluc_time_distribution

##### Enteral feeding -----
df_nut_enteralrate <- df_nut_enteral %>% mutate(difftime = endtime_nut_seq - starttime_nut_seq,
                                                rate_nut = amount_nut / (endtime_nut_seq - starttime_nut_seq),
                                                starttime_dec = starttime_nut_seq %% 24,
                                                endtime_dec = endtime_nut_seq %% 24) %>%
  filter(stay_id %in% stayids_incl_subjects)

df_nut_enteralrate <- df_nut_enteralrate %>% 
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


dfhour_rate <- df_nut_enteral %>% ungroup() %>% #filter(stay_id %in% sample_pat) %>% 
  select(stay_id, outtime_seq, intime_seq) %>% 
  ungroup() %>% distinct(stay_id, outtime_seq, intime_seq) %>% 
  group_by(stay_id) %>% 
  mutate(time_seq = floor(intime_seq)) %>% 
  #create dataframe with a row for every hour during ICU stay
  complete(time_seq = seq(from=floor(intime_seq), to=ceiling(outtime_seq), by=1)) %>% 
  mutate(outtime_seq = ceiling(max(outtime_seq, na.rm=T)),
         intime_seq = floor(max(intime_seq, na.rm=T ))) %>% 
  #add nutrition intervals to indicate whther there was feeding for each hour during ICU stay
  left_join(df_nut_enteralrate %>% ungroup() %>% select(stay_id, starttime_nut_seq, endtime_nut_seq, rate_nut_cho), by="stay_id", relationship="many-to-many") %>% 
  arrange(stay_id, time_seq) %>%  
  filter(time_seq >= ceiling(starttime_nut_seq) & time_seq < ceiling(endtime_nut_seq)) %>% #this line is different from dfhour_feeding to account for consecutive periods of feeding
  filter(!is.na(rate_nut_cho)) # Remove rows with NaN values in rate_nut_cho column


dfhour_ratesumm <- dfhour_rate %>% mutate(time_dec = time_seq %% 24,
                                          rate_nut_cho_cat = cut(rate_nut_cho, breaks= c(-Inf, 4.5, 6.5, 8.5, Inf))) %>% 
  group_by(time_dec, rate_nut_cho_cat) %>% 
  summarize(n = n()) %>% 
  mutate(prop = n/sum(n)) %>% 
  mutate(rate_nut_cho_cat = factor(rate_nut_cho_cat, levels=c("(-Inf,4.5]", "(4.5,6.5]", "(6.5,8.5]", "(8.5, Inf]"),
                                   labels=c("\u2264 4.5 grams/hour", "4.5-6.5 grams/hour", "6.5-8.5 grams/hour", "\u2265 8.5 grams/hour" )))

plt_entrate_timeofday <- ggplot(dfhour_ratesumm, aes(x=time_dec, y=prop, fill=rate_nut_cho_cat))+ theme_bw()+
  geom_bar(stat="identity", color="#6BAED6", position = position_stack(reverse = T), size=0.2)+
  scale_x_continuous(breaks=seq(0,24, by=4), expand=c(0,0), limits=c(-0.5,23.5))+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values = (c(brewer.pal(7, "Blues")[2:7])), name="Carbohydrate administration rate") +
  xlab("Time of day (h)")+ylab("Proportion of glucose measurements")+
  coord_cartesian(clip="off")+
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5, ncol = 2)) +  # Set ncol to 2 for two lines
  theme_plot ##+ 
  ##ggtitle("")
plt_entrate_timeofday


##### Insulin administration in included dataset -----

color_insu_plot <- "#117733"

p_insu_per_hour = df_gluc_nut_variables %>% group_by(time_hbin) %>% 
  summarise(insu_yes = sum(insu_yn == "y"),
            insu_no = sum(insu_yn == "n")) %>% 
  mutate(prop_insu = insu_yes / (insu_yes + insu_no))

plt_insu_timeofday <- ggplot(p_insu_per_hour, aes(x=time_hbin, y=prop_insu))+ theme_bw()+
  geom_bar(stat="identity", color=color_insu_plot, width=0.8, size=0.2, fill=color_insu_plot) +
  xlab("Time of day (h)")+
  ylab("Proportion of glucose measurements") +
  scale_x_continuous(expand=c(0.01,0), breaks=seq(0, 24, by = 4))+
  scale_y_continuous(expand=c(0.01,0), breaks=seq(0, 1, by = 0.2), limits=c(0,1)) +
  theme_bw() + # theme_bw will remove the background
  theme_plot + 
  ggtitle("Insulin administration") +
  guides(fill = guide_legend(title.position="top", title.hjust=0.5, label.hjust=0.5))+coord_cartesian(clip = "off")
plt_insu_timeofday

##### Dextrose administration in included dataset -----

color_dextrose_plot <- "#882255" 

p_dextrose_per_hour = df_gluc_nut_variables %>% group_by(time_hbin) %>% 
  summarise(dext_yes = sum(dext_yn == "y"),
            dext_no = sum(dext_yn == "n")) %>% 
  mutate(prop_dext = dext_yes / (dext_yes + dext_no))

plt_dext_timeofday <- ggplot(p_dextrose_per_hour, aes(x=time_hbin, y=prop_dext))+ theme_bw()+
  geom_bar(stat="identity", color=color_dextrose_plot, width=0.8, size=0.2, fill=color_dextrose_plot) +
  xlab("Time of day (h)")+
  ylab("Proportion of glucose measurements") +
  scale_x_continuous(expand=c(0.01,0), breaks=seq(0, 24, by = 4))+
  scale_y_continuous(expand=c(0.01,0), breaks=seq(0, 1, by = 0.2), limits=c(0,1)) +
  theme_bw() + # theme_bw will remove the background
  theme_plot + 
  ggtitle("Dextrose administration") +
  guides(fill = guide_legend(title.position="top", title.hjust=0.5, label.hjust=0.5))+coord_cartesian(clip = "off")
plt_dext_timeofday

##### Glucocorticoids administration in included dataset -----

color_glucocorticoid_plot <- "#88CCEE" #"#BE398D"

p_glucocorticoid_per_hour = df_gluc_nut_variables %>% group_by(time_hbin) %>% 
  summarise(glucocorticoid_yes = sum(glucocorticoid_yn == "y"),
            glucocorticoid_no = sum(glucocorticoid_yn == "n")) %>% 
  mutate(prop_glucocorticoid = glucocorticoid_yes / (glucocorticoid_yes + glucocorticoid_no))

plt_glucocorticoid_timeofday <- ggplot(p_glucocorticoid_per_hour, aes(x=time_hbin, y=prop_glucocorticoid))+ theme_bw()+
  geom_bar(stat="identity", color=color_glucocorticoid_plot, width=0.8, size=0.2, fill=color_glucocorticoid_plot) +
  xlab("Time of day (h)")+
  ylab("Proportion of glucose measurements") +
  scale_x_continuous(expand=c(0.01,0), breaks=seq(0, 24, by = 4))+
  scale_y_continuous(expand=c(0.01,0), breaks=seq(0, 1, by = 0.2), limits=c(0,1)) +
  theme_bw() + # theme_bw will remove the background
  theme_plot + 
  ggtitle("Glucocorticoid administration") +
  guides(fill = guide_legend(title.position="top", title.hjust=0.5, label.hjust=0.5))+coord_cartesian(clip = "off")
plt_glucocorticoid_timeofday




##### Combine plots 
plts_time_distribution_gluc_feeding <- cowplot::plot_grid(plt_gluc_time_distribution, plt_entrate_timeofday, 
                                                          ncol=2, nrow=1, align="hv", axis="tblr", labels=c("a","b"), label_size = 10)
plts_time_distribution_gluc_feeding

plts_time_distribution_administrations <- cowplot::plot_grid(plt_insu_timeofday, plt_dext_timeofday, plt_glucocorticoid_timeofday,
                                                             ncol=3, nrow=1, align="hv", axis="tblr", labels=c("a","b", "c"), label_size = 10)
plts_time_distribution_administrations


##### Save plots -----
ggsave(plot = plt_entrate_timeofday, width = 4, height = 4, dpi = 500, filename = "2-output/mimic-iv_24h-glucose_variation_02_time-of-day-distribution-enteral-feeding-rates.png")
ggsave(plot = plt_gluc_time_distribution, width = 4, height = 4, dpi = 500, filename = "2-output/mimic-iv_24h-glucose_variation_02_time-of-day-distribution-glucose.png")

ggsave(plot = plt_entrate_timeofday, width = 4, height = 4, dpi = 500, filename = "2-output/mimic-iv_24h-glucose_variation_02_time-of-day-distribution-enteral-feeding-rates.svg")
ggsave(plot = plt_gluc_time_distribution, width = 4, height = 4, dpi = 500, filename = "2-output/mimic-iv_24h-glucose_variation_02_time-of-day-distribution-glucose.svg")

ggsave(plot = plts_time_distribution_administrations, width = 8, height = 3, dpi = 500, filename = "2-output/mimic-iv_24h-glucose_variation_02_time-of-day-distribution-administrations-combined.png")
ggsave(plot = plts_time_distribution_administrations, width = 8, height = 3, dpi = 500, filename = "2-output/mimic-iv_24h-glucose_variation_02_time-of-day-distribution-administrations-combined.svg")




########################################################################################################################
##### Example of timing of enteral feeding, glucose measurements and medication administration -----
#SELECTION <- df_gluc_nut_variables %>% filter(dext_yn == "y" & glucocorticoid_yn == "y" & insu_yn == "y") 
example_stayid <- 30045625     
dfexall <- df_gluc_nut_variables %>% filter(stay_id == example_stayid) %>% filter(time_seq > 72 & time_seq < 144)
dfexall_nut <- dfexall %>% distinct(stay_id, starttime_nut_seq , endtime_nut_seq, intime) %>% drop_na
dfexall_insu <- dfexall %>% distinct(stay_id, starttime_insu_seq_corr, endtime_insu_seq_corr, intime) %>% drop_na
dfexall_dext <- dfexall %>% distinct(stay_id, starttime_dext_seq , endtime_dext_seq, intime) %>% drop_na
dfexall_gluccort <- dfexall %>% distinct(stay_id, cort_interval_starttime_seq , cort_interval_endtime_seq, intime) %>% drop_na
dfexall_gluc_all <- df_gluc_all %>% filter(stay_id == example_stayid)%>% 
  mutate(time_seq = as.numeric(difftime(charttime, as.Date(min(intime)), units="hours")))%>% filter(time_seq > 72 & time_seq < 144)
  
plt_example_patient <- ggplot(dfexall , aes(x=time_seq, y="glucose"))+theme_bw()+
  geom_vline(xintercept=seq(48, 168, by=24), linetype="dotted", color="grey80", size=0.2)+
  geom_point(data=dfexall_gluc_all, aes(x=time_seq, y="glucose"), shape="x", color="grey70", size=1.1) + 
  geom_point(shape="x", color="black", size=1.1)+
  geom_segment(data=dfexall_nut , aes(x=starttime_nut_seq, xend=endtime_nut_seq), y="feeding", yend="feeding", color="#117733", size=2)+
  geom_segment(data=dfexall_insu , aes(x=starttime_insu_seq_corr, xend=endtime_insu_seq_corr), y="insulin", yend="insulin", color="#192e27", size=2)+
  geom_segment(data=dfexall_dext , aes(x=starttime_dext_seq, xend=endtime_dext_seq), y="dextrose", yend="dextrose", color="#9f93d1", size=2)+
  geom_segment(data=dfexall_gluccort , aes(x=cort_interval_starttime_seq, xend=cort_interval_endtime_seq), y="glucocorticoids", yend="glucocorticoids", color="#88CCEE", size=2)+
  scale_y_discrete(limits=c("day","glucocorticoids", "insulin", "dextrose","feeding", "glucose"), labels=c("", "Glucocorticoids", "Insulin", "Dextrose","Enteral\nfeeding", "Glucose\nmeasurement"), expand=expansion(mult=c(0, 0.1)))+
  scale_x_continuous(breaks=seq(0, 240, by=24))+
  coord_cartesian(xlim=c(72, 144), clip="on")+
  xlab("")+ylab("")+
  annotate(geom="text", x=seq(72, 168, by=24)+12, y="day", label=paste("Day", 4:8), vjust=1.4, size=2, fontface="bold")+
  theme_plot
plt_example_patient

ggsave(plt_example_patient, filename=paste0("2-output/mimic-iv_glucose-project_02_time-of-day-example_pat30045625.png"), dpi=600, width=6, height=2)
ggsave(plt_example_patient, filename=paste0("2-output/mimic-iv_glucose-project_02_time-of-day-example_pat30045625.SVG"), dpi=600, width=6, height=2)


########################################################################################################################
##### Example of ICU feeding pattern -----

samplepat <- c(32795641, 38870879)

df_nut_enteral_intervals <- df_nut_enteral %>% filter(stay_id %in% samplepat) %>% 
  group_by(stay_id, nut_intervalnr) %>% 
  summarize(nut_interval_start = min(starttime_nut_seq ),
            nut_interval_end = max(endtime_nut_seq ))

dfhour_feeding <- df_nut_enteral %>% ungroup() %>% filter(stay_id %in% samplepat) %>% 
  ungroup() %>% distinct(stay_id, outtime_seq, intime_seq) %>% 
  group_by(stay_id) %>% 
  mutate(time_seq = floor(intime_seq)) %>% 
  #create dataframe with a row for every hour during ICU stay
  complete(time_seq = seq(from=floor(intime_seq), to=ceiling(outtime_seq), by=1)) %>% 
  mutate(outtime_seq = ceiling(max(outtime_seq, na.rm=T)),
         intime_seq = floor(max(intime_seq, na.rm=T ))) %>% 
  #add nutrition intervals to indicate whther there was feeding for each hour during ICU stay
  left_join(df_nut_enteral_intervals %>% ungroup() %>% select(stay_id, nut_interval_start, nut_interval_end), by="stay_id", relationship="many-to-many") %>% 
  arrange(stay_id, time_seq) %>%  
  filter(time_seq >= floor(nut_interval_start) & time_seq < ceiling(nut_interval_end)) %>% 
  group_by(stay_id) %>% 
  complete(time_seq = seq(from=min(intime_seq,na.rm=T), max(outtime_seq, na.rm=T), by=1)) %>% 
  mutate(feeding_yn = ifelse(!is.na(nut_interval_start), "y","n"),
         day = time_seq %/% 24,
         time_dec = time_seq %% 24,
         patnumplot = ifelse(stay_id == c(32795641), "#1", "#2"))

pl_feedingacto <- dfhour_feeding %>% 
  group_by(stay_id, day, feeding_yn) %>% 
  summarize(time_dec_start = min(time_dec), time_dec_end=max(time_dec+1), .groups = "drop") %>% 
  arrange(stay_id, day, time_dec_start) %>% 
  ggplot(aes(xmin=time_dec_start,xmax=time_dec_end, ymin=day-0.5, ymax=day+0.5))+theme_bw()+
  ggforce::facet_col(vars(paste("Patient", stay_id)), scales = "free", space = "free")+
  geom_rect(color="grey90", aes(fill=feeding_yn))+
  scale_y_continuous(breaks=0:50, expand=c(0,0), trans="reverse")+
  scale_x_continuous(breaks=seq(0, 24, by=4), expand=expansion(add=0.01,0.01), limits=c(-0,24))+
  scale_fill_manual(values=c("n"="grey90", "y"="#440154"))+
  xlab("Time of day (h)")+ylab("Day in ICU")+
  coord_cartesian(clip="off")+
  theme(legend.position="none", axis.text=element_text(colour="black", size=7), panel.grid = element_blank(),
        axis.title = element_text(colour="black", size=7), legend.text = element_text(colour="black", size=7), #plot.background=element_rect(fill="transparent", color=NA),
        plot.title=element_text(colour="black", size=7, hjust=0.5, face="bold"), legend.key.width = unit(7, units="pt"), legend.key.height = unit(6, units="pt"), legend.title = element_text(colour="black", size=7),
        strip.background = element_blank(), strip.text=element_text(face="bold"), legend.background = element_blank(), panel.border=element_rect(fill=NA, colour="black", size=0.25),  axis.ticks = element_line(size = 0.25, color="black"))

ggsave(pl_feedingacto, filename = paste0("2-output/mimic-iv_glucose-variation_02_feeding-actograms.pdf"), width=3, height = 5 )




















































