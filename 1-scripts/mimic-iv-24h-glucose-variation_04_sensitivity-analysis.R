# Run-Time Environment: 	R version 4.3.1 & R Studio 1.4.1106
# Authors:				        Floor Hiemstra & Laura Kervezee
# Date created:           January 2024
# Project:      			    24-hour glucose variation in ICU patients
# Filename:               mimic-iv-24h-glucose-variation_04_sensitivity_analysis
# Purpose:  			        Sensitivity analysis of subgroups based on: 
#                         - Ventilation mode
#                         - Survivor status 
#                         - Sedation depth
#                         - Days in the ICU
#                         - Sample type 
#                         - Ventilation mode
#                         - Sample frequency
#                         - Insulin resistance
#                         Sensitivity analysis is performed by fitting linear mixed-effects models to each subgroup
#                         Visualization is done by plotting estimated marignal means for each fitted model.

########################################################################################################################
##### Packages -----
library(tidyverse)
library(stringr)
library(lubridate)
library(lmerTest)
library(emmeans)
library(sjPlot)
library(car)
library(cowplot)


########################################################################################################################
##### Load data from csv -----
df_gluc_nut_variables <- read_csv("2-output/mimic-iv_24h-glucose-variation_01_glucose-incl-covariates.csv", lazy=F, show_col_types=FALSE)

## OR ##

########################################################################################################################
##### Load data from script -----
#source("C:\Users\fwhiemstra\Documents\mimic-glucose\1-scripts\R script versie 4\mimic-iv-24h-glucose-variation_01_preprocessing.R")


## Select and prepare variables for LME model and sensitivity analysis 
df_gluc_data_analysis <- df_gluc_nut_variables %>% 
  mutate(time_hbin_factor = as.factor(time_hbin),
         age_cat = as.factor(age_cat), 
         sex = as.factor(sex),
         diabetes_status = as.factor(diabetes_status), 
         rate_insu_cat = as.factor(rate_insu_cat),  
         rate_dext_cat = as.factor(rate_dext_cat), 
         glucocorticoid_yn = as.factor(glucocorticoid_yn), 
         rate_nut_cho_cat = as.factor(rate_nut_cho_cat)
  ) %>% 
  select(glucose_value_mmolL, 
         time_hbin, time_hbin_factor, 
         stay_id, 
         age_cat,  
         sex, 
         diabetes_status, 
         rate_insu_cat, 
         rate_dext_cat, 
         glucocorticoid_yn, 
         rate_nut_cho_cat,
         died_hospstay_yn, 
         ventilation_status_cat, 
         time_icudays, 
         time_next_measurement, 
         sample_type,
         ras_score,
         unitsperday_insu) 

# Set reference categories
df_gluc_data_analysis$diabetes_status = relevel(df_gluc_data_analysis$diabetes_status, ref = "n")
df_gluc_data_analysis$sex = relevel(df_gluc_data_analysis$sex, ref = "F")
df_gluc_data_analysis$glucocorticoid_yn = relevel(df_gluc_data_analysis$glucocorticoid_yn, ref = "n")
df_gluc_data_analysis$age_cat = relevel(df_gluc_data_analysis$age_cat, ref = "[-Inf,55]")
df_gluc_data_analysis$rate_insu_cat = relevel(df_gluc_data_analysis$rate_insu_cat, ref = "[-Inf,0]")
df_gluc_data_analysis$rate_dext_cat = relevel(df_gluc_data_analysis$rate_dext_cat, ref = "[-Inf,0]")
df_gluc_data_analysis$rate_nut_cho_cat = relevel(df_gluc_data_analysis$rate_nut_cho_cat, ref = "[-Inf,4.5]")


## Define variables 
dependent_variable <- "glucose_value_mmolL"
time_variable <- "time_hbin_factor"
patient_level_variables <- c("diabetes_status", "age_cat", "sex")
sample_level_variables <- c("rate_insu_cat", "rate_dext_cat", "glucocorticoid_yn", "rate_nut_cho_cat")
random_variable <- "(1 | stay_id)"

## LME model formula
formula_model_final <-  as.formula(paste(dependent_variable, "~", time_variable, " + ", paste(patient_level_variables, collapse=" + "), " + ", paste(sample_level_variables, collapse=" + "), " + ", random_variable))

  
###############################################################################
##### Sensitivity analyses - Ventilation mode -----

df_sensanalysis_group <- df_gluc_data_analysis %>% filter(ventilation_status_cat == "invasive")
lme_sensanalysis_group <- lmer(formula_model_final, data=df_sensanalysis_group, REML=TRUE)
emm_sensanalysis_ventilation_invasive <- data.frame(emmeans(lme_sensanalysis_group, ~time_hbin_factor, rg.limit = 60000, weights="prop"), 
                                                      group = "Invasive",
                                                      npat = length(unique(df_sensanalysis_group$stay_id)),
                                                      nmeas = nrow(df_sensanalysis_group))

df_sensanalysis_group <- df_gluc_data_analysis %>% filter(ventilation_status_cat == "noninvasive")
lme_sensanalysis_group <- lmer(formula_model_final, data=df_sensanalysis_group, REML=TRUE)
emm_sensanalysis_ventilation_noninvasive <- data.frame(emmeans(lme_sensanalysis_group, ~time_hbin_factor, rg.limit = 60000, weights="prop"), 
                                                    group = "Non-Invasive",
                                                    npat = length(unique(df_sensanalysis_group$stay_id)),
                                                    nmeas = nrow(df_sensanalysis_group))
emm_sensanalysis_ventilation <- bind_rows(emm_sensanalysis_ventilation_invasive, emm_sensanalysis_ventilation_noninvasive)

plt_text_ventilation_mode <- emm_sensanalysis_ventilation %>% group_by(group, npat, nmeas) %>% 
  summarize(mean = mean(emmean), .groups="drop") %>% 
  mutate(y = case_when(mean == max(mean) ~ Inf, 
                       mean == min(mean) ~ -Inf,
                       TRUE ~ -Inf),
         vjust = case_when(mean == max(mean) ~ 1.2,
                           mean == min(mean) ~ -0.2,
                           TRUE ~ -1.4))

###############################################################################
##### Sensitivity analyses - Survival status -----
df_sensanalysis_group <- df_gluc_data_analysis %>% filter(died_hospstay_yn == "y")
lme_sensanalysis_group <- lmer(formula_model_final, data=df_sensanalysis_group, REML=TRUE)
emm_sensanalysis_survival_died <- data.frame(emmeans(lme_sensanalysis_group, ~time_hbin_factor, rg.limit = 60000, weights="prop"), 
                                                    group = "Non-survivors",
                                                    npat = length(unique(df_sensanalysis_group$stay_id)),
                                                    nmeas = nrow(df_sensanalysis_group))

df_sensanalysis_group <- df_gluc_data_analysis %>% filter(died_hospstay_yn == "n")
lme_sensanalysis_group <- lmer(formula_model_final, data=df_sensanalysis_group, REML=TRUE)
emm_sensanalysis_survival_survived <- data.frame(emmeans(lme_sensanalysis_group, ~time_hbin_factor, rg.limit = 60000, weights="prop"), 
                                                       group = "Survivors",
                                                       npat = length(unique(df_sensanalysis_group$stay_id)),
                                                       nmeas = nrow(df_sensanalysis_group))
emm_sensanalysis_survival <- bind_rows(emm_sensanalysis_survival_died, emm_sensanalysis_survival_survived)

plt_text_survival <- emm_sensanalysis_survival %>% group_by(group, npat, nmeas) %>% 
  summarize(mean = mean(emmean), .groups="drop") %>% 
  mutate(y = case_when(mean == max(mean) ~ Inf, 
                       mean == min(mean) ~ -Inf,
                       TRUE ~ -Inf),
         vjust = case_when(mean == max(mean) ~ 1.2,
                           mean == min(mean) ~ -0.2,
                           TRUE ~ -1.4))


###############################################################################
##### Sensitivity analyses - Sedation depth -----

df_sensanalysis_group <- df_gluc_data_analysis %>% filter(ras_score > -2) 
lme_sensanalysis_group <- lmer(formula_model_final, data=df_sensanalysis_group, REML=TRUE)
emm_sensanalysis_rass_high <- data.frame(emmeans(lme_sensanalysis_group, ~time_hbin_factor, rg.limit = 60000, weights="prop"), 
                                             group = "RASS ≥ -1",
                                             npat = length(unique(df_sensanalysis_group$stay_id)),
                                             nmeas = nrow(df_sensanalysis_group))

df_sensanalysis_group <- df_gluc_data_analysis %>% filter(ras_score <= -2)
lme_sensanalysis_group <- lmer(formula_model_final, data=df_sensanalysis_group, REML=TRUE)
emm_sensanalysis_rass_low <- data.frame(emmeans(lme_sensanalysis_group, ~time_hbin_factor, rg.limit = 60000, weights="prop"), 
                                                 group = "RASS ≤ -2",
                                                 npat = length(unique(df_sensanalysis_group$stay_id)),
                                                 nmeas = nrow(df_sensanalysis_group))
emm_sensanalysis_sedation <- bind_rows(emm_sensanalysis_rass_high, emm_sensanalysis_rass_low)


emm_sensanalysis_sedation$group <- factor(emm_sensanalysis_sedation$group, levels = c(
  "RASS ≥ -1",
  "RASS ≤ -2"
))


plt_text_sedation <- emm_sensanalysis_sedation %>% group_by(group, npat, nmeas) %>% 
  summarize(mean = mean(emmean), .groups="drop") %>% 
  mutate(y = case_when(mean == max(mean) ~ Inf, 
                       mean == min(mean) ~ -Inf,
                       TRUE ~ -Inf),
         vjust = case_when(mean == max(mean) ~ 1.2,
                           mean == min(mean) ~ -0.2,
                           TRUE ~ -1.4))

###############################################################################
##### Sensitivity analyses - Sample type -----

df_sensanalysis_group <- df_gluc_data_analysis %>% filter(sample_type == "Fingerstick")  
lme_sensanalysis_group <- lmer(formula_model_final, data=df_sensanalysis_group, REML=TRUE)
emm_sensanalysis_sample_type_poc <- data.frame(emmeans(lme_sensanalysis_group, ~time_hbin_factor, rg.limit = 60000, weights="prop"), 
                                         group = "Point-of-care test",
                                         npat = length(unique(df_sensanalysis_group$stay_id)),
                                         nmeas = nrow(df_sensanalysis_group))

df_sensanalysis_group <- df_gluc_data_analysis %>% filter(sample_type %in% c("Lab (serum)", "Lab (whole blood)"))  
lme_sensanalysis_group <- lmer(formula_model_final, data=df_sensanalysis_group, REML=TRUE)
emm_sensanalysis_sample_type_lab <- data.frame(emmeans(lme_sensanalysis_group, ~time_hbin_factor, rg.limit = 60000, weights="prop"), 
                                        group = "Lab test",
                                        npat = length(unique(df_sensanalysis_group$stay_id)),
                                        nmeas = nrow(df_sensanalysis_group))
emm_sensanalysis_sample_type <- bind_rows(emm_sensanalysis_sample_type_poc, emm_sensanalysis_sample_type_lab)

plt_text_sample_type <- emm_sensanalysis_sample_type %>% group_by(group, npat, nmeas) %>% 
  summarize(mean = mean(emmean), .groups="drop") %>% 
  mutate(y = case_when(mean == max(mean) ~ Inf, 
                       mean == min(mean) ~ -Inf,
                       TRUE ~ -Inf),
         vjust = case_when(mean == max(mean) ~ 1.2,
                           mean == min(mean) ~ -0.2,
                           TRUE ~ -1.4))



###############################################################################
##### Sensitivity analyses - Time in ICU days -----

df_sensanalysis_group <- df_gluc_data_analysis %>% filter(time_icudays <= 2)  
lme_sensanalysis_group <- lmer(formula_model_final, data=df_sensanalysis_group, REML=TRUE)
emm_sensanalysis_time_in_icu_0_to_3 <- data.frame(emmeans(lme_sensanalysis_group, ~time_hbin_factor, rg.limit = 60000, weights="prop"), 
                                               group = "Day 0-2",
                                               npat = length(unique(df_sensanalysis_group$stay_id)),
                                               nmeas = nrow(df_sensanalysis_group))

df_sensanalysis_group <- df_gluc_data_analysis %>% filter(time_icudays > 2 & time_icudays <= 7)  
lme_sensanalysis_group <- lmer(formula_model_final, data=df_sensanalysis_group, REML=TRUE)
emm_sensanalysis_time_in_icu_4_to_8 <- data.frame(emmeans(lme_sensanalysis_group, ~time_hbin_factor, rg.limit = 60000, weights="prop"), 
                                                  group = "Day 3-7",
                                                  npat = length(unique(df_sensanalysis_group$stay_id)),
                                                  nmeas = nrow(df_sensanalysis_group))

df_sensanalysis_group <- df_gluc_data_analysis %>% filter(time_icudays > 7 )  
lme_sensanalysis_group <- lmer(formula_model_final, data=df_sensanalysis_group, REML=TRUE)
emm_sensanalysis_time_in_icu_9_to_inf <- data.frame(emmeans(lme_sensanalysis_group, ~time_hbin_factor, rg.limit = 60000, weights="prop"), 
                                                  group = "Day ≥ 8",
                                                  npat = length(unique(df_sensanalysis_group$stay_id)),
                                                  nmeas = nrow(df_sensanalysis_group))

emm_sensanalysis_time_in_icu <- bind_rows(emm_sensanalysis_time_in_icu_0_to_3, emm_sensanalysis_time_in_icu_4_to_8, emm_sensanalysis_time_in_icu_9_to_inf)
emm_sensanalysis_time_in_icu$group <- factor(emm_sensanalysis_time_in_icu$group, levels = c(
  "Day 0-2",
  "Day 3-7",
  "Day ≥ 8"))

plt_text_time_in_icu <- emm_sensanalysis_time_in_icu %>% group_by(group, npat, nmeas) %>% 
  summarize(mean = mean(emmean), .groups="drop") %>% 
  mutate(y = case_when(mean == max(mean) ~ Inf, 
                       mean == min(mean) ~ -Inf,
                       TRUE ~ -Inf),
         vjust = case_when(mean == max(mean) ~ 1.2,
                           mean == min(mean) ~ -0.2,
                           TRUE ~ -1.4))



###############################################################################
##### Sensitivity analyses - Time to next measurement -----

df_sensanalysis_group <- df_gluc_data_analysis %>% filter(time_next_measurement < 4)  
lme_sensanalysis_group <- lmer(formula_model_final, data=df_sensanalysis_group, REML=TRUE)
emm_sensanalysis_time_to_next_0_to_4 <- data.frame(emmeans(lme_sensanalysis_group, ~time_hbin_factor, rg.limit = 60000, weights="prop"), 
                                                  group = "< 4 hours",
                                                  npat = length(unique(df_sensanalysis_group$stay_id)),
                                                  nmeas = nrow(df_sensanalysis_group))

df_sensanalysis_group <- df_gluc_data_analysis %>% filter(time_next_measurement >= 4 & time_next_measurement <= 8 )  
lme_sensanalysis_group <- lmer(formula_model_final, data=df_sensanalysis_group, REML=TRUE)
emm_sensanalysis_time_to_next_4_to_8 <- data.frame(emmeans(lme_sensanalysis_group, ~time_hbin_factor, rg.limit = 60000, weights="prop"), 
                                                   group = "4-8 hours",
                                                   npat = length(unique(df_sensanalysis_group$stay_id)),
                                                   nmeas = nrow(df_sensanalysis_group))

df_sensanalysis_group <- df_gluc_data_analysis %>% filter(time_next_measurement > 8)  
lme_sensanalysis_group <- lmer(formula_model_final, data=df_sensanalysis_group, REML=TRUE)
emm_sensanalysis_time_to_next_8_and_more <- data.frame(emmeans(lme_sensanalysis_group, ~time_hbin_factor, rg.limit = 60000, weights="prop"), 
                                                   group = "> 8 hours",
                                                   npat = length(unique(df_sensanalysis_group$stay_id)),
                                                   nmeas = nrow(df_sensanalysis_group))

emm_sensanalysis_time_to_next <- bind_rows(emm_sensanalysis_time_to_next_0_to_4, emm_sensanalysis_time_to_next_4_to_8, emm_sensanalysis_time_to_next_8_and_more)

emm_sensanalysis_time_to_next$group <- factor(emm_sensanalysis_time_to_next$group, levels = c(
  "< 4 hours",
  "4-8 hours",
  "> 8 hours"))


plt_text_time_to_next <- emm_sensanalysis_time_to_next %>% group_by(group, npat, nmeas) %>% 
  summarize(mean = mean(emmean), .groups="drop") %>% 
  mutate(y = case_when(mean == max(mean) ~ Inf, 
                       mean == min(mean) ~ -Inf,
                       TRUE ~ -Inf),
         vjust = case_when(mean == max(mean) ~ 1.2,
                           mean == min(mean) ~ -0.2,
                           TRUE ~ -1.4))


###############################################################################
##### Define plotting themes -----

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

plotcode_sensanalysis <- list(
  theme_plot,
  geom_hline(yintercept=c(7.7, 10), linetype="dashed", size=0.2, color="grey30"),
  geom_ribbon(alpha=0.45, color=NA),
  geom_line(size=0.2),
  geom_point(size=2),
  xlab("Time of day (h)"),
  ylab("Model-predicted glucose (mmol/L)"),
  coord_cartesian(ylim=c(7, 10.5), clip="off"),
  guides(fill="none"),
  scale_y_continuous(breaks=seq(7, 10.5, by = 0.5)),
  scale_x_continuous(expand=c(0,0.5), breaks=seq(0, 24, by = 4), limits=c(-0.2,24)))

color_3 <- "#117733" #"#185D8C"
color_1 <- "#882255"#228B22" #"#89B2D3"
color_2 <- "#88CCEE" #"#DC143C" #"#3489BE"
color_4 <- "#185D8C"

###############################################################################
##### Plots sensitivity analysis -----

plt_sensanalysis_ventilation <- ggplot(emm_sensanalysis_ventilation, aes(x=as.numeric(as.character(time_hbin_factor)), y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, 
                                                                         shape=group, color=group, fill=group, group=group))+
  theme_bw()+
  scale_fill_manual(values=c("Non-Invasive"=color_1,"Invasive"=color_2), name=NULL)+
  scale_shape_manual(values=c("Non-Invasive"=16,"Invasive"=15), name=NULL)+
  scale_color_manual(values=c("Non-Invasive"=color_1,"Invasive"=color_2), name=NULL)+
  ggtitle("Ventilation mode") +
  plotcode_sensanalysis
  #geom_text(data=plt_text_ventilation_mode, x=24,  aes(color=group, y=y, vjust=vjust, label=paste0(nmeas, " glucose measurements\n", npat, " patients")), hjust=1, size=3, inherit.aes = F, show.legend=F)
plt_sensanalysis_ventilation

plt_sensanalysis_survival <- ggplot(emm_sensanalysis_survival, aes(x=as.numeric(as.character(time_hbin_factor)), y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, 
                                                                      shape=group, color=group, fill=group, group=group))+theme_bw()+
  scale_fill_manual(values=c("Non-survivors"=color_1,"Survivors"=color_2), name=NULL)+
  scale_shape_manual(values=c("Non-survivors"=16,"Survivors"=15), name=NULL)+
  scale_color_manual(values=c("Non-survivors"=color_1,"Survivors"=color_2), name=NULL)+
  ggtitle("In-hospital mortality")+
  plotcode_sensanalysis
  #geom_text(data=plt_text_survival, x=24,  aes(color=group, y=y, vjust=vjust, label=paste0(nmeas, " glucose measurements\n", npat, " patients")), hjust=1, size=3, inherit.aes = F, show.legend=F)
plt_sensanalysis_survival

plt_sensanalysis_sedation <- ggplot(emm_sensanalysis_sedation, aes(x=as.numeric(as.character(time_hbin_factor)), y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, 
                                                                   shape=group, color=group, fill=group, group=group))+theme_bw()+
  scale_fill_manual(values=c("RASS ≥ -1"=color_1,"RASS ≤ -2"=color_2), name=NULL)+
  scale_shape_manual(values=c("RASS ≥ -1"=16,"RASS ≤ -2"=15), name=NULL)+
  scale_color_manual(values=c("RASS ≥ -1"=color_1,"RASS ≤ -2"=color_2), name=NULL)+
  ggtitle("Sedation depth")+
  plotcode_sensanalysis
  #geom_text(data=plt_text_sedation, x=24,  aes(color=group, y=y, vjust=vjust, label=paste0(nmeas, " glucose measurements\n", npat, " patients")), hjust=1, size=3, inherit.aes = F, show.legend=F)
plt_sensanalysis_sedation

plt_sensanalysis_sample_type <- ggplot(emm_sensanalysis_sample_type, aes(x=as.numeric(as.character(time_hbin_factor)), y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, 
                                                                         shape=group, color=group, fill=group, group=group))+theme_bw()+
  scale_fill_manual(values=c("Point-of-care test"=color_1,"Lab test"=color_2), name=NULL)+
  scale_shape_manual(values=c("Point-of-care test"=16,"Lab test"=15), name=NULL)+
  scale_color_manual(values=c("Point-of-care test"=color_1,"Lab test"=color_2), name=NULL)+
  ggtitle("Sample type")+
  plotcode_sensanalysis
  #geom_text(data=plt_text_sample_type, x=24,  aes(color=group, y=y, vjust=vjust, label=paste0(nmeas, " glucose measurements\n", npat, " patients")), hjust=1, size=3, inherit.aes = F, show.legend=F)
plt_sensanalysis_sample_type

plt_sensanalysis_time_in_icu <- ggplot(emm_sensanalysis_time_in_icu, aes(x=as.numeric(as.character(time_hbin_factor)), y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, 
                                                                         shape=group, color=group, fill=group, group=group))+theme_bw()+
  scale_fill_manual(values=c("Day 0-2"=color_1,"Day 3-7"=color_2,"Day ≥ 8"=color_3), name=NULL)+
  scale_shape_manual(values=c("Day 0-2"=16,"Day 3-7"=15,"Day ≥ 8"=18), name=NULL)+
  scale_color_manual(values=c("Day 0-2"=color_1,"Day 3-7"=color_2,"Day ≥ 8"=color_3), name=NULL)+
  ggtitle("Days in ICU")+
  plotcode_sensanalysis
  #geom_text(data=plt_text_time_in_icu, x=24,  aes(color=group, y=y, vjust=vjust, label=paste0(nmeas, " glucose measurements\n", npat, " patients")), hjust=1, size=3, inherit.aes = F, show.legend=F)
plt_sensanalysis_time_in_icu

plt_sensanalysis_time_to_next <- ggplot(emm_sensanalysis_time_to_next, aes(x=as.numeric(as.character(time_hbin_factor)), y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, 
                                                                           shape=group, color=group, fill=group, group=group))+theme_bw()+
  scale_fill_manual(values=c("< 4 hours"=color_1,"4-8 hours"=color_2,"> 8 hours"=color_3), name=NULL)+
  scale_shape_manual(values=c("< 4 hours"=16,"4-8 hours"=15,"> 8 hours"=18), name=NULL)+
  scale_color_manual(values=c("< 4 hours"=color_1,"4-8 hours"=color_2,"> 8 hours"=color_3), name=NULL)+
  ggtitle("Time to next glucose sample")+
  plotcode_sensanalysis
  #geom_text(data=plt_text_time_to_next, x=24,  aes(color=group, y=y, vjust=vjust, label=paste0(nmeas, " glucose measurements\n", npat, " patients")), hjust=1, size=3, inherit.aes = F, show.legend=F)
plt_sensanalysis_time_to_next


###############################################################################
##### Combine plots sensitivity analysis -----
plt_sensitivity_analysis <- cowplot::plot_grid(plt_sensanalysis_ventilation,
                                               plt_sensanalysis_survival,
                                               plt_sensanalysis_sedation,
                                               plt_sensanalysis_time_in_icu,
                                               plt_sensanalysis_sample_type,
                                               plt_sensanalysis_time_to_next,
                                               ncol=3, nrow=2, axis="tblr", labels=c("a","b","c","d", "e","f"), label_size = 10)
plt_sensitivity_analysis

ggsave(plt_sensitivity_analysis, filename=paste0("2-output/mimic-iv_24h-glucose-variation_04_sensitivity_analyses.png"), dpi=600, width=12, height=7, unit='in')
ggsave(plt_sensitivity_analysis, filename=paste0("2-output/mimic-iv_24h-glucose-variation_04_sensitivity_analyses.svg"), dpi=600, width=12, height=7, unit='in')
##ggsave(plt_sensitivity_analysis, filename=paste0("2-output/mimic-iv_24h-glucose-variation_04_sensitivity_analyses.pdf"), dpi=600, width=12, height=7, unit='in')
## NB: ≤ & ≥ signs do not work when saving to PDF --> First save to svg, than open in Illustrator, then save as PDF 




###############################################################################
##### Sensitivity analysis - Impaired glucose tolerance (Insulin administration) -----
df_gluc_data_analysis <- df_gluc_data_analysis %>% group_by(stay_id) %>% 
  mutate(insu_requirement = case_when(max(unitsperday_insu) == 0 ~ 'no', 
                                      mean(unitsperday_insu) > 0 & mean(unitsperday_insu) <= 35 ~ '0-35 units/day',
                                      mean(unitsperday_insu) > 35 & mean(unitsperday_insu) <= 70 ~ '35-70 units/day',
                                      mean(unitsperday_insu) > 70 ~ '>70 units/day')) %>% ungroup()


df_sensanalysis_group <- df_gluc_data_analysis %>% filter(insu_requirement == "no")  
lme_sensanalysis_group <- lmer(formula_model_final, data=df_sensanalysis_group, REML=TRUE)
emm_sensanalysis_insu_requirement_no <- data.frame(emmeans(lme_sensanalysis_group, ~time_hbin_factor, rg.limit = 60000, weights="prop"), 
                                                   group = "No insulin administration during ICU stay",
                                                   npat = length(unique(df_sensanalysis_group$stay_id)),
                                                   nmeas = nrow(df_sensanalysis_group))

df_sensanalysis_group <- df_gluc_data_analysis %>% filter(insu_requirement == "0-35 units/day")  
lme_sensanalysis_group <- lmer(formula_model_final, data=df_sensanalysis_group, REML=TRUE)
emm_sensanalysis_insu_requirement_0to35 <- data.frame(emmeans(lme_sensanalysis_group, ~time_hbin_factor, rg.limit = 60000, weights="prop"), 
                                                    group = "Mean insulin administration of 0<x≤35 units/day over ICU stay",
                                                    npat = length(unique(df_sensanalysis_group$stay_id)),
                                                    nmeas = nrow(df_sensanalysis_group))

df_sensanalysis_group <- df_gluc_data_analysis %>% filter(insu_requirement == "35-70 units/day")  
lme_sensanalysis_group <- lmer(formula_model_final, data=df_sensanalysis_group, REML=TRUE)
emm_sensanalysis_insu_requirement_35to70 <- data.frame(emmeans(lme_sensanalysis_group, ~time_hbin_factor, rg.limit = 60000, weights="prop"), 
                                                      group = "Mean insulin administration of 35<x≤70 units/day over ICU stay",
                                                      npat = length(unique(df_sensanalysis_group$stay_id)),
                                                      nmeas = nrow(df_sensanalysis_group))

df_sensanalysis_group <- df_gluc_data_analysis %>% filter(insu_requirement == ">70 units/day")  
lme_sensanalysis_group <- lmer(formula_model_final, data=df_sensanalysis_group, REML=TRUE)
emm_sensanalysis_insu_requirement_70plus <- data.frame(emmeans(lme_sensanalysis_group, ~time_hbin_factor, rg.limit = 60000, weights="prop"), 
                                                      group = "Mean insulin administration of >70 units/day over ICU stay",
                                                      npat = length(unique(df_sensanalysis_group$stay_id)),
                                                      nmeas = nrow(df_sensanalysis_group))
emm_sensanalysis_insu_requirement <- bind_rows(emm_sensanalysis_insu_requirement_no, emm_sensanalysis_insu_requirement_0to35, 
                                               emm_sensanalysis_insu_requirement_35to70, emm_sensanalysis_insu_requirement_70plus)

plt_text_insu_requirement <- emm_sensanalysis_insu_requirement %>% group_by(group, npat, nmeas) %>% 
  summarize(mean = mean(emmean), .groups="drop") %>% 
  mutate(y = case_when(mean == max(mean) ~ Inf, 
                       mean == min(mean) ~ -Inf,
                       TRUE ~ -Inf),
         vjust = case_when(mean == max(mean) ~ 1.2,
                           mean == min(mean) ~ -0.2,
                           TRUE ~ -1.4))


emm_sensanalysis_insu_requirement$group <- factor(emm_sensanalysis_insu_requirement$group, levels = c(
  "No insulin administration during ICU stay",
  "Mean insulin administration of 0<x≤35 units/day over ICU stay",
  "Mean insulin administration of 35<x≤70 units/day over ICU stay",
  "Mean insulin administration of >70 units/day over ICU stay"
))

plt_sensanalysis_insu_requirement <- ggplot(emm_sensanalysis_insu_requirement, aes(x=as.numeric(as.character(time_hbin_factor)), y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, shape=group, color=group, fill=group, group=group))+theme_bw()+
  scale_fill_manual(values=c("No insulin administration during ICU stay"=color_1,
                             "Mean insulin administration of 0<x≤35 units/day over ICU stay"=color_2,
                             "Mean insulin administration of 35<x≤70 units/day over ICU stay"=color_3,
                             "Mean insulin administration of >70 units/day over ICU stay"=color_4), name=NULL) +
  scale_shape_manual(values=c("No insulin administration during ICU stay"=17,
                              "Mean insulin administration of 0<x≤35 units/day over ICU stay"=18,
                              "Mean insulin administration of 35<x≤70 units/day over ICU stay"=16,
                              "Mean insulin administration of >70 units/day over ICU stay"=15), 
                              name=NULL)+
  scale_color_manual(values=c("No insulin administration during ICU stay"=color_1,
                              "Mean insulin administration of 0<x≤35 units/day over ICU stay"=color_2,
                              "Mean insulin administration of 35<x≤70 units/day over ICU stay"=color_3,
                              "Mean insulin administration of >70 units/day over ICU stay"=color_4), name=NULL) +
  #ggtitle("Insulin resistance")+
  plotcode_sensanalysis +
  coord_cartesian(ylim=c(6.5, 12), clip="off") + 
  scale_y_continuous(breaks=seq(6, 12, by = 1)) + 
  guides(
    color = guide_legend(nrow = 4, byrow = TRUE),
    shape = guide_legend(nrow = 4, byrow = TRUE),  # Add this line for shape legend
    fill = guide_legend(nrow = 4, byrow = TRUE),    # Add this line for fill legend
    labels = c(
      "No insulin administration during ICU stay",
      "Mean insulin administration of 0<x≤35 units/day over ICU stay",
      "Mean insulin administration of 35<x≤70 units/day over ICU stay",
      "Mean insulin administration of >70 units/day over ICU stay"
    ))
plt_sensanalysis_insu_requirement

ggsave(plt_sensanalysis_insu_requirement, filename=paste0("2-output/mimic-iv_24h-glucose-variation_04_sensitivity_analyses_insulin.png"), dpi=600, width=6, height=6, unit='in')
ggsave(plt_sensanalysis_insu_requirement, filename=paste0("2-output/mimic-iv_24h-glucose-variation_04_sensitivity_analyses_insulin.svg"), dpi=600, width=6, height=6, unit='in')
#ggsave(plt_sensanalysis_insu_requirement, filename=paste0("2-output/mimic-iv_24h-glucose-variation_04_sensitivity_analyses_insulin.pdf"), dpi=600, width=6, height=6, unit='in')
## NB: ≤ & ≥ signs do not work when saving to PDF --> First save to svg, than open in Illustrator, then save as PDF 



