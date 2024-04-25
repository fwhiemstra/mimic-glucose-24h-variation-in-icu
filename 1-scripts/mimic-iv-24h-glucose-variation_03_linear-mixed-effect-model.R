# Run-Time Environment: 	R version 4.3.1
# Authors:				        Floor Hiemstra & Laura Kervezee
# Date created:           January 2024
# Project:      			    24-hour glucose variation in ICU patients
# Filename:               mimic-iv-24h-glucose-variation_03_linear-mixed-effect-model
# Purpose:  			        Analysis of 24-h variation in glucose levels: 
#                             - Visualization of unadjusted 24-hour pattern of glucose levels
#                             - Linear mixed effect models to determine effect of time of day on glucose levels during enteral tube feeding
#                             - Estimated marginal means are used for visualization
#                             - Model diagnostics of linear mixed effect models 
# Dependencies:          mimic-iv-24h-glucose-variation_01_preprocessing			    

########################################################################################################################
##### Packages -----
library(tidyverse)
library(stringr)
library(lubridate)
library(lmerTest)
library(emmeans)
library(sjPlot)
library(car)
library(svglite)
library(cowplot)

########################################################################################################################
##### Load data from csv -----
df_gluc_nut_variables <- read_csv("2-output/mimic-iv_24h-glucose-variation_01_glucose-incl-covariates.csv", lazy=F, show_col_types=FALSE)

## OR ##

########################################################################################################################
##### Load data from script -----
#source("C:\Users\fwhiemstra\Documents\mimic-glucose\1-scripts\R script versie 4\mimic-iv-24h-glucose-variation_01_preprocessing.R")

########################################################################################################################
##### STEP 1: Define plotting theme  -----

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


########################################################################################################################
##### STEP 2: Prepare data  -----

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
         ras_score) 

# Set reference categories
df_gluc_data_analysis$diabetes_status = relevel(df_gluc_data_analysis$diabetes_status, ref = "n")
df_gluc_data_analysis$sex = relevel(df_gluc_data_analysis$sex, ref = "F")
df_gluc_data_analysis$glucocorticoid_yn = relevel(df_gluc_data_analysis$glucocorticoid_yn, ref = "n")
df_gluc_data_analysis$age_cat = relevel(df_gluc_data_analysis$age_cat, ref = "[-Inf,55]")
df_gluc_data_analysis$rate_insu_cat = relevel(df_gluc_data_analysis$rate_insu_cat, ref = "[-Inf,0]")
df_gluc_data_analysis$rate_dext_cat = relevel(df_gluc_data_analysis$rate_dext_cat, ref = "[-Inf,0]")
df_gluc_data_analysis$rate_nut_cho_cat = relevel(df_gluc_data_analysis$rate_nut_cho_cat, ref = "[-Inf,4.5]")


df_gluc_data_model <- df_gluc_data_analysis %>% 
  select(glucose_value_mmolL, 
         time_hbin, time_hbin_factor, 
         stay_id, 
         age_cat,  
         sex, 
         diabetes_status, 
         rate_insu_cat, 
         rate_dext_cat, 
         glucocorticoid_yn, 
         rate_nut_cho_cat)

########################################################################################################################
##### STEP 3: Plot raw 24-h pattern (per-patient normalized) -----

## Raw plot normalized per patient
plt_raw_normalized <-  df_gluc_data_analysis %>% group_by(stay_id) %>% 
  # Group by stay_id and calculate change from baseline per patientgroup_by(stay_id) %>% 
  mutate(value_bslcorr = glucose_value_mmolL - mean(glucose_value_mmolL)) %>% 
  ungroup() %>% 
  # Plot by time of day
  ggplot(aes(x=time_hbin, y=value_bslcorr))+
  theme_bw()+
  theme_plot+
  stat_summary(fun.data="mean_cl_normal", fun.args=list(conf.int=.95), geom="ribbon", alpha=0.3, color=NA)+
  stat_summary(fun.data="mean_cl_normal", fun.args=list(conf.int=.95), geom="line")+
  stat_summary(fun.data="mean_cl_normal", fun.args=list(conf.int=.95), geom="point")+
  xlab("Time of day (h)")+
  ylab("Glucose (mmol/L) (change from patient average)") + 
  scale_x_continuous(expand=c(0,0.5), breaks=seq(0, 24, by = 4), limits=c(-0.2,24))

plt_raw_normalized

ggsave(plt_raw_normalized, filename=paste0("2-output/mimic-iv_24h-glucose-variation_03_normalized-24h-pattern.png"), dpi=600, width=4.5, height=4)
ggsave(plt_raw_normalized, filename=paste0("2-output/mimic-iv_24h-glucose-variation_03_normalized-24h-pattern.svg"), dpi=600, width=4.5, height=4)

###############################################################################
##### STEP 4: Linear mixed effect models -----

#######################################
## Linear mixed effect models ---

## Define variables 
dependent_variable <- "glucose_value_mmolL"
time_variable <- "time_hbin_factor"
patient_level_variables <- c("diabetes_status", "age_cat", "sex")
sample_level_variables <- c("rate_insu_cat", "rate_dext_cat", "glucocorticoid_yn", "rate_nut_cho_cat")
random_variable <- "(1 | stay_id)"

## Model 1
formula_model_1 <- as.formula(paste(dependent_variable, "~", "1 +", random_variable))
lmemodel_gluc_time_model_1 <- lmer(formula_model_1, data=df_gluc_data_model, REML=TRUE) # NB: REML=TRUE is the default setting

## Model 2 (Patient-level variables)
formula_model_2 <- as.formula(paste(dependent_variable, "~", paste(patient_level_variables, collapse=" + "), " + ", random_variable))
lmemodel_gluc_time_model_2 <- lmer(formula_model_2, data=df_gluc_data_model, REML=TRUE)

## Model 3 (Patient-level + sample-level variables)
formula_model_3 <-  as.formula(paste(dependent_variable, "~", paste(patient_level_variables, collapse=" + "), " + ", paste(sample_level_variables, collapse=" + "), " + ", random_variable))
lmemodel_gluc_time_model_3 <- lmer(formula_model_3, data=df_gluc_data_model, REML=TRUE)

## Model 4 (Patient-level + sample-level + time variables)
formula_model_4 <-  as.formula(paste(dependent_variable, "~", time_variable, " + ", paste(patient_level_variables, collapse=" + "), " + ", paste(sample_level_variables, collapse=" + "), " + ", random_variable))
lmemodel_gluc_time_model_4 <- lmer(formula_model_4, data=df_gluc_data_model, REML=TRUE)

#######################################
## Goodness of fit ---

# Calculate goodness of fit (log-likelihood, AIC, BIC) and compare model fits with log-likelihood ratio test 
aov_model_1_model_2 <- anova(lmemodel_gluc_time_model_1, lmemodel_gluc_time_model_2, refit = TRUE) # Refitting with ML instead of REML
aov_model_2_model_3 <- anova(lmemodel_gluc_time_model_2, lmemodel_gluc_time_model_3, refit = TRUE)
aov_model_3_model_4 <- anova(lmemodel_gluc_time_model_3, lmemodel_gluc_time_model_4, refit = TRUE)

## Root mean squared error 
rmse <- data.frame(model_1 = sqrt(mean(residuals(lmemodel_gluc_time_model_1)^2)),
                   model_2 = sqrt(mean(residuals(lmemodel_gluc_time_model_2)^2)),
                   model_3 = sqrt(mean(residuals(lmemodel_gluc_time_model_3)^2)),
                   model_4 = sqrt(mean(residuals(lmemodel_gluc_time_model_4)^2)))


###############################################################################
##### STEP 5: Final linear mixed effect models -----

lmemodel_gluc_time_model_final <- lmemodel_gluc_time_model_4

#######################################
## Goodness of fit ---
summary(lmemodel_gluc_time_model_final)
AIC_final_model <- AIC(lmemodel_gluc_time_model_final)
BIC_final_model <- BIC(lmemodel_gluc_time_model_final)
logLik_final_model <- logLik(lmemodel_gluc_time_model_final)
rmse_final_model <- sqrt(mean(residuals(lmemodel_gluc_time_model_final)^2))

#######################################
## Weight coefficients ---
coefs_model_final <- data.frame(summary(lmemodel_gluc_time_model_final)$coefficients)
coefs_model_final$variables <- rownames(coefs_model_final)


plt_coefs_final <- ggplot(data=coefs_model_final, aes(x=variables, y=Estimate)) + 
  theme_bw() + 
  geom_bar(stat='identity', position='dodge') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
plt_coefs_final

#######################################
## Estimated marginal means ---

emm_gluc_time_final_model <- data.frame(emmeans(lmemodel_gluc_time_model_final, ~time_hbin_factor, rg.limit = 90000, weights="prop"), 
                                 model = "Linear-mixed effect model",
                                 group="all data", 
                                 npat = length(unique(df_gluc_data_model$stay_id)),
                                 nmeas = nrow(df_gluc_data_model))

plt_gluc_time_final_emm <- ggplot(emm_gluc_time_final_model, aes(x=as.numeric(as.character(time_hbin_factor)), 
                                                                 y=emmean, ymin=asymp.LCL, ymax=asymp.UCL)) +
  theme_bw()+
  theme_plot + 
  geom_ribbon(alpha=0.3, color=NA) + 
  geom_line() + 
  geom_point(size=1)+ 
  xlab("Time of day (h)")+ 
  ylab("Model-predicted glucose (mmol/L)")+ 
  coord_cartesian(ylim=c(8.3, 9.6), clip="off") + 
  guides(fill="none") + 
  scale_y_continuous(breaks=seq(8.3, 9.6, by = 0.25))+ 
  scale_x_continuous(expand=c(0,0.5), breaks=seq(0, 24, by = 4), limits=c(-0.2,24))
  #ggtitle("Estimated marginal means \n Final linear-mixed effects model")
plt_gluc_time_final_emm

ggsave(plt_gluc_time_final_emm, filename=paste0("2-output/mimic-iv_24h-glucose-variation_03_lmemodel_gluc_time_final_emm.png"),
       dpi=600, width=4.5, height=4, units="in")
ggsave(plt_gluc_time_final_emm, filename=paste0("2-output/mimic-iv_24h-glucose-variation_03_lmemodel_gluc_time_final_emm.svg"),
       dpi=600, width=4.5, height=4, units="in")

#######################################
## Model diagnostics ---
## Plot distribution of residuals 
plt_hist_residuals <- ggplot(data.frame(residuals=residuals(lmemodel_gluc_time_model_final)), aes(x=residuals)) + 
  theme_bw()+
  theme_plot+
  #ggtitle("Distribution of residuals (LME model)") +
  geom_histogram(color='black', fill='darkgrey', binwidth=0.5) + 
  #scale_x_continuous(breaks=seq(-20, 45, by = 10))+ 
  #scale_y_continuous(breaks=seq(0, 30000, by = 5000)) + 
  xlab('Residuals') + ylab('Count')
plt_hist_residuals

## Residuals versus fitted values
plt_residuals_versus_fitted <- ggplot(data.frame(residuals=residuals(lmemodel_gluc_time_model_final), fitted=predict(lmemodel_gluc_time_model_final)), 
                                      aes(x=fitted, y=residuals)) + 
  theme_bw()+
  theme_plot+
  #ggtitle("Residuals versus fitted values (LME model)") +
  geom_point(size=2, shape=1, color="#104060") + 
  geom_hline(yintercept=0, color='black', linetype="dashed") + 
  #scale_y_continuous(breaks=seq(-20, 45, by = 10))+ 
  #scale_x_continuous(breaks=seq(5, 25, by = 5)) + 
  ylab('Residuals') + xlab('Fitted values')
plt_residuals_versus_fitted

plt_hist_and_residuals_versus_fitted <- cowplot::plot_grid(plt_residuals_versus_fitted, plt_hist_residuals,
                                                           ncol=2, nrow=1, axis="tblr", labels=c("a","b"), label_size = 10)

ggsave(plt_hist_and_residuals_versus_fitted, filename=paste0("2-output/mimic-iv_24h-glucose-variation_03_lmemodel_gluc_time_final_residuals_versus_fitted_and_hist_residuals.png"), 
       dpi=600, width=9, height=4)
ggsave(plt_hist_and_residuals_versus_fitted, filename=paste0("2-output/mimic-iv_24h-glucose-variation_03_lmemodel_gluc_time_final_residuals_versus_fitted_and_hist_residuals.svg"), 
       dpi=600, width=9, height=4)


#######################################
## Save LME model ---
tab_model(lmemodel_gluc_time_model_final, file=paste0("2-output/mimic-iv_24h-glucose-variation_03_lmemodel_gluc_time_final.html"))



