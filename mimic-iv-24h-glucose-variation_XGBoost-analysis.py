"""
# Run-Time Environment: 	    Python version 3.10.2
# Authors:				        Floor Hiemstra
# Date created:                 August 2023
# Project:      			    24-hour glucose variation in ICU patients
# Filename:                     mimic-iv-24h-glucose-variation_XGBoost-analysis
# Purpose:                      24-hour analysis of glucose levels in ICU patients using XGBoost and SHAP
# Datafiles used:               mimic-iv_24h-glucose-variation_01_glucose-incl-covariates.csv (output of R script:
#                               mimic-iv-24h-glucose-variation_01_preprocessing.R
"""

#%% #################################################################
##### Import packages (see requirements.txt file for used package versions) -----
import pandas as pd
import numpy as np
import shap
import xgboost
from sklearn.model_selection import GridSearchCV
import os
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')

#%% #################################################################
##### Load data -----
file_path = r'C:\Users\fwhiemstra\Documents\mimic-glucose\2-output'
file_name = 'mimic-iv_24h-glucose-variation_01_glucose-incl-covariates.csv'
df_gluc = pd.read_csv(os.path.join(file_path, file_name))

#%% #################################################################
##### Select variables and prepare data for XGBoost model -----
## Select variables
gluc_level = df_gluc['glucose_value_mmolL']
time_variable = df_gluc['time_dec']
static_variables = df_gluc[['age_cat', 'gender', 'diabetes_status']]
dynamic_variables = df_gluc[["insu_yn", "dext_yn", "glucocorticoid_yn", "rate_nut_cat"]]

## Create dummy variables
static_variables_dummies = pd.get_dummies(static_variables, drop_first=True)
dynamic_variables_dummies = pd.get_dummies(dynamic_variables, drop_first=True)

## Combine
df_xgboost = pd.concat([gluc_level, time_variable, static_variables_dummies, dynamic_variables_dummies], axis=1).dropna()

## Prepare data
y, features = df_xgboost.iloc[:, 0], df_xgboost.iloc[:, 1:]

## Remove {, }, <, > from feature names (XGBoost cannot deal with this characters in feature names)
features_names = []
for i in range(0, len(features.columns)):
    new_feature_name = features.columns[i].replace("[", "(").replace("]", ")").replace(">", "→").replace("<", "←")
    features_names.append(new_feature_name)
column_mapping = dict(zip(features.columns, features_names))
features = features.rename(columns=column_mapping)

#%% #################################################################
##### Hyperparameter tuning -----
repeat_hyperparameter_tuning = 0

if repeat_hyperparameter_tuning == 1:
    import datetime
    starttime=datetime.datetime.now()
    print('Start time hyperparameter tuning: {}'.format(starttime))

    xgb_model = xgboost.XGBRegressor()

    param_grid = {'max_depth': [5, 10, 25, 50],
                  'eta': [0.05, 0.1, 0.3],
                  'n_estimators': [25, 50, 100, 150],
                  'min_child_weight': [1, 5, 10],
                  'colsample_bytree': [0.75, 1],
                  'subsample': [0.75, 1]}

    grid_search = GridSearchCV(estimator=xgb_model,
                               param_grid=param_grid,
                               scoring='neg_root_mean_squared_error',
                               cv=5,
                               verbose=3)
    grid_search.fit(features, y)

    params = grid_search.best_params_
    print('Best params: {}'.format(params))

    print('End time hyperparameter tuning: {}'.format(datetime.datetime.now()))
    print('Total time hyperparameter tuning: {}'.format(datetime.datetime.now()-starttime))

else:
    params = {'max_depth': 5,
              'eta': 0.05,
              'n_estimators': 150,
              'min_child_weight': 5,
              'colsample_bytree': 0.75,
              'subsample': 0.75}

#%% #################################################################
##### Model training -----

xgb_model = xgboost.XGBRegressor(n_estimators=params['n_estimators'], max_depth=params['max_depth'], eta=params['eta'],
                                 colsample_bytree=params['colsample_bytree'], min_child_weight=params['min_child_weight'],
                                 subsample=params['subsample'])
xgb_model.fit(features, y)


#%% #################################################################
##### Determine train score -----

y_pred = xgb_model.predict(features)
rmse_train = np.sqrt(np.mean((y - y_pred) ** 2))
print('RMSE = {}'.format(rmse_train))
print(params)

#%% #################################################################
##### SHAP analysis -----

## Set folder to save figures
folder_figures = r'C:\Users\fwhiemstra\Documents\mimic-glucose\1-scripts\mimic-glucose-python-xgboost\results'

## Some setting for plot lay-outs
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10

# Change feature names for plots:
features = features.rename(columns={'time_dec': 'Time (decimals)',
                                    'age_cat_(50,70)': 'Age 50-70 years',
                                    'age_cat_(70, Inf)': 'Age > 70 years',
                                    'gender_M': 'Male sex',
                                    'diabetes_status_y': 'Diabetes',
                                    'insu_yn_y': 'Insulin administration',
                                    'dext_yn_y': 'Dextrose administration',
                                    'glucocorticoid_yn_y': 'Glucocorticoids administration',
                                    'rate_nut_cat_→50 mL/h': 'Enteral nutrition rate > 50 mL/h'})

# Randomly select a sample of features for SHAP analysis
features_shapsample = features.sample(n=5000, random_state=42)
explainer = shap.TreeExplainer(xgb_model)
shap_values = explainer.shap_values(features_shapsample)
shap_outcome = explainer(features_shapsample)

plt.close('all')

##### SHAP over time -----
shap.dependence_plot(0, shap_values, features_shapsample, interaction_index=None,
                     color='black', dot_size=10)

fig2, ax = plt.gcf(), plt.gca()
fig2.set_size_inches(4.5, 4)  # Adjust the size as needed
ax.tick_params(labelsize=8)
ax.set_xlabel("Time of day (h)", fontsize=10)
ax.set_ylabel("SHAP value", fontsize=10)
plt.xticks([0,4,8,12,16,20,24])
plt.xlim([0, 24])

plt.tight_layout()
plt.show()

plt.savefig(os.path.join(folder_figures, 'XGboost-analysis_SHAP_time_of_day.png'), dpi=600)
plt.savefig(os.path.join(folder_figures, 'XGboost-analysis_SHAP_time_of_day.eps'), dpi=600)


##### SHAP summary plot -----
plt.figure()
shap.summary_plot(shap_values, features_shapsample, show=False)

fig1, ax = plt.gcf(), plt.gca()
fig1.set_size_inches(6, 4)  # Adjust the size as needed
ax.tick_params(labelsize=8)
ax.set_xlabel("SHAP value (impact on model output)", fontsize=10)

# Get colorbar
cb_ax = fig1.axes[1]
# Modifying color bar parameters
cb_ax.tick_params(labelsize=10)
cb_ax.set_ylabel("Feature value", fontsize=10)

plt.tight_layout()
plt.show()

plt.savefig(os.path.join(folder_figures, 'XGboost-analysis_SHAP_summary_plot.png'), dpi=600)
plt.savefig(os.path.join(folder_figures, 'XGboost-analysis_SHAP_summary_plot.eps'), dpi=600)


##### SHAP Global bar plot -----
plt.figure()
shap.plots.bar(shap_outcome)
fig3, ax = plt.gcf(), plt.gca()
fig3.set_size_inches(5, 4)  # Adjust the size as needed
ax.tick_params(labelsize=8)
ax.set_xlabel("mean(|SHAP value|)", fontsize=10)
ax.set_ylabel("", fontsize=10)

plt.tight_layout()
plt.show()

plt.savefig(os.path.join(folder_figures, 'XGboost-analysis_SHAP_global_bar_plot.png'), dpi=600)
plt.savefig(os.path.join(folder_figures, 'XGboost-analysis_SHAP_global_bar_plot.eps'), dpi=600)





