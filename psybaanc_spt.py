# -*- coding: utf-8 -*-
"""
Created on Thu Dec  4 13:10:16 2025
psybaanc sucrose preference test analysis

@author: olu

Create an excel file with your data, and the following columns:
    - mouse
    - sex
    - treatment
    - stress
    - sucrose_consumed
    - water_consumed
The path of this file goes in DATA_PATH.

"""
# %% Import packages
import os
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import functions.psybaanc_stats as psy_stats
import functions.psybaanc_plot as psy_plot

matplotlib.use("Agg")
# Dataframe print settings (do not change)
pd.set_option('display.max_rows', None)  # Display all rows
pd.set_option('display.max_columns', None)  # Display all columns
pd.set_option('display.width', 1000)  # Adjust display width to prevent line breaks
pd.set_option('display.max_colwidth', None)  # Display full content of each column

# %% Variables to change
INSTITUTION = "Berkeley 1"  # Stanford, Berkeley 1, Berkeley 2, UCSF 1, UCSF 2
DATA_PATH = r"Y:\PsyBAANC\paperExperiments\chronic CORT\SPT\SPT_edited_251204.xlsx"

# %% Get data
mouse_key = pd.read_excel(DATA_PATH)
mouse_id = mouse_key["mouse"].tolist()
sex_key = mouse_key["sex"].tolist()
treatment_key = mouse_key["treatment"].tolist()
stress_key = mouse_key["stress"].tolist()
sucrose_consumed = mouse_key["sucrose_consumed"]
water_consumed = mouse_key["water_consumed"]
total_consumed = sucrose_consumed + water_consumed
pref_score = sucrose_consumed/total_consumed

# %% Create data dictionary for summary data.
data_summary_dict = {
    'Animal_ID': mouse_id,
    'Sex': sex_key,
    'Treatment': treatment_key,
    'Stress': stress_key,
    'Institution': INSTITUTION,
    'Water_consumed': water_consumed,
    'Sucrose_consumed': sucrose_consumed,
    'Total_consumed': total_consumed,
    'Preference_score': pref_score
    }

data_summary_df = pd.DataFrame(data_summary_dict)
column_names = data_summary_df.columns[-4:].tolist()
column_labels = ["Water consumed (g)", "Sucrose consumed (g)", "Total consumed (g)",
                 "Preference Score"]
data_all = data_summary_df

# %% Do statistical tests.
stats_summary_results = {}
for col in column_names:
    stats_summary_results[col] = psy_stats.stats_treatment_sex_stress(data_summary_df, col)

# %% Plot all data
BASE_FOLDER = os.path.split(DATA_PATH)[0]
os.makedirs(os.path.join(BASE_FOLDER, "saved_data"), exist_ok=True)
# ymin = [0, 0, 0, 0]
# ymax = [4, 15, 20, 1]
fig, ax = plt.subplots(1, len(column_names), figsize=(1.2*len(column_names), 1.5))
for col_i, col in enumerate(column_names):
    psy_plot.plot_bars_thirdfactor(ax[col_i], data_summary_df, col, column_labels[col_i],
                                   INSTITUTION, stars=stats_summary_results[col]["significance"],
                                   # ymax=ymax[col_i], ymin=ymin[col_i]
                                   )
plt.savefig(os.path.join(BASE_FOLDER, "saved_data", "SPT_summary_data_plots.png"))
plt.close()

# %% Save raw data if wanted.
os.makedirs(os.path.join(BASE_FOLDER, "saved_data"), exist_ok=True)
data_all.to_csv(os.path.join(BASE_FOLDER, "saved_data", 'SPT_data.csv'), index=False)
