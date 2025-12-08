# -*- coding: utf-8 -*-
"""
Created on Tue Dec  2 12:38:47 2025

@author: olu
"""
# %% Import packages
import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import cv2
import functions.psybaanc_behavior as psy_beh
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
FOLDER_PATH = r"Y:/PsyBAANC/paperExperiments/chronic CORT/FST/cropped"  # path to folder with data
VIDEO_TYPE = "mp4"  # options: "mp4", "avi", others also likely OK.
COORDINATE_FILE_TYPE = "csv"  # options: "csv", "xlsx"

# ####### Create keys for animal conditions ##########
# EITHER write out condition in list format OR import data from excel file.

# Option to write out conditions in list format, or set variables as "None".
sex_key = None  # List of sex of the animals, i.e., ["M", "F", "M", etc.]
treatment_key = None  # List of treatment of animals.
stress_key = None  # List of animal's stress condition. Options: np.nan, "Ctrl", "Stress"

# Option to get list of mouse conditions from excel file. Or set mouse_key_path to "None".
mouse_key_path = r"Y:\PsyBAANC\paperExperiments\chronic CORT\Mouse_key.xlsx"  # edit as needed
if mouse_key_path is not None:
    mouse_key = pd.read_excel(mouse_key_path)
    sex_key = mouse_key["Sex"].tolist()
    treatment_key = mouse_key["Treatment"].tolist()
    stress_key = mouse_key["Stress"].tolist()

# %% Get paths for data
paths_csv = psy_beh.get_file_paths(FOLDER_PATH, "csv")
paths_vid = psy_beh.get_file_paths(FOLDER_PATH, VIDEO_TYPE)

# confirm accuracy of files and key information
csv_ids = [os.path.split(string)[-1] for string in paths_csv]
vid_ids = [os.path.split(string)[-1] for string in paths_vid]

files_analyzed = pd.DataFrame({"vid ID": vid_ids, "csv ID": csv_ids,
                               "sex": sex_key, "treatment": treatment_key})
print(f"""
      Please check that data is correct:  \n
      - vid ID and csv ID correspond to the same animal. \n
      - sex and treatment are assigned correctly.\n
      {files_analyzed}
      """)

while True:
    user_input = input("Type 'yes'/'no' to indicate whether data is accurate. \n")
    if user_input.lower() == 'yes':
        break  # Exit the loop
    elif user_input.lower() == 'no':
        raise ValueError("Please correct IDs or group assignments as needed and re-run script")

# %% Get framerate of all your video files
FRAMERATE = np.empty((len(paths_vid)))
for video_idx, video_path in enumerate(paths_vid):
    cap = cv2.VideoCapture(video_path)
    FRAMERATE[video_idx] = cap.get(cv2.CAP_PROP_FPS)

print(f"Note: the framerate of your files is {FRAMERATE}.")

# %% Main logic
immobility = np.empty((len(paths_csv)))
struggling_bouts_n = np.empty((len(paths_csv)))
struggling_bouts_time_avg = np.empty((len(paths_csv)))
for path_i, path in enumerate(paths_csv):
    data = pd.read_csv(path)
    data = data.iloc[:, 1]
    total_time_scored = np.sum(~np.isnan(data))/FRAMERATE[path_i]
    struggling_time = np.nansum(data)/FRAMERATE[path_i]
    immobility[path_i] = total_time_scored - struggling_time

    struggling_frames = np.where(data == 1)[0]
    struggling_bouts = psy_beh.split_into_visits(struggling_frames)
    struggling_bouts_n[path_i] = len(struggling_bouts)
    struggling_bouts_time = [len(bout)/FRAMERATE[path_i] for bout in struggling_bouts]
    struggling_bouts_time_avg[path_i] = np.mean(struggling_bouts_time)

    print(f"Video no. {path_i+1}: Scoring time: {total_time_scored}."
          f"Immobility: {immobility[path_i]}")

# %% Create data dictionary for summary data.
data_summary_dict = {
    'Animal_ID': paths_vid,
    'Sex': sex_key,
    'Treatment': treatment_key,
    'Stress': stress_key,
    'Institution': INSTITUTION,
    'Immobility': immobility,
    'Struggling_bouts': struggling_bouts_n,
    'Struggling_avg_time': struggling_bouts_time_avg,
    }

data_summary_df = pd.DataFrame(data_summary_dict)
column_names = data_summary_df.columns[-3:].tolist()
column_labels = ["Immobility (s)", "Number of struggling bouts", "Average time struggling (s)"]
data_all = data_summary_df

# %% Do statistical tests.
stats_summary_results = {}
for col in column_names:
    stats_summary_results[col] = psy_stats.stats_treatment_sex_stress(data_summary_df, col)

# %% Plot all data
os.makedirs(os.path.join(FOLDER_PATH, "saved_data"), exist_ok=True)
# ymin = [100, 0, 0]
# ymax = [300, 40, 15]
fig, ax = plt.subplots(1, len(column_names), figsize=(1.2*len(column_names), 1.5))
for col_i, col in enumerate(column_names):
    psy_plot.plot_bars_thirdfactor(ax[col_i], data_summary_df, col, column_labels[col_i],
                                   INSTITUTION, stats_summary_results[col]["significance"],
                                   # ymin=ymin[col_i], ymax=ymax[col_i]
                                   )
plt.savefig(os.path.join(FOLDER_PATH, "saved_data", "FST_summary_data_plots.png"))
plt.close()

# %% Save raw data if wanted.
os.makedirs(os.path.join(FOLDER_PATH, "saved_data"), exist_ok=True)
data_all.to_csv(os.path.join(FOLDER_PATH, "saved_data", 'FST_data.csv'), index=False)
