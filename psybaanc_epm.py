# -*- coding: utf-8 -*-
"""
Created on Wed Jul  2 16:01:40 2025

@author: olu

EPM Analysis
Instructions:
    1. Put all video and coordinate files from a single experiment in a folder.
    2. Adjust the variables in the "Variables to change" section to match your settings.
    3. Run the script. Follow command-line instructions as prompted.
    4. Update Supplementary table 1 (meta data sheet) with data in the "data_all" variable.

The variables in data_all include summary data for the following measures:
    'Time_Open': (%) time in open arms
    'Time_Center': (%) time in center of epm
    'Latency_Open': (%) time it took for mouse to first enter open arm
    'Distance': (cm) distance travelled in epm

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
FOLDER_PATH = r"Y:/PsyBAANC/paperExperiments/EPMAcute/Videos/all"  # path to folder with data
VIDEO_TYPE = "mp4"  # options: "mp4", "avi", others also likely OK.
COORDINATE_FILE_TYPE = "csv"  # options: "csv", "xlsx"
START_SEC = 0  # the time in seconds that you wish to begin analysis.
END_SEC = 10*60  # the time in seconds that you wish to end analysis.

LENGTH_CM = 75.5  # true size (cm) of one dimension of your epm (used for pixel to cm calibration)

X_COORDINATE_INDEX = 13  # Index of x-coord column in coordinates file (index starts at 0).
Y_COORDINATE_INDEX = 14  # Index of y-coord column in coordinates file (index starts at 0).
ROW_INDEX = 4  # Row number that position data starts in coordinate files
DATA_DLC = True  # Is your data from deeplabcut (DLC)? true or false.
COORDINATES_CM = False  # Are your coordinates in centimeters? (And not pixels)

sex_key = ["M"]*20 + ["F"]*20  # List of sex of the animals, i.e., ["M", "F", "M", etc.]
treatment_key = (["P"]*5 + ["S"]*5 + ["S"]*5 + ["P"]*5 +
                 ["S"]*5 + ["P"]*5 + ["P"]*5 + ["S"]*5)  # List of treatment of animals.
stress_key = [np.nan]*40

# %% Get paths of videos and csv files
paths_coordinates = psy_beh.get_file_paths(FOLDER_PATH, COORDINATE_FILE_TYPE)
paths_vid = psy_beh.get_file_paths(FOLDER_PATH, VIDEO_TYPE)

# confirm accuracy of files and key information
coordinate_ids = [os.path.split(string)[-1] for string in paths_coordinates]
vid_ids = [os.path.split(string)[-1] for string in paths_vid]

files_analyzed = pd.DataFrame({"vid ID": vid_ids, "coordinate ID": coordinate_ids,
                               "sex": sex_key, "treatment": treatment_key})
print(f"""
      Please check that data is correct:  \n
      - vid ID and coordinate ID correspond to the same animal. \n
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

# %% Read in all of the coordinate data.
print("""
      READING COORDINATE DATA
      """)
body_coords = psy_beh.read_coordinates(paths_coordinates,
                                       X_COORDINATE_INDEX, Y_COORDINATE_INDEX, ROW_INDEX,
                                       DATA_DLC=DATA_DLC)

if COORDINATES_CM:
    cm_to_pixels = np.empty(len(paths_coordinates))
    for file_idx, coordinates in enumerate(body_coords):
        cm_to_pixels[file_idx] = psy_beh.calibrate_pixels_to_cm(paths_vid[file_idx],
                                                                real_world_cm=LENGTH_CM,
                                                                frame_number=100)
        body_coords[file_idx] = body_coords[file_idx]/cm_to_pixels[file_idx]

# %% Define all relevant ROIs if you have not already.
# Note: ROIs are saved for each video after defining them.
# If you want to re-define an ROI, you must delete the existing one.

# Preallocate lists.
epm_roi_open_arm = [None]*(len(paths_vid))
epm_roi_closed_arm = [None]*(len(paths_vid))
epm_roi_center = [None]*(len(paths_vid))

# Define open and closed arms.
print("""
      Draw a polygon ROI to define the open/closed arms or read in existing ROIs.
      """)
for video_idx, video_file in enumerate(paths_vid):
    COORDINATES_PLOTTED = body_coords[video_idx]
    epm_roi_open_arm[video_idx] = psy_beh.get_roi_mask(video_file, "open_arm", "polygon",
                                                       COORDINATES_PLOTTED)
    epm_roi_closed_arm[video_idx] = psy_beh.get_roi_mask(video_file, "closed_arm", "polygon",
                                                         COORDINATES_PLOTTED)

    # define overlap between closed and open arms as center.
    epm_roi_center[video_idx] = cv2.bitwise_and(epm_roi_open_arm[video_idx],
                                                epm_roi_closed_arm[video_idx])
    epm_roi_open_arm[video_idx] = cv2.subtract(epm_roi_open_arm[video_idx],
                                               epm_roi_center[video_idx])
    epm_roi_closed_arm[video_idx] = cv2.subtract(epm_roi_closed_arm[video_idx],
                                                 epm_roi_center[video_idx])

# %% Main analysis
LENGTH_EXPERIMENT = END_SEC-START_SEC

frames_open = [None]*len(paths_vid)
frames_open_visits = [None]*len(paths_vid)

time_open = np.empty(len(paths_vid))
visits_open = np.empty(len(paths_vid))
latency_open = np.empty(len(paths_vid))
time_center = np.empty(len(paths_vid))

distance_travelled = np.empty(len(paths_vid))
velocity = np.empty(len(paths_vid))
cm_to_pixels = np.empty(len(paths_vid))

for video_idx, video_path in enumerate(paths_vid):
    START_FRAME = int(round(START_SEC*FRAMERATE[video_idx]))
    END_FRAME = int(round(END_SEC*FRAMERATE[video_idx])) - 1
    frames_open[video_idx] = psy_beh.get_timepoints_in_mask(body_coords[video_idx],
                                                            epm_roi_open_arm[video_idx])
    frames_open[video_idx] = frames_open[video_idx][(frames_open[video_idx] >= START_FRAME) &
                                                    (frames_open[video_idx] <= END_FRAME)]
    frames_open_visits[video_idx] = psy_beh.split_into_visits(frames_open[video_idx],
                                                              min_length=FRAMERATE[video_idx]/6)

    time_open[video_idx] = len(frames_open[video_idx])/FRAMERATE[video_idx]/LENGTH_EXPERIMENT*100
    visits_open[video_idx] = len(frames_open_visits[video_idx])
    if len(frames_open_visits[video_idx]) > 0:
        latency_open[video_idx] = frames_open_visits[video_idx][0][0]/FRAMERATE[video_idx]
    else:
        latency_open[video_idx] = END_FRAME/FRAMERATE[video_idx]
    frames_center = psy_beh.get_timepoints_in_mask(body_coords[video_idx],
                                                   epm_roi_center[video_idx])
    frames_center = frames_center[(frames_center >= START_FRAME) & (frames_center <= END_FRAME)]
    time_center[video_idx] = len(frames_center)/FRAMERATE[video_idx]/LENGTH_EXPERIMENT*100

    # Get distance travelled
    cm_to_pixels[video_idx] = psy_beh.calibrate_pixels_to_cm(video_path,
                                                             real_world_cm=LENGTH_CM,
                                                             frame_number=100)

    # get distance array: distance travelled per frame.
    diff_array = np.diff(body_coords[video_idx], axis=0)
    diff_array = diff_array * cm_to_pixels[video_idx]
    distance_array = np.sqrt((diff_array[:, 0]**2) + (diff_array[:, 1]**2))
    distance_array = distance_array[START_FRAME:END_FRAME+1]

    # Summary locomotion metrics
    distance_travelled[video_idx] = sum(distance_array)
    velocity[video_idx] = distance_travelled[video_idx]/LENGTH_EXPERIMENT

# %% Create data dictionary for summary data.
data_summary_dict = {
    'Animal_ID': paths_vid,
    'Sex': sex_key,
    'Treatment': treatment_key,
    'Stress': stress_key,
    'Institution': INSTITUTION,
    'Time_Open': time_open,
    'Time_Center': time_center,
    'Latency_Open': latency_open,
    'Distance': distance_travelled,
    }

data_summary_df = pd.DataFrame(data_summary_dict)
column_names = data_summary_df.columns[-4:].tolist()
column_labels = ["Time in open arms (%)", "Time in center (%)", "Latency open arm entry (s)",
                 "Distance travelled (cm)"]
data_all = data_summary_df

# %% Do statistical tests.
stats_summary_results = {}
for col in column_names:
    if "Stress" not in stress_key:
        stats_summary_results[col] = psy_stats.stats_treatment_sex(data_summary_df, col)
    else:
        stats_summary_results[col] = psy_stats.stats_treatment_sex_stress(data_summary_df, col)
# %% Plot all data
os.makedirs(os.path.join(FOLDER_PATH, "saved_data"), exist_ok=True)

fig, ax = plt.subplots(1, len(column_names), figsize=(1.2*len(column_names), 1.5))
for col_i, col in enumerate(column_names):
    if "Stress" not in stress_key:
        psy_plot.plot_individual_lab(ax[col_i], data_summary_df, col, column_labels[col_i],
                                     INSTITUTION, stats_summary_results[col]["significance"])
    else:
        psy_plot.plot_bars_thirdfactor(ax[col_i], data_summary_df, col, column_labels[col_i],
                                       INSTITUTION, stats_summary_results[col]["significance"])
plt.savefig(os.path.join(FOLDER_PATH, "saved_data", "EPM_summary_data_plots.png"))
plt.close()

# %% Plot animal traces for visualization.
for video_idx, video_path in enumerate(paths_vid):
    mouse_name = os.path.splitext(vid_ids[video_idx])[0]
    path_trace = os.path.join(FOLDER_PATH, "saved_data", f"{mouse_name}_trace.png")
    if os.path.exists(path_trace) is False:
        fig, ax = plt.subplots()
        psy_plot.plot_traces(ax, video_path, body_coords[video_idx])
        ax.imshow(epm_roi_center[video_idx], cmap='Greens', alpha=0.25)
        ax.imshow(epm_roi_open_arm[video_idx], cmap='Reds', alpha=0.25)
        ax.imshow(epm_roi_closed_arm[video_idx], cmap='Purples', alpha=0.25)
        plt.savefig(path_trace)
        plt.close()

# %% Save raw data if wanted.
os.makedirs(os.path.join(FOLDER_PATH, "saved_data"), exist_ok=True)
data_all.to_csv(os.path.join(FOLDER_PATH, "saved_data", 'EPM_data.csv'), index=False)
