# -*- coding: utf-8 -*-
"""
Created on Mon Jun 30 09:41:34 2025 / Edited by cliu to autosave csv files
@author: olu

OFT Analysis
Instructions:
    1. Put all video and coordinate files from a single experiment in a folder.
    2. Adjust the variables in the "Variables to change" section to match your settings.
    3. Run the script. Follow command-line instructions as prompted.
    4. Update Supplementary table 1 (meta data sheet) with data in the OFT_data saved csv file.

Results are stored in the experiment's folder path, within the saved_data directory and include:
    - plots of results
    - plots of behavior traces
    - OFT_data.csv file with animal data

The variables saved include summary and time-binned data for the following measures:
    'Time_Center': (%) time an animal spends in the center
    'Time_Corners': (%) time an animal spends in the corners
    'Distance_in_center': (cm) distance travelled in the center
    'Velocity': (cm/s) average velocity in open field
    'Time_moving': (%) time an animal spends moving (defined as movement >= 5 cm/s)
    'Velocity_while_moving': (cm/s) average velocity during time points considered as moving

    Additional (optional) variables outputted from script:
    'Time_Edges': (%) time an animal spends in the edge
    'Distance': (cm) distance travelled in open field
    'Time_running': (%) time an animal spends running (defined as movement > 20 cm/s)
    'Time_walking': (%) time an animal spends walking (defined as movement >= 5 or < 20 cm/s)
    'Time_freezing': (%) time an animal spends freezing (defined as movement < 5 cm/s)

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

matplotlib.use('Agg')

# %% Variables to change
INSTITUTION = "Berkeley 1"  # Stanford, Berkeley 1, Berkeley 2, UCSF 1, UCSF 2
FOLDER_PATH = r"Y:/PsyBAANC/paperExperiments/OFT_EPM/OFT_all"  # path to folder with data
VIDEO_TYPE = "avi"  # options: "mp4", "avi", others also likely OK.
COORDINATE_FILE_TYPE = "csv"  # options: "csv", "xlsx"
START_SEC = 0  # the time in seconds that you wish to begin analysis.
END_SEC = 30*60  # the time in seconds that you wish to end analysis.
INTERVAL_SEC = 5*60  # for time-binned measures, what is the bin interval, in seconds?

LENGTH_CM = [50]*40  # list of true size of your open field box in cm, for each video
N_BOXES = 25  # How many squares do you want to divide OF into? Keep at 25.

X_COORDINATE_INDEX = 25  # Index of x-coord column in coordinates file (index starts at 0).
Y_COORDINATE_INDEX = 26  # Index of y-coord column in coordinates file (index starts at 0).
ROW_INDEX = 4  # Row number that position data starts in coordinate files
DATA_DLC = True  # Is your data from deeplabcut (DLC)? true or false.
COORDINATES_CM = False  # Are your coordinates in centimeters? (And not pixels)

sex_key = ["M"]*20 + ["F"]*20  # List of sex of the animals, i.e., ["M", "F", "M", etc.]
treatment_key = (["P"]*5 + ["S"]*5 + ["P"]*5 + ["S"]*5 +
                 ["S"]*5 + ["P"]*5 + ["S"]*5 + ["P"]*5)  # List of treatment of animals.

# Dataframe print settings (do not change)
pd.set_option('display.max_rows', None)  # Display all rows
pd.set_option('display.max_columns', None)  # Display all columns
pd.set_option('display.width', 1000)  # Adjust display width to prevent line breaks
pd.set_option('display.max_colwidth', None)  # Display full content of each column

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
                                                                real_world_cm=LENGTH_CM[file_idx],
                                                                frame_number=0)
        body_coords[file_idx] = body_coords[file_idx]/cm_to_pixels[file_idx]

# %% Define all relevant ROIs if you have not already.
# Note: ROIs are saved for each video after defining them.
# If you want to re-define an ROI, you must delete the existing one.

# Preallocate lists.
open_field_base = [None]*(len(paths_vid))
edge_rois = [None]*(len(paths_vid))
center_rois = [None]*(len(paths_vid))
corner_rois = [None]*(len(paths_vid))

# First, draw a rectangle that defines the base of the open field.
print("""
      FIRST: Draw a rectangle ROI to define the base of the open field or read in existing ROI.
      """)
for video_idx, video_file in enumerate(paths_vid):
    open_field_base[video_idx] = psy_beh.get_roi_coords(video_file, "roi", body_coords[video_idx])

# Then, divide the open field base into N_BOXES number of equal portions.
# From the boxes, define edge, center, and corner ROIs.
for video_idx, video_file in enumerate(paths_vid):
    open_field_width = open_field_base[video_idx][1] - open_field_base[video_idx][0]
    open_field_length = open_field_base[video_idx][3] - open_field_base[video_idx][2]
    cm_to_pixels_x = LENGTH_CM[video_idx]/open_field_width
    cm_to_pixels_y = LENGTH_CM[video_idx]/open_field_length

    n_rows = int(np.sqrt(N_BOXES))
    sub_width = open_field_width / n_rows  # width of each box.
    sub_length = open_field_length / n_rows  # length of each box.

    sub_rois = []  # splits the roi into n_boxes of equal rois.
    for row in range(n_rows):
        for col in range(n_rows):
            sub_x = open_field_base[video_idx][0] + col*sub_width
            sub_y = open_field_base[video_idx][2] + row*sub_length
            sub_x_max, sub_y_max = sub_x + sub_width, sub_y + sub_length
            sub_rois.append((sub_x, sub_x_max, sub_y, sub_y_max))

    # split the sub_rois into edges and centers.
    edge_rois_temp = []
    center_rois_temp = []
    for sub_roi in sub_rois:
        # if any of the sub_rois include a dimension of the full roi, it means it's an edge!
        if any(item in sub_roi for item in open_field_base[video_idx]):
            edge_rois_temp.append(sub_roi)
        else:
            center_rois_temp.append(sub_roi)
    edge_rois[video_idx] = edge_rois_temp
    center_rois[video_idx] = center_rois_temp

    # define corner rois.
    corner_roi_idxs = [0, n_rows-1, N_BOXES-n_rows, N_BOXES-1]
    corner_rois[video_idx] = [sub_rois[index] for index in corner_roi_idxs]

# %% Anxiety-psy_behaviors analysis
# Function to define frames that animal is in various ROIs.
START_FRAME = np.round(START_SEC*FRAMERATE).astype(int)
END_FRAME = np.round(END_SEC*FRAMERATE).astype(int)


def timestamps_in_roi(roi_type, percent_time=True):
    """
    Returns timestamps when animal is in a given roi type (edge, center, corner)

    Parameters:
    - roi_type: list of list of tuples
        Each element corresponds to a video. Each video element is a list of rois of that type.
        Each roi is defined as a tuple of (x_min, x_max, y_min, y_max).
    - percent_time: bool
        If True, returns time in roi as percent of total analysis time.

    Returns:
    - frames_in_roi: list of np.arrays
        Each element corresponds to a video. Each video element is an array of frame indices
        when the animal is in any of the rois of that type.
    """
    frames_in_roi = [None]*len(paths_vid)
    time_in_roi = np.empty(len(paths_vid))
    for vid_i in range(len(paths_vid)):
        frames_in_roi_temp = []
        # cycle through all the rois in the roi type
        for roi_no in range(len(roi_type[vid_i])):
            frames_in_roi_temp.append(
                psy_beh.find_timepoints_in_coords(body_coords[vid_i],
                                                  roi_type[vid_i][roi_no]))
            frames_in_roi[vid_i] = [item for sublist in frames_in_roi_temp for item in sublist]
            frames_in_roi[vid_i] = np.sort(np.array(frames_in_roi[vid_i]))
            # trim frames in roi so that it only includes analysis frames of interest
            frames_in_roi[vid_i] = frames_in_roi[vid_i][(
                frames_in_roi[vid_i] >= START_FRAME[vid_i]) &
                (frames_in_roi[vid_i] <= END_FRAME[vid_i])]
            time_in_roi[vid_i] = len(frames_in_roi[vid_i])/FRAMERATE[vid_i]
            if percent_time:
                time_in_roi[vid_i] = time_in_roi[vid_i]/(END_SEC-START_SEC)*100

    return frames_in_roi, time_in_roi


# main analysis
frames_in_edges, time_in_edges = timestamps_in_roi(edge_rois)
frames_in_center, time_in_center = timestamps_in_roi(center_rois)
frames_in_corner, time_in_corners = timestamps_in_roi(corner_rois)

# Get data in interval chunks
n_intervals = int(round((END_SEC-START_SEC)/INTERVAL_SEC))
time_in_edges_bins = np.empty((len(paths_vid), n_intervals))
time_in_corners_bins = np.empty((len(paths_vid), n_intervals))
time_in_center_bins = np.empty((len(paths_vid), n_intervals))
for video_idx in range(len(paths_vid)):
    i = 0
    for sec_start in range(START_SEC, END_SEC, INTERVAL_SEC):
        sec_end = sec_start + INTERVAL_SEC - 1
        frame_start = int(round(sec_start*FRAMERATE[video_idx]))
        frame_end = int(round(sec_end*FRAMERATE[video_idx] + (FRAMERATE[video_idx]-1)))

        time_in_edges_bins[video_idx, i] = np.sum((frames_in_edges[video_idx] >= frame_start) &
                                                  (frames_in_edges[video_idx] <= frame_end)
                                                  )/FRAMERATE[video_idx]/INTERVAL_SEC*100
        time_in_corners_bins[video_idx, i] = np.sum((frames_in_corner[video_idx] >= frame_start) &
                                                    (frames_in_corner[video_idx] <= frame_end)
                                                    )/FRAMERATE[video_idx]/INTERVAL_SEC*100
        time_in_center_bins[video_idx, i] = np.sum((frames_in_center[video_idx] >= frame_start) &
                                                   (frames_in_center[video_idx] <= frame_end)
                                                   )/FRAMERATE[video_idx]/INTERVAL_SEC*100
        i = i+1

# %% Locomotor measures
# Preallocate variables.
distance_travelled = np.empty(len(paths_vid))
velocity = np.empty(len(paths_vid))
distance_in_center = np.empty(len(paths_vid))
time_moving = np.empty(len(paths_vid))
time_running = np.empty(len(paths_vid))
time_walking = np.empty(len(paths_vid))
time_freezing = np.empty(len(paths_vid))
velocity_while_moving = np.empty(len(paths_vid))

distance_travelled_bins = np.empty((len(paths_vid), n_intervals))
velocity_bins = np.empty((len(paths_vid), n_intervals))
distance_in_center_bins = np.empty((len(paths_vid), n_intervals))
time_moving_bins = np.empty((len(paths_vid), n_intervals))
time_running_bins = np.empty((len(paths_vid), n_intervals))
time_walking_bins = np.empty((len(paths_vid), n_intervals))
time_freezing_bins = np.empty((len(paths_vid), n_intervals))
velocity_while_moving_bins = np.empty((len(paths_vid), n_intervals))

for video_idx, vid_path in enumerate(paths_vid):
    # get distance array: distance travelled per frame.
    diff_array = np.diff(body_coords[video_idx], axis=0)
    diff_array[:, 0] = diff_array[:, 0] * cm_to_pixels_x
    diff_array[:, 1] = diff_array[:, 1] * cm_to_pixels_y
    distance_array = np.sqrt((diff_array[:, 0]**2) + (diff_array[:, 1]**2))
    if len(distance_array) < END_FRAME[video_idx]-1:
        print(
            "WARNING:", vid_path, "does not have the full analysis period present."
            " Zero-padding distance array."
            )
        distance_array = np.pad(
            distance_array,
            (0, END_FRAME[video_idx] - len(distance_array)),
            mode='constant')  # pad distance array with 0s to make divisible by framerate.

    if len(distance_array) % FRAMERATE[video_idx] != 0:
        distance_array = np.pad(
            distance_array,
            (0, int(FRAMERATE[video_idx] - len(distance_array) % FRAMERATE[video_idx])),
            mode='constant')  # pad distance array with 0s to make divisible by framerate.

    distance_array = distance_array[START_FRAME[video_idx]:END_FRAME[video_idx]+1]

    # Summary locomotion metrics
    distance_travelled[video_idx] = sum(distance_array)
    velocity[video_idx] = distance_travelled[video_idx]/(END_SEC-START_SEC)
    distance_in_center[video_idx] = sum(distance_array[frames_in_center[video_idx]-1])

    # Fine grained movement data. First, transform distance array from per frame to per second.
    n_frames = distance_array.shape[0]
    n_seconds = n_frames/FRAMERATE[video_idx]
    t_boundaries = np.arange(n_frames+1) / FRAMERATE[video_idx]
    cumulative_distance = np.concatenate(([0.0], np.cumsum(distance_array)))
    t_samples = np.arange(0, int(np.ceil(n_seconds)) + 1)
    cum_at_seconds = np.interp(t_samples, t_boundaries, cumulative_distance)
    distance_per_second = np.diff(cum_at_seconds)

    time_moving[video_idx] = sum(distance_per_second >= 5)/n_seconds*100
    time_running[video_idx] = sum(distance_per_second > 20)/n_seconds*100
    time_walking[video_idx] = time_moving[video_idx] - time_running[video_idx]
    time_freezing[video_idx] = sum(distance_per_second < 5)/n_seconds*100
    velocity_while_moving[video_idx] = np.mean(distance_per_second[distance_per_second >= 5])

    # Get all measures in interval chunks
    i = 0
    for sec_start in range(START_SEC, END_SEC, INTERVAL_SEC):
        sec_end = sec_start + INTERVAL_SEC - 1
        frame_start = int(round(sec_start*FRAMERATE[video_idx]))
        frame_end = int(round(sec_end*FRAMERATE[video_idx] + (FRAMERATE[video_idx]-1)))

        distance_travelled_bins[video_idx, i] = np.sum(distance_array[frame_start:frame_end+1])
        velocity_bins[video_idx, i] = distance_travelled_bins[video_idx, i]/(INTERVAL_SEC)
        frames_in_center_bin = (
            frames_in_center[video_idx][(frames_in_center[video_idx] >= frame_start) &
                                        (frames_in_center[video_idx] <= frame_end)])
        distance_in_center_bins[video_idx, i] = sum(distance_array[frames_in_center_bin])
        time_moving_bins[video_idx, i] = (
            sum(distance_per_second[sec_start:sec_end] >= 5)/INTERVAL_SEC*100)
        time_running_bins[video_idx, i] = (
            sum(distance_per_second[sec_start:sec_end] > 20)/INTERVAL_SEC*100)
        time_walking_bins[video_idx, i] = (
            time_moving_bins[video_idx, i] - time_running_bins[video_idx, i])
        time_freezing_bins[video_idx, i] = (
            sum(distance_per_second[sec_start:sec_end] < 5)/INTERVAL_SEC*100)
        velocity_while_moving_bins[video_idx, i] = np.mean(
            distance_per_second[sec_start:sec_end][distance_per_second[sec_start:sec_end] >= 5])
        i = i+1

# %% Create data dictionary and tables for data of interest.
# Summary data
data_summary_dict = {
    'Animal_ID': paths_vid,
    'Sex': sex_key,
    'Treatment': treatment_key,
    'Institution': INSTITUTION,
    'Time_Center': time_in_center,
    'Time_Corners': time_in_corners,
    'Distance_in_center': distance_in_center,
    'Velocity': velocity,
    'Time_moving': time_moving,
    'Velocity_while_moving': velocity_while_moving
}

data_summary_df = pd.DataFrame(data_summary_dict)
column_names = data_summary_df.columns[4:].tolist()
column_labels = ["Time in center (%)", "Time in corners (%)", "Distance in center (cm)",
                 "Velocity (cm/s)", "Time spent moving (%)", "Velocity while moving (cm/s)"]


def make_binned_data_df(data_for_df, bin_name, col_name):
    """
    Returns a dataframe with binned data in long format.
    Parameters:
    - data_for_df: np.array
        2D array where each row corresponds to an animal and each column to a time bin.
    - bin_name: str
        Name of the bin variable (e.g., "Interval").
    - col_name: str
        Name of the data variable (e.g., "Time_Center").

    Returns:
    - df_final: pd.DataFrame
        DataFrame in long format with columns for animal info, bin variable, and data variable.
    """
    df = pd.DataFrame(data_for_df, columns=[f'{i}' for i in range(data_for_df.shape[1])])
    melted_df = df.melt(var_name=bin_name, value_name=col_name)
    df_final = pd.concat([animal_info, melted_df], axis=1)
    df_final[bin_name] = df_final[bin_name].astype(float)
    return df_final


# Make data dictionary for binned data
animal_info = pd.DataFrame({"Animal_ID": paths_vid*n_intervals,
                            "Sex": sex_key*n_intervals,
                            "Treatment": treatment_key*n_intervals,
                            "Institution": INSTITUTION})

data_to_bin = [time_in_center_bins, time_in_corners_bins,
               distance_in_center_bins, velocity_bins,
               time_moving_bins, velocity_while_moving_bins]
binned_data_df = []
for col_i, binned_data in enumerate(data_to_bin):
    binned_data_df.append(make_binned_data_df(binned_data, "Interval", column_names[col_i]))
data_binned_dict = dict(zip(column_names, binned_data_df))

data_binned_dict_flat = {}
for col, data in zip(column_names, data_to_bin):
    for interval in range(n_intervals):
        col_interval_name = col + "_" + str(interval+1)
        data_binned_dict_flat[col_interval_name] = data[:, interval]

data_binned_df = pd.DataFrame(data_binned_dict_flat)

data_all = pd.concat([data_summary_df, data_binned_df], axis=1)

# %% Get stats
stats_summary_results = {}
stats_binned_results = {}
for col in column_names:
    stats_summary_results[col] = psy_stats.stats_treatment_sex(data_summary_df, col)
    stats_binned_results[col] = psy_stats.stats_treatment_sex_third(data_binned_dict[col],
                                                                    col, 'Interval')

# %% Plot all results
os.makedirs(os.path.join(FOLDER_PATH, "saved_data"), exist_ok=True)

fig, ax = plt.subplots(1, len(column_names), figsize=(1.2*len(column_names), 1.5))
for col_i, col in enumerate(column_names):
    psy_plot.plot_individual_lab(ax[col_i], data_summary_df, col, column_labels[col_i],
                                 INSTITUTION, stats_summary_results[col]["significance"])
plt.savefig(os.path.join(FOLDER_PATH, "saved_data", "summary_data_plots.png"))

fig, ax = plt.subplots(1, len(column_names), figsize=(1.5*len(column_names), 1.5))
for col_i, col in enumerate(column_names):
    psy_plot.plot_over_time(ax[col_i], data_binned_dict, col, column_labels[col_i],
                            [5+START_SEC, END_SEC/60+5])
plt.savefig(os.path.join(FOLDER_PATH, "saved_data", "time_binned_plots.png"))

# %% Plot animal traces for visualization.
for video_idx, video_path in enumerate(paths_vid):
    fig, ax = plt.subplots()
    psy_plot.plot_traces(ax, video_path, body_coords[video_idx])
    # show open field base rectangle
    psy_plot.plot_roi_coords(ax, open_field_base[video_idx])
    plt.savefig(os.path.join(FOLDER_PATH, "saved_data", f"{video_idx:02}_trace.png"))

# %% Save csv files with the relevant output data in the same directory as the script.
data_all.to_csv(os.path.join(FOLDER_PATH, "saved_data", 'OFT_data.csv'), index=False)
