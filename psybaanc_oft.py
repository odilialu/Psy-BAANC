# -*- coding: utf-8 -*-
"""
Created on Mon Jun 30 09:41:34 2025
@author: olu

OFT Analysis 
Instructions:
    1. Put all video and coordinate files from a single experiment in a folder.
    2. Adjust the variables in the "Variables to change" section to match your settings.
    3. Run the script.
    4. Update Supplementary table 1 [meta data sheet] with your raw data.
    5. Update Supplementary table 2 [stats table] with the statistics.
    6. Update your Prism files with the raw data. Ensure figure format is that of the paper.
    7. Send prism files back to Odilia. 
    
The output data is stored in: 
    data_final: all summary measures
    data_final_zscored: all summary measures zscored to sal conditions
    data_final_binned: all measures reported in 5 minute intervals. 
    levenes_results: statistical results from levene's test for homogeneity of variances
    stats_results: statistical results for every measure. 
        - If equal variances: 
            Two-way ANOVA. If interaction significant, t-test post-hoc comparisons with Sidak's corrections
        - If unequal variances: 
            One-way Kruskal-Wallis. If significant, Dunn's post hoc comparisons with Holms corrections
    levenes_results_zscored: for the zscored data
    stats_results_zscored: for the zscored data
    
The output variables include: 
    'Time_Edges': (%) time an animal spends in the edge
    'Time_Center': (%) time an animal spends in the center
    'Time_Corners': (%) time an animal spends in the corners
    'Distance_in_center': (cm) distance travelled in the center
    'Distance': (cm) distance travelled in open field
    'Velocity': (cm/s) average velocity in open field
    'Time_moving': (%) time an animal spends moving (defined as movement >= 5 cm/s)
    'Velocity_while_moving': (cm/s) average velocity during time points considered as moving
    'Time_running': (%) time an animal spends running (defined as movement > 20 cm/s)
    'Time_walking': (%) time an animal spends walking (defined as movement >= 5 or < 20 cm/s)
    'Time_freezing': (%) time an animal spends freezing (defined as movement < 5 cm/s)
    
"""

#%% Import packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cv2
import psybaanc_behavior as psy_beh
import psybaanc_stats as psy_stats

#%% Variables to change
FOLDER_PATH = r"Y:\PsyBAANC\paperExperiments\OFT_EPM\OFT_all" # path to the folder with all video and coordinate data
VIDEO_TYPE = "avi" # options: "mp4", "avi", others also likely OK.
COORDINATE_FILE_TYPE = "csv" # options: "csv", "xlsx"
START_SEC = 0 # the time in seconds that you wish to begin analysis.
END_SEC = 30*60 # the time in seconds that you wish to end analysis. 
INTERVAL_SEC = 5*60 # for time-binned measures, what is the bin interval, in seconds? 

LENGTH_CM = 50 # true size of your open field box in cm
N_BOXES = 25 # How many squares do you want to divide OF into? number must be square of an integer. PsyBAANC keep at 25. 

X_COORDINATE_INDEX = 25 # Index of your x-coordinate column in your coordinates file. Note, index starts at 0. 
Y_COORDINATE_INDEX = 26 # Index of your y-coordinate column in your coordinates file. Note, index starts at 0. 
ROW_INDEX = 4 # what row do you start to see data appear in your coordinate files? For DLC, usually 4. 
DATA_DLC = True # Is your data from deeplabcut (DLC)? true or false. If true, linear interpolation based on likelihood is done on coordinates.

sex_key = ["M"]*20 + ["F"]*20 # Create a list indicating the sex of the animals, i.e., ["M", "F", "M"]
treatment_key = ["P"]*5 + ["S"]*5 + ["P"]*5 + ["S"]*5 + ["S"]*5 + ["P"]*5 + ["S"]*5 + ["P"]*5 #Create a list indicating the treatment of the animals, i.e., ["S", "S", "P"]

zscore = True # Do you want to get z-scored data? True or False
stats = True # Do you want to return stats? True or False
plot_traces = False # Do you want to plot animal traces? True or False

# matplotlib plotting parameters
plt.rcParams['figure.dpi'] = 600
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 8

#%% Get paths of videos and csv files.
paths_coordinates = psy_beh.get_file_paths(FOLDER_PATH, COORDINATE_FILE_TYPE)
paths_vid = psy_beh.get_file_paths(FOLDER_PATH, VIDEO_TYPE)

#%% Get framerate of all your video files
FRAMERATE = np.empty((len(paths_vid)))
for video_idx, video_path in enumerate(paths_vid):
    cap = cv2.VideoCapture(video_path)
    FRAMERATE[video_idx] = (cap.get(cv2.CAP_PROP_FPS))

#%% Read in all of the coordinate data.
print("reading coordinate data. please be patient")
body_coords = [None]*len(paths_vid)
for file_idx, coordinate_file in enumerate(paths_coordinates):
    if DATA_DLC == True: # use likelihood estimates to do linear interpolation on uncertain coordinates. 
        coordinates = pd.read_csv(coordinate_file, skiprows=list(range(0, ROW_INDEX-2)))
        body_coords[file_idx] = psy_beh.get_body_parts(coordinates, body_part_idx_x = X_COORDINATE_INDEX, likelihood=0.6)

    else: # otherwise read in coordinates as is. 
        if COORDINATE_FILE_TYPE == "xlsx":
            body_coords[file_idx] = pd.read_excel(coordinate_file, 
                                                  usecols=[X_COORDINATE_INDEX,Y_COORDINATE_INDEX], 
                                                  skiprows=list(range(0,(ROW_INDEX-2))))
        elif COORDINATE_FILE_TYPE == "csv":
            body_coords[file_idx] = pd.read_csv(coordinate_file, 
                                                usecols=[X_COORDINATE_INDEX,Y_COORDINATE_INDEX], 
                                                skiprows=list(range(0,(ROW_INDEX-2))))
print("done reading coordinate data")

#%% Define all relevant ROIs. 
# Preallocate lists. 
open_field_base = [None]*(len(paths_vid))
edge_rois = [None]*(len(paths_vid))
center_rois = [None]*(len(paths_vid))
corner_rois = [None]*(len(paths_vid))

# First, loop through all the videos and draw a rectangle that defines the base of the open field if you haven't already. 
# Function saves the ROI tuple into a folder named roi in your video path so you only define it for each video once.
# If you want to re-define the ROI, you must delete the ROI pickle that is saved in the folder.
print("If not already done, draw rectangle ROI to define the base of the open field. Hit enter. Repeat for each video.")
for video_idx, video_file in enumerate(paths_vid):
    if plot_traces:
        open_field_base[video_idx] = psy_beh.get_roi_from_frame(video_file, "roi", body_coords[video_idx])
    else:
        open_field_base[video_idx] = psy_beh.get_roi_from_frame(video_file, "roi")

# Then, divide the open field base into N_BOXES number of equal portions. 
# From the boxes, define edge, center, and corner ROIs.
for video_idx, video_file in enumerate(paths_vid):
    open_field_width = open_field_base[video_idx][1] - open_field_base[video_idx][0] # x_max - x_min
    open_field_length = open_field_base[video_idx][3] - open_field_base[video_idx][2] # y_max - y_min
    cm_to_pixels_x = LENGTH_CM/open_field_width # length in centimeters of open field / size in pixels. 
    cm_to_pixels_y = LENGTH_CM/open_field_length
    
    n_rows = int(np.sqrt(N_BOXES))
    sub_width = open_field_width / n_rows # width of each box.
    sub_length = open_field_length / n_rows # length of each box.
    
    sub_rois = [] # splits the roi into n_boxes of equal rois. 
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
    
#%% Anxiety-psy_behaviors analysis
# Function to define frames that animal is in various ROIs. 
START_FRAME = (START_SEC*FRAMERATE).astype(int)
END_FRAME = (END_SEC*FRAMERATE).astype(int)
def timestamps_in_roi(roi_type, percent_time=True):
    frames_in_roi = [None]*len(paths_vid)
    time_in_roi = np.empty(len(paths_vid))
    for video_idx in range(len(paths_vid)):

        timestamps = np.array(list(range(1,len(body_coords[video_idx])+1)))
        frames_in_roi_temp = []
        # cycle through all the rois in the roi type
        for roi_no in range(len(roi_type[video_idx])):
            frames_in_roi_temp.append(psy_beh.find_time_in_location(body_coords[video_idx], timestamps, roi_type[video_idx][roi_no]))
            frames_in_roi[video_idx] = [item for sublist in frames_in_roi_temp for item in sublist]
            frames_in_roi[video_idx] = np.sort(np.array(frames_in_roi[video_idx]))
            # trim frames in roi so that it only includes analysis frames of interest 
            frames_in_roi[video_idx] = frames_in_roi[video_idx][(frames_in_roi[video_idx]>=START_FRAME[video_idx]) & (frames_in_roi[video_idx]<=END_FRAME[video_idx])]
            time_in_roi[video_idx] = len(frames_in_roi[video_idx])/FRAMERATE[video_idx] # in seconds.
            if percent_time:
                time_in_roi[video_idx] = time_in_roi[video_idx]/(END_SEC-START_SEC)*100
            
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
        
        time_in_edges_bins[video_idx, i] = np.sum((frames_in_edges[video_idx]>=frame_start) & (frames_in_edges[video_idx]<=frame_end))/FRAMERATE[video_idx]/INTERVAL_SEC*100
        time_in_corners_bins[video_idx, i] = np.sum((frames_in_corner[video_idx]>=frame_start) & (frames_in_corner[video_idx]<=frame_end))/FRAMERATE[video_idx]/INTERVAL_SEC*100
        time_in_center_bins[video_idx, i] = np.sum((frames_in_center[video_idx]>=frame_start) & (frames_in_center[video_idx]<=frame_end))/FRAMERATE[video_idx]/INTERVAL_SEC*100
        i = i+1
        
#%% Locomotor measures
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

for video_idx in range(len(paths_vid)):
    # get distance array: distance travelled per frame. 
    diff_array = np.diff(body_coords[video_idx], axis=0)
    diff_array[:, 0] = diff_array[:, 0] * cm_to_pixels_x # transform x and y differently since they may have different scales. 
    diff_array[:, 1] = diff_array[:, 1] * cm_to_pixels_y
    distance_array = np.sqrt((diff_array[:, 0]**2) + (diff_array[:, 1]**2))
    if len(distance_array) < END_FRAME[video_idx]-1:
        print("WARNING:", paths_vid[video_idx], "does not have the full analysis period present. Zero-padding distance array.")
        distance_array = np.pad(distance_array, (0, END_FRAME[video_idx] - len(distance_array)), mode='constant') # pad distance array with 0s to make divisible by 30. 

    if len(distance_array) % FRAMERATE[video_idx] != 0:
        distance_array = np.pad(distance_array, (0, int(FRAMERATE[video_idx] - len(distance_array) % FRAMERATE[video_idx])), mode='constant') # pad distance array with 0s to make divisible by 30. 

    distance_array = distance_array[START_FRAME[video_idx]:END_FRAME[video_idx]]
    
    # Summary locomotion metrics
    distance_travelled[video_idx] = sum(distance_array)
    velocity[video_idx] = distance_travelled[video_idx]/(END_SEC-START_SEC)
    distance_in_center[video_idx] = sum(distance_array[frames_in_center[video_idx]-1])
    
    # Fine grained movement data
    n_seconds = END_SEC - START_SEC
    distance_per_second = distance_array.reshape(n_seconds, int(FRAMERATE[video_idx])).sum(axis=1)
    time_moving[video_idx] = sum(distance_per_second>=5)/n_seconds*100
    time_running[video_idx] = sum(distance_per_second>20)/n_seconds*100
    time_walking[video_idx] = time_moving[video_idx] - time_running[video_idx]
    time_freezing[video_idx] = sum(distance_per_second<5)/n_seconds*100
    velocity_while_moving[video_idx] = np.mean(distance_per_second[distance_per_second>=5])
    
    # Get all measures in interval chunks
    i = 0
    for sec_start in range(START_SEC, END_SEC, INTERVAL_SEC):
        sec_end = sec_start + INTERVAL_SEC - 1
        frame_start = int(round(sec_start*FRAMERATE[video_idx]))
        frame_end = int(round(sec_end*FRAMERATE[video_idx] + (FRAMERATE[video_idx]-1)))

        distance_travelled_bins[video_idx, i] = np.sum(distance_array[frame_start:frame_end+1])
        velocity_bins[video_idx, i] = distance_travelled_bins[video_idx, i]/(INTERVAL_SEC)
        frames_in_center_bin = (frames_in_center[video_idx][(frames_in_center[video_idx] >=frame_start) & (frames_in_center[video_idx] <=frame_end)])-1
        distance_in_center_bins[video_idx, i] = sum(distance_array[frames_in_center_bin])
        time_moving_bins[video_idx, i] = sum(distance_per_second[sec_start:sec_end]>=5)/INTERVAL_SEC*100
        time_running_bins[video_idx, i] = sum(distance_per_second[sec_start:sec_end]>20)/INTERVAL_SEC*100
        time_walking_bins[video_idx, i] = time_moving_bins[video_idx, i] - time_running_bins[video_idx, i]
        time_freezing_bins[video_idx, i] = sum(distance_per_second[sec_start:sec_end]<5)/INTERVAL_SEC*100
        velocity_while_moving_bins[video_idx, i] = np.mean(distance_per_second[sec_start:sec_end][distance_per_second[sec_start:sec_end]>=5])
        i = i+1
        
#%% Create data dictionary for summary data. 
data_dict = {
    'Animal_ID': paths_vid,
    'Sex': sex_key,
    'Treatment': treatment_key,
    'Time_Edges': time_in_edges,
    'Time_Center': time_in_center,
    'Time_Corners': time_in_corners,
    'Distance_in_center': distance_in_center,
    'Distance': distance_travelled,
    'Velocity': velocity,
    'Time_moving': time_moving,
    'Velocity_while_moving': velocity_while_moving,
    'Time_running': time_running,
    'Time_walking': time_walking,
    'Time_freezing': time_freezing
    }

data_final = pd.DataFrame(data_dict)
column_names = data_final.columns[3:].tolist()

#%% Make data dictionary for binned data 
animal_info = pd.DataFrame({"Animal_ID": paths_vid, "Sex": sex_key, "Treatment": treatment_key})
def make_binned_data_df(binned_data):
    df = pd.DataFrame(binned_data)
    df_final = pd.concat([animal_info, df], axis=1)
    return df_final

binned_data_df = []
for binned_data in [time_in_edges_bins, time_in_center_bins, time_in_corners_bins,
                    distance_in_center_bins, distance_travelled_bins, velocity_bins,
                    time_moving_bins, velocity_while_moving_bins,
                    time_running_bins, time_walking_bins, time_freezing_bins]:
    binned_data_df.append(make_binned_data_df(binned_data))

data_final_binned = dict(zip(column_names, binned_data_df))
    
#%% Z-scored data of summary data
if zscore:
    data_final_zscored = psy_beh.zscore_dataframe(data_final)

#%% Do statistical tests.
def get_stats(data_final):
    levenes = psy_stats.levenes_test_dataframe(data_final)
    
    results = [None]*len(column_names)
    for col_idx, col_name in enumerate(column_names):
        if levenes[levenes["Measure"]==col_name]["P_value"].tolist()[0] < 0.05:
            kruskal_results = np.array(psy_stats.kruskal_wallis(data_final, col_name)).T
            keys = ["h_stat", "p-value"]
            if kruskal_results[1] < 0.05:
                pvals_holm = psy_stats.dunns(data_final, col_name)
                kruskal_results = np.concatenate((kruskal_results, pvals_holm)).reshape(1, -1)
                keys = ["h-stat", "p-value", "m,sal v. m,psi", "f,sal v. f,psi", "m,sal v. f,sal", "m,psi v. f,psi"]
                
            results[col_idx] = pd.DataFrame(data = kruskal_results, columns = keys)
                
        else:
            formula = col_name + ' ~ Treatment*Sex'
            results[col_idx] = (psy_stats.sm_ANOVA(formula=formula, data=data_final))
            if results[col_idx].loc["Treatment:Sex", "PR(>F)"] < 0.05:
                results_sidaks = psy_stats.sidaks(data_final, col_name)
                results[col_idx] = pd.concat([results[col_idx], results_sidaks])
                
    stats_results =dict(zip(column_names, results))
    levenes_results = levenes.set_index('Measure').to_dict(orient='index')
    
    return levenes_results, stats_results

if stats:
    levenes_results, stats_results = get_stats(data_final)
    levenes_results_zscored, stats_results_zscored = get_stats(data_final_zscored)
    
#%% Plot animal traces for visualization if you like. 
if plot_traces:
    for video_idx, video_path in enumerate(paths_vid):
        psy_beh.plot_traces(video_path, body_coords[video_idx])
