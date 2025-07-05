# -*- coding: utf-8 -*-
"""
Created on Wed Jul  2 16:01:40 2025

@author: olu

EPM Analysis 
Instructions:
    1. Put all video and coordinate files from a single experiment in a folder.
    2. Adjust the variables in the "Variables to change" section to match your settings.
    3. Run the script. Follow command-line instructions when prompted. 
    4. Update Supplementary table 1 [meta data sheet] with your raw data.
    5. Update Supplementary table 2 [stats table] with the statistics.
    6. Update your Prism files with the raw data. Ensure figure format is that of the paper.
    7. Send prism files back to Odilia. 
    
The output data is stored in: 
    data_final: all summary measures
    data_final_zscored: all summary measures zscored to sal conditions
    levenes_results: statistical results from levene's test for homogeneity of variances
    stats_results: statistical results for every measure. 
        - If equal variances: 
            Two-way ANOVA. If interaction significant, t-test post-hoc comparisons with Sidak's corrections
        - If unequal variances: 
            One-way Kruskal-Wallis. If significant, Dunn's post hoc comparisons with Holms corrections
    
The output variables include: 
    'Time_Open': (%) time in open arms
    'Time_Center': (%) time in center of epm
    'Latency_Open': (%) time it took for mouse to first enter open arm
    'Distance': (cm) distance travelled in epm
  
"""

#%% Import packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cv2
import psybaanc_behavior as psy_beh
import psybaanc_stats as psy_stats

#%% Variables to change
FOLDER_PATH = r"Y:/PsyBAANC/paperExperiments/EPMPers" # path to the folder with all video and coordinate data
VIDEO_TYPE = "avi" # options: "mp4", "avi", others also likely OK.
COORDINATE_FILE_TYPE = "csv" # options: "csv", "xlsx"
START_SEC = 0 # the time in seconds that you wish to begin analysis.
END_SEC = 10*60 # the time in seconds that you wish to end analysis. 

LENGTH_CM = 75.5 # true size (cm) of one dimension of your epm (used for pixel to cm calibration)

X_COORDINATE_INDEX = 13 # Index of your x-coordinate column in your coordinates file. Note, index starts at 0. 
Y_COORDINATE_INDEX = 14 # Index of your y-coordinate column in your coordinates file. Note, index starts at 0. 
ROW_INDEX = 4 # what row do you start to see data appear in your coordinate files? For DLC, usually 4. 
DATA_DLC = True # Is your data from deeplabcut (DLC)? true or false. If true, linear interpolation based on likelihood is done on coordinates.

sex_key = ["M"]*27 + ["F"]*30 # Create a list indicating the sex of the animals, i.e., ["M", "F", "M"]
treatment_key = (["S", "S"] + ["P"]*5 + ["S"]*5 + ["P", "P"] + ["S"]*3 + ["P"]*6 + ["S"]*2 + ["P", "P"]
                 + ["P"]*3 + ["S"]*3 + ["P"]*3 + ["S"]*2 + ["P"]*2 + ["S"]*5 + ["P"]*3 + ["S"]*3 + ["P"]*2 + ["S"]*2 + ["P"]*2)

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
epm_roi_open_arm = [None]*(len(paths_vid))
epm_roi_closed_arm = [None]*(len(paths_vid))
epm_roi_center = [None]*(len(paths_vid))

# First, loop through all the videos and draw a polygon that defines the roi of interest. 
# Function saves the ROI mask into a folder with the roi name in your video path so you only define it for each video once.
# If you want to re-define the ROI, you must delete the ROI pickle that is saved in the folder.
print("If not already done, draw polygon ROIs to define the open and closed arms of the epm.")
for video_idx, video_file in enumerate(paths_vid):
    epm_roi_open_arm[video_idx] = psy_beh.get_roi_flexible(video_file, "open_arm", "polygon", coordinates = None)
    epm_roi_closed_arm[video_idx] = psy_beh.get_roi_flexible(video_file, "closed_arm", "polygon", coordinates = None)
    
    # define overlap between closed and open arms as center. Remove center from open and closed arm masks. 
    epm_roi_center[video_idx] = cv2.bitwise_and(epm_roi_open_arm[video_idx], epm_roi_closed_arm[video_idx])
    epm_roi_open_arm[video_idx] = cv2.subtract(epm_roi_open_arm[video_idx], epm_roi_center[video_idx])
    epm_roi_closed_arm[video_idx] = cv2.subtract(epm_roi_closed_arm[video_idx], epm_roi_center[video_idx])
    
#%% Main analysis: time in open arms (%), time in center (%), number of open arm entries, latency to enter open arm, distance travelled
length_experiment = END_SEC-START_SEC

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
    END_FRAME = int(round(END_SEC*FRAMERATE[video_idx]))
    frames_open[video_idx] = psy_beh.get_timepoints_in_mask(body_coords[video_idx], epm_roi_open_arm[video_idx])
    frames_open[video_idx] = frames_open[video_idx][(frames_open[video_idx]>=START_FRAME) & (frames_open[video_idx] <= END_FRAME)]
    frames_open_visits[video_idx] = psy_beh.split_into_visits(frames_open[video_idx], min_length = (FRAMERATE[video_idx]/6))
    
    time_open[video_idx] = len(frames_open[video_idx])/FRAMERATE[video_idx]/length_experiment*100
    visits_open[video_idx] = len(frames_open_visits[video_idx])
    latency_open[video_idx] = frames_open_visits[video_idx][0][0]/FRAMERATE[video_idx] # in seconds.
    
    frames_center = psy_beh.get_timepoints_in_mask(body_coords[video_idx], epm_roi_center[video_idx])
    frames_center = frames_center[(frames_center>=START_FRAME) & (frames_center<=END_FRAME)]
    time_center[video_idx] = len(frames_center)/FRAMERATE[video_idx]/length_experiment*100
    
    # Get distance travelled
    cm_to_pixels[video_idx] = psy_beh.calibrate_pixels_to_cm(video_path, real_world_cm=LENGTH_CM, frame_number=0)
    
    # get distance array: distance travelled per frame. 
    diff_array = np.diff(body_coords[video_idx], axis=0)
    diff_array = diff_array * cm_to_pixels[video_idx]
    distance_array = np.sqrt((diff_array[:, 0]**2) + (diff_array[:, 1]**2))
    distance_array = distance_array[START_FRAME:END_FRAME]
    
    # Summary locomotion metrics
    distance_travelled[video_idx] = sum(distance_array)
    velocity[video_idx] = distance_travelled[video_idx]/length_experiment

#%% Create data dictionary for summary data. 
data_dict = {
    'Animal_ID': paths_vid,
    'Sex': sex_key,
    'Treatment': treatment_key,
    'Time_Open': time_open,
    'Time_Center': time_center,
    'Latency_Open': latency_open,
    'Distance': distance_travelled,
    }

data_final = pd.DataFrame(data_dict)
column_names = data_final.columns[3:].tolist()

#%% Z-scored data of summary data
if zscore:
    data_final_zscored = psy_beh.zscore_dataframe(data_final)

#%% Do statistical tests.
def get_stats(data_final):
    levenes = psy_stats.levenes_test_dataframe(data_final)
    print(psy_stats.sample_size(data_final))
    
    results = [None]*len(column_names)
    for col_idx, col_name in enumerate(column_names):
        if levenes[levenes["Measure"]==col_name]["P_value"].tolist()[0] < 0.05:
            kruskal_results = np.array(psy_stats.kruskal_wallis(data_final, col_name)).reshape(1, -1)
            keys = ["h_stat", "p-value"]
            if kruskal_results[0, 1] < 0.05:
                pvals_holm = psy_stats.dunns(data_final, col_name).reshape(1, -1)
                kruskal_results = np.concatenate((kruskal_results, pvals_holm), axis=1)
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
    
#%% Plot animal traces for visualization if you like. 
if plot_traces:
    for video_idx, video_path in enumerate(paths_vid):
        psy_beh.plot_traces(video_path, body_coords[video_idx])
    
    