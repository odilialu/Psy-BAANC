# -*- coding: utf-8 -*-
"""
Created on Sat Jul  5 18:42:30 2025

@author: olu

SIT Analysis:
Instructions:
    1. Put all video and coordinate files from a single experiment in a folder.
    2. Adjust the variables in the "Variables to change" section to match your settings.
    3. Run the script. Follow command-line instructions when prompted.
    4. Update Supplementary table 1 [meta data sheet] with your raw data.
    5. Update Supplementary table 2 [stats table] with the statistics.
    6. Update your Prism files with the raw data. Ensure figure format is that of the paper.
    7. Upload your Prism file to the google drive, and notify Odilia when complete.

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
    'Empty_cup_time': time spent exploring empty cup (s),
    'Social_cup_time': time spent exploring social cup (s),
    'Social_index': social preference score = (social time - empty time)/(social time + empty time),
    'Latency_empty': latency to visit empty cup (s),
    'Latency_social': latency to visit social cup (s),
    'Visits_empty': number of visits to empty cup,
    'Visits_social': number of visits to social cup,
    'Empty_chamber_time': time in empty chamber (%),
    'Social_chamber_time': time in social chamber (%),
    'Center_chamber_time': time in center chamber (%),
    'Distance_empty_chamber': distance moved in empty chamber (cm),
    'Distance_social_chamber': distance moved in social chamber (cm),
    'Distance_travelled': distance travelled total (cm),
    'Velocity': average velocity (cm/s)

"""

#%% Import packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cv2
import psybaanc_behavior as psy_beh
import psybaanc_stats as psy_stats

#%% Variables to change
FOLDER_PATH = r"Y:/PsyBAANC/paperExperiments/SM_SIT/SIT_Pers_all" # path to the folder with all video and coordinate data
VIDEO_TYPE = "avi" # options: "mp4", "avi", others also likely OK.
COORDINATE_FILE_TYPE = "csv" # options: "csv", "xlsx"
START_SEC = 0 # the time in seconds that you wish to begin analysis.
END_SEC = 10*60 # the time in seconds that you wish to end analysis. 

LENGTH_CM = 40 # true size of your open field box in cm

X_NOSE_COORDINATE_INDEX = 1
Y_NOSE_COORDINATE_INDEX = 2
X_BODY_COORDINATE_INDEX = 4 # Index of your x-coordinate column in your coordinates file. Note, index starts at 0. 
Y_BODY_COORDINATE_INDEX = 5 # Index of your y-coordinate column in your coordinates file. Note, index starts at 0. 
ROW_INDEX = 4 # what row do you start to see data appear in your coordinate files? For DLC, usually 4. 
DATA_DLC = True # Is your data from deeplabcut (DLC)? true or false. If true, linear interpolation based on likelihood is done on coordinates.
COORDINATES_CM = False # Are your coordinates in centimeters? (And not pixels)

OBJECT_ONE_SHAPE = "ellipse" # options: "circle", "ellipse", "rectangle", "polygon"
OBJECT_TWO_SHAPE = "ellipse" # options: "circle", "ellipse", "rectangle", "polygon"
sex_key = ["M"]*16 + ["F"]*16 # Create a list indicating the sex of the animals, i.e., ["M", "F", "M"]
treatment_key = ["P"]*5 + ["S"]*5 + ["P"]*3 + ["S"]*3 + ["P"]*5 + ["S"]*8 + ["P"]*3

zscore = True # Do you want to get z-scored data? True or False
stats = True # Do you want to return stats? True or False
plot_traces = True # Do you want to plot animal traces? True or False

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
nose_coords = [None]*len(paths_vid)
cm_to_pixels = np.empty(len(paths_vid))
for file_idx, coordinate_file in enumerate(paths_coordinates):
    if DATA_DLC == True: # use likelihood estimates to do linear interpolation on uncertain coordinates. 
        coordinates = pd.read_csv(coordinate_file, skiprows=list(range(0, ROW_INDEX-2)))
        body_coords[file_idx] = psy_beh.get_body_parts(coordinates, body_part_idx_x = X_BODY_COORDINATE_INDEX, likelihood=0.6)
        nose_coords[file_idx] = psy_beh.get_body_parts(coordinates, body_part_idx_x = X_NOSE_COORDINATE_INDEX, likelihood=0.6)
        
    else: # otherwise read in coordinates as is. 
        if COORDINATE_FILE_TYPE == "xlsx":
            body_coords[file_idx] = pd.read_excel(coordinate_file, 
                                                  usecols=[X_BODY_COORDINATE_INDEX,Y_BODY_COORDINATE_INDEX], 
                                                  skiprows=list(range(0,(ROW_INDEX-2))))
            nose_coords[file_idx] = pd.read_excel(coordinate_file, 
                                                  usecols=[X_NOSE_COORDINATE_INDEX,Y_NOSE_COORDINATE_INDEX], 
                                                  skiprows=list(range(0,(ROW_INDEX-2))))
        elif COORDINATE_FILE_TYPE == "csv":
            body_coords[file_idx] = pd.read_csv(coordinate_file, 
                                                usecols=[X_BODY_COORDINATE_INDEX,Y_BODY_COORDINATE_INDEX], 
                                                skiprows=list(range(0,(ROW_INDEX-2))))
            nose_coords[file_idx] = pd.read_csv(coordinate_file, 
                                                  usecols=[X_NOSE_COORDINATE_INDEX,Y_NOSE_COORDINATE_INDEX], 
                                                  skiprows=list(range(0,(ROW_INDEX-2))))
            
        body_coords[file_idx] = body_coords[0].to_numpy().astype(float)
        
    if COORDINATES_CM:
        cm_to_pixels[file_idx] = psy_beh.calibrate_pixels_to_cm(paths_vid[file_idx], real_world_cm=LENGTH_CM, frame_number=0)
        body_coords[file_idx] = body_coords[file_idx]/cm_to_pixels[file_idx]

print("done reading coordinate data")

#%% Calibrate pixels to cm
cm_to_pixels = np.empty(len(paths_vid))
for video_idx, video_path in enumerate(paths_vid):
    cm_to_pixels[video_idx] = psy_beh.calibrate_pixels_to_cm(video_path, real_world_cm=LENGTH_CM, frame_number=0)

#%% Define all relevant ROIs. 
# Preallocate lists. 
empty_chamber = [None]*(len(paths_vid))
social_chamber = [None]*(len(paths_vid))

empty_cup = [None]*len(paths_vid)
social_cup = [None]*len(paths_vid)
empty_cup_roi = [None]*len(paths_vid)
social_cup_roi = [None]*len(paths_vid)

print("For each video, you will define in this order: 1) empty chamber 2) social chamber 3) empty cup and 4) social cup.")
for video_idx, video_file in enumerate(paths_vid):
    empty_chamber[video_idx] = psy_beh.get_roi_from_frame(video_file, "empty_chamber_base")
    social_chamber[video_idx] = psy_beh.get_roi_from_frame(video_file, "social_chamber_base")
    empty_cup[video_idx] = psy_beh.get_roi_flexible(video_file, "empty_cup", shape=OBJECT_ONE_SHAPE)
    social_cup[video_idx] = psy_beh.get_roi_flexible(video_file, "social_cup", shape=OBJECT_TWO_SHAPE)
    empty_cup_roi[video_idx] = psy_beh.resize_mask_by_cm(empty_cup[video_idx], cm_offset=3, cm_per_pixel = cm_to_pixels[video_idx])
    social_cup_roi[video_idx] = psy_beh.resize_mask_by_cm(social_cup[video_idx], cm_offset=3, cm_per_pixel = cm_to_pixels[video_idx])

    
#%%
# Object exploration analysis
# Exploration is defined as: 
# - nose within 3 cm of cup, body is not within cup mask. (animal cannot be sitting on object)
# - animal must also be looking at object. 
def get_exploration_metrics(empty_cup_roi, empty_cup):
    empty_cup_time = np.empty(len(paths_vid))
    empty_cup_latency = np.empty(len(paths_vid))
    empty_cup_visits = np.empty(len(paths_vid))
    for video_idx in range(len(paths_vid)):
        START_FRAME = START_SEC*FRAMERATE[video_idx]
        END_FRAME = END_SEC*FRAMERATE[video_idx]
        empty_cup_nose_frames = psy_beh.get_timepoints_in_mask(nose_coords[video_idx], empty_cup_roi[video_idx])
        empty_cup_body_frames = psy_beh.get_timepoints_in_mask(body_coords[video_idx], empty_cup[video_idx])
        common_frames = np.intersect1d(empty_cup_nose_frames, empty_cup_body_frames)
        empty_cup_exploration = np.setdiff1d(empty_cup_nose_frames, common_frames)
        
        # Animal must also be looking towards the object. 
        M = cv2.moments(empty_cup[video_idx])
        cx = M["m10"] / M["m00"]
        cy = M["m01"] / M["m00"]
        object_coord = np.array([cx, cy]).reshape(1, 2)
        time_looking_empty_cup = psy_beh.timepoints_looking_at_object(nose_coords[video_idx], body_coords[video_idx], object_coord, angle_thresh_deg=30)
        
        empty_cup_frames = np.intersect1d(empty_cup_exploration, time_looking_empty_cup)
        empty_cup_frames = empty_cup_frames[(empty_cup_frames>=START_FRAME) & (empty_cup_frames<=END_FRAME)] # filter to only include analysis time of interest
        
        empty_cup_time[video_idx] = len(empty_cup_frames)/FRAMERATE[video_idx]
        
        empty_cup_frames_visits = psy_beh.split_into_visits(empty_cup_frames, min_length = FRAMERATE[video_idx]/6)
        empty_cup_visits[video_idx] = len(empty_cup_frames_visits)
        if empty_cup_visits[video_idx] > 0:
            empty_cup_latency[video_idx] = empty_cup_frames_visits[0][0]/FRAMERATE[video_idx]
        else:
            empty_cup_latency[video_idx] = END_SEC-START_SEC
        
    return empty_cup_time, empty_cup_latency, empty_cup_visits

empty_cup_time, empty_cup_latency, empty_cup_visits = get_exploration_metrics(empty_cup_roi, empty_cup)
social_cup_time, social_cup_latency, social_cup_visits = get_exploration_metrics(social_cup_roi, social_cup)
    
total_time = empty_cup_time + social_cup_time
social_index = (social_cup_time-empty_cup_time)/total_time

#%% Time spent in chambers.
length_experiment = END_SEC-START_SEC

frames_empty = [None]*len(paths_vid)
frames_social = [None]*len(paths_vid)

time_empty = np.empty(len(paths_vid))
time_social = np.empty(len(paths_vid))
time_center = np.empty(len(paths_vid))

distance_empty = np.empty(len(paths_vid))
distance_social = np.empty(len(paths_vid))

for video_idx, video_path in enumerate(paths_vid):
    timestamps = np.array(list(range(1,len(body_coords[video_idx])+1)))
    START_FRAME = int(round(START_SEC*FRAMERATE[video_idx]))
    END_FRAME = int(round(END_SEC*FRAMERATE[video_idx]))
    frames_empty[video_idx] = psy_beh.find_time_in_location(body_coords[video_idx], timestamps, empty_chamber[video_idx])
    frames_empty[video_idx] = frames_empty[video_idx][(frames_empty[video_idx]>=START_FRAME) & (frames_empty[video_idx] <= END_FRAME)]
    time_empty[video_idx] = len(frames_empty[video_idx])/FRAMERATE[video_idx]
    
    frames_social[video_idx] = psy_beh.find_time_in_location(body_coords[video_idx], timestamps, social_chamber[video_idx])
    frames_social[video_idx] = frames_social[video_idx][(frames_social[video_idx]>=START_FRAME) & (frames_social[video_idx] <= END_FRAME)]
    time_social[video_idx] = len(frames_social[video_idx])/FRAMERATE[video_idx]
    
    time_center[video_idx] = length_experiment - time_empty[video_idx] - time_social[video_idx]
    
time_center = time_center/length_experiment*100
time_empty = time_empty/length_experiment*100
time_social = time_social/length_experiment*100

#%% Distance travelled and distance travelled in each chamber.
# Preallocate variables. 
distance_travelled = np.empty(len(paths_vid))
velocity = np.empty(len(paths_vid))

distance_empty = np.empty(len(paths_vid))
distance_social = np.empty(len(paths_vid))

for video_idx in range(len(paths_vid)):
    START_FRAME = int(round(START_SEC*FRAMERATE[video_idx]))
    END_FRAME = int(round(END_SEC*FRAMERATE[video_idx]))
    # get distance array: distance travelled per frame. 
    diff_array = np.diff(body_coords[video_idx], axis=0)
    diff_array[:, 0] = diff_array[:, 0] * cm_to_pixels[video_idx] # transform x and y differently since they may have different scales. 
    diff_array[:, 1] = diff_array[:, 1] * cm_to_pixels[video_idx]
    distance_array = np.sqrt((diff_array[:, 0]**2) + (diff_array[:, 1]**2))
    if len(distance_array) < END_FRAME-1:
        print("WARNING:", paths_vid[video_idx], "does not have the full analysis period present. Zero-padding distance array.")
        distance_array = np.pad(distance_array, (0, END_FRAME - len(distance_array)), mode='constant') # pad distance array with 0s to make divisible by 30. 

    if len(distance_array) % FRAMERATE[video_idx] != 0:
        distance_array = np.pad(distance_array, (0, int(FRAMERATE[video_idx] - len(distance_array) % FRAMERATE[video_idx])), mode='constant') # pad distance array with 0s to make divisible by 30. 

    distance_array = distance_array[START_FRAME:END_FRAME]

    # Summary locomotion metrics
    distance_travelled[video_idx] = sum(distance_array)
    velocity[video_idx] = distance_travelled[video_idx]/(END_SEC-START_SEC)
    
    distance_array_empty = distance_array[(frames_empty[video_idx])-1]
    distance_empty[video_idx] = sum(distance_array_empty)
    
    distance_array_social = distance_array[(frames_social[video_idx])-1]
    distance_social[video_idx] = sum(distance_array_social)

#%% Create data dictionary for summary data. 
data_dict = {
    'Animal_ID': paths_vid,
    'Sex': sex_key,
    'Treatment': treatment_key,
    'Empty_cup_time': empty_cup_time,
    'Social_cup_time': social_cup_time,
    'Social_index': social_index,
    'Latency_empty': empty_cup_latency,
    'Latency_social': social_cup_latency,
    'Visits_empty': empty_cup_visits,
    'Visits_social': social_cup_visits,
    'Empty_chamber_time': time_empty,
    'Social_chamber_time': time_social,
    'Center_chamber_time': time_center,
    'Distance_empty_chamber': distance_empty,
    'Distance_social_chamber': distance_social,
    'Distance_travelled': distance_travelled,
    'Velocity': velocity
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
        psy_beh.plot_traces(video_path, nose_coords[video_idx])
