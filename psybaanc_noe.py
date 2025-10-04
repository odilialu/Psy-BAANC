# -*- coding: utf-8 -*-
"""
Created on Sat Jul  5 14:17:43 2025 / Edited by cliu to export csv files

@author: olu

NOE Analysis: 
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
    'Avg_time': (s) sum of exploration time of both objects
    'Latency_explore': (s) time it took for mouse to first explore an object
    'Object_one_time': (s) exploration time of object 1
    'Object_two_time': (s) exploration time of object 
    'Proportion_time': (fraction) = object_one_time/object_two_time
    'Time_edges': (%) Time animal spent in the edges
    'Time_corners': (%) Time animal spent in the corners
    'Distance': (cm) distance travelled in open field
    'Velocity': (cm/s) average velocity in open field

"""

#%% Import packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cv2
import psybaanc_behavior as psy_beh
import psybaanc_stats as psy_stats

#%% Variables to change
FOLDER_PATH = r"C:\Users\sohal\Box\Sohal lab psilocybin data\Revisions\Revision Analyses\pNOE\videos and dlc data" # path to the folder with all video and coordinate data
VIDEO_TYPE = "avi" # options: "mp4", "avi", others also likely OK.
COORDINATE_FILE_TYPE = "csv" # options: "csv", "xlsx"
START_SEC = 0 # the time in seconds that you wish to begin analysis.
END_SEC = 10*60 # the time in seconds that you wish to end analysis. 

LENGTH_CM = 49.53 # true size of your open field box in cm
N_BOXES = 25 # How many squares do you want to divide OF into? number must be square of an integer. PsyBAANC keep at 25. 

X_NOSE_COORDINATE_INDEX = 1
Y_NOSE_COORDINATE_INDEX = 2
X_BODY_COORDINATE_INDEX = 4# Index of your x-coordinate column in your coordinates file. Note, index starts at 0.
Y_BODY_COORDINATE_INDEX = 5# Index of your y-coordinate column in your coordinates file. Note, index starts at 0.
ROW_INDEX = 4 # what row do you start to see data appear in your coordinate files? For DLC, usually 4. 
DATA_DLC = True # Is your data from deeplabcut (DLC)? true or false. If true, linear interpolation based on likelihood is done on coordinates.
COORDINATES_CM = False # Are your coordinates in centimeters? (And not pixels)

OBJECT_ONE_SHAPE = "circle" # options: "circle", "ellipse", "rectangle", "polygon"
OBJECT_TWO_SHAPE = "circle" # options: "circle", "ellipse", "rectangle", "polygon"
sex_key = ["F"]*5 + ["M"]*5 + ["F"]*5 + ["M"]*5 + ["F"]*5 + ["M"]*5 + ["F"]*5 + ["M"]*5 # Create a list indicating the sex of the animals, i.e., ["M", "F", "M"]
treatment_key = ["P","S","S","P","P","P","S","S","P","P","S","S","P","S","P","S","S","P","S","P","S","P","P","S","S","S","P","P","S","S","P","P","S","S","P","P","P","S","S","P"]

zscore = True # Do you want to get z-scored data? True or False
stats = True # Do you want to return stats? True or False
plot_traces = False # Do you want to plot animal traces? True or False
export_to_csv = True # Do you want to export csv files of relevant data outputs? True or False

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
            
        body_coords[file_idx].replace('-', np.nan, inplace=True) 
        body_coords[file_idx] = body_coords[file_idx].interpolate()
        body_coords[file_idx] = body_coords[file_idx].to_numpy().astype(float)
        
        nose_coords[file_idx].replace('-', np.nan, inplace=True) 
        nose_coords[file_idx] = nose_coords[file_idx].interpolate()
        nose_coords[file_idx] = nose_coords[file_idx].to_numpy().astype(float)
        
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
open_field_base = [None]*(len(paths_vid))
edge_rois = [None]*(len(paths_vid))
center_rois = [None]*(len(paths_vid))
corner_rois = [None]*(len(paths_vid))

object_one = [None]*len(paths_vid)
object_two = [None]*len(paths_vid)
object_one_roi = [None]*len(paths_vid)
object_two_roi = [None]*len(paths_vid)

# First, loop through all the videos and draw a rectangle that defines the base of the open field if you haven't already. 
# Function saves the ROI tuple into a folder named roi in your video path so you only define it for each video once.
# If you want to re-define the ROI, you must delete the ROI pickle that is saved in the folder.
print("If not already done, for each video, draw ROIs to define the left object, then the right object.")
for video_idx, video_file in enumerate(paths_vid):
    open_field_base[video_idx] = psy_beh.get_roi_from_frame(video_file, "open_field_base", coordinates = body_coords[video_idx])

    object_one[video_idx] = psy_beh.get_roi_flexible(video_file, "object_one", shape=OBJECT_ONE_SHAPE)
    object_two[video_idx] = psy_beh.get_roi_flexible(video_file, "object_two", shape=OBJECT_TWO_SHAPE)
    object_one_roi[video_idx] = psy_beh.resize_mask_by_cm(object_one[video_idx], cm_offset=3, cm_per_pixel = cm_to_pixels[video_idx])
    object_two_roi[video_idx] = psy_beh.resize_mask_by_cm(object_two[video_idx], cm_offset=3, cm_per_pixel = cm_to_pixels[video_idx])

    
#%%
# Object exploration analysis
# Exploration is defined as: 
# - nose within 3 cm of object, body is not within object mask. (animal cannot be sitting on object)
# - animal must also be looking at object. 
START_FRAME = np.round(START_SEC*FRAMERATE).astype(int)
END_FRAME = np.round(END_SEC*FRAMERATE).astype(int)
def get_exploration_metrics(object_one_roi, object_one):
    object_one_time = np.empty(len(paths_vid))
    object_one_latency = np.empty(len(paths_vid))
    object_one_visits = np.empty(len(paths_vid))
    for video_idx in range(len(paths_vid)):

        object_one_nose_frames = psy_beh.get_timepoints_in_mask(nose_coords[video_idx], object_one_roi[video_idx])
        object_one_body_frames = psy_beh.get_timepoints_in_mask(body_coords[video_idx], object_one[video_idx])
        common_frames = np.intersect1d(object_one_nose_frames, object_one_body_frames)
        object_one_exploration = np.setdiff1d(object_one_nose_frames, common_frames)
        
        # Animal must also be looking towards the object. 
        M = cv2.moments(object_one[video_idx])
        cx = M["m10"] / M["m00"]
        cy = M["m01"] / M["m00"]
        object_coord = np.array([cx, cy]).reshape(1, 2)
        time_looking_object_one = psy_beh.timepoints_looking_at_object(nose_coords[video_idx], body_coords[video_idx], object_coord, angle_thresh_deg=30)
        
        object_one_frames = np.intersect1d(object_one_exploration, time_looking_object_one)
        object_one_frames = object_one_frames[(object_one_frames>=START_FRAME[video_idx]) & (object_one_frames<=END_FRAME[video_idx])] # filter to only include analysis time of interest
        
        object_one_time[video_idx] = len(object_one_frames)/FRAMERATE[video_idx]
        
        object_one_frames_visits = psy_beh.split_into_visits(object_one_frames, min_length = int(FRAMERATE[video_idx]/6))
        object_one_visits[video_idx] = len(object_one_frames_visits)
        if object_one_visits[video_idx] > 0:
            object_one_latency[video_idx] = object_one_frames_visits[0][0]/FRAMERATE[video_idx]
        else:
            object_one_latency[video_idx] = END_SEC-START_SEC
        
    return object_one_time, object_one_latency, object_one_visits

object_one_time, object_one_latency, object_one_visits = get_exploration_metrics(object_one_roi, object_one)
object_two_time, object_two_latency, object_two_visits = get_exploration_metrics(object_two_roi, object_two)
    
total_time = object_one_time + object_two_time
average_time = total_time/2
total_visits = object_one_visits + object_two_visits
latency_exploration = np.min(np.vstack((object_one_latency, object_two_latency)).T, axis=1)

#%% Anxiety-psy_behaviors analysis
# divide the open field base into N_BOXES number of equal portions. 
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

# Function to define frames that animal is in various ROIs. 
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
frames_in_corner, time_in_corners = timestamps_in_roi(corner_rois)
        
#%% Locomotor measures
# Preallocate variables. 
distance_travelled = np.empty(len(paths_vid))
velocity = np.empty(len(paths_vid))

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
        
#%% Create data dictionary for summary data. 
data_dict = {
    'Animal_ID': paths_vid,
    'Sex': sex_key,
    'Treatment': treatment_key,
    'Avg_time': average_time,
    'Latency_explore': latency_exploration,
    'Object_one_time': object_one_time,
    'Object_two_time': object_two_time,
    'Proportion_time': object_one_time/object_two_time,
    'Time_Edges': time_in_edges,
    'Time_Corners': time_in_corners,
    'Distance': distance_travelled,
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
    


if export_to_csv:
    data_final.to_csv('NOE_dataFinal.csv')
    data_final_zscored.to_csv('NOE_dataFinal_zscored.csv')
    
    stats_results['Avg_time'].to_csv('NOE_stats_avgTime.csv') # Figure 3b

    ## Stats for extended figures
    stats_results['Latency_explore'].to_csv('NOE_stats_latencyExplore.csv')
    stats_results['Time_Edges'].to_csv('NOE_stats_timeEdges.csv')
    stats_results['Time_Corners'].to_csv('NOE_stats_timeCorners.csv')
    stats_results['Distance'].to_csv('NOE_stats_distance.csv')


# %% Plot animal traces for visualization if you like.
if plot_traces:
    for video_idx, video_path in enumerate(paths_vid):
        psy_beh.plot_traces(video_path, body_coords[video_idx])

