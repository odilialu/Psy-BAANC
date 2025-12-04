# -*- coding: utf-8 -*-
"""
Created on Sat Jul  5 14:17:43 2025 / Edited by cliu to export csv files

@author: olu

NOE Analysis:
Instructions:
    1. Put all video and coordinate files from a single experiment in a folder.
    2. Adjust the variables in the "Variables to change" section to match your settings.
    3. Run the script. Follow command-line instructions as prompted.
    4. Update Supplementary table 1 (meta data sheet) with data in the "data_all" variable.

The variables in data_all include summary data for the following measures:
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
FOLDER_PATH = r"Y:/PsyBAANC/paperExperiments/NOE_NOR/NOE/Analyze"  # path to folder with data
VIDEO_TYPE = "mp4"  # options: "mp4", "avi", others also likely OK.
COORDINATE_FILE_TYPE = "csv"  # options: "csv", "xlsx"
START_SEC = 0  # the time in seconds that you wish to begin analysis.
END_SEC = 10*60  # the time in seconds that you wish to end analysis.

LENGTH_CM = 50  # true size of your open field box in cm
N_BOXES = 25  # How many squares do you want to divide OF into? Keep at 25.

X_NOSE_COORDINATE_INDEX = 1  # Index of nose x-coord column in coordinates file (index starts at 0)
Y_NOSE_COORDINATE_INDEX = 2  # Index of nose y-coord column in coordinates file (index starts at 0)
X_BODY_COORDINATE_INDEX = 16  # Index of body x-coord column in coordinates file (index starts at 0)
Y_BODY_COORDINATE_INDEX = 17  # Index of body y-coord column in coordinates file (index starts at 0)
ROW_INDEX = 4  # Row number that position data starts in coordinate files
DATA_DLC = True  # Is your data from deeplabcut (DLC)? true or false.
COORDINATES_CM = False  # Are your coordinates in centimeters? (And not pixels)

OBJECT_ONE_SHAPE = "circle"  # options: "circle", "ellipse", "rectangle", "polygon"
OBJECT_TWO_SHAPE = "circle"  # options: "circle", "ellipse", "rectangle", "polygon"

sex_key = ["M"]*20 + ["F"]*20  # List of sex of the animals, i.e., ["M", "F", "M", etc.]
treatment_key = (["P"]*5 + ["S"]*5 + ["S"]*5 + ["P"]*5 +
                 ["S"]*5 + ["S"]*5 + ["P"]*5 + ["P"]*5)  # List of treatment of animals.

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
                                       X_BODY_COORDINATE_INDEX, Y_BODY_COORDINATE_INDEX, ROW_INDEX,
                                       DATA_DLC=DATA_DLC)
nose_coords = psy_beh.read_coordinates(paths_coordinates,
                                       X_NOSE_COORDINATE_INDEX, Y_NOSE_COORDINATE_INDEX, ROW_INDEX,
                                       DATA_DLC=DATA_DLC)

if COORDINATES_CM:
    cm_to_pixels = np.empty(len(paths_coordinates))
    for file_idx, coordinates in enumerate(body_coords):
        cm_to_pixels[file_idx] = psy_beh.calibrate_pixels_to_cm(paths_vid[file_idx],
                                                                real_world_cm=LENGTH_CM,
                                                                frame_number=0)
        body_coords[file_idx] = body_coords[file_idx]/cm_to_pixels[file_idx]
        nose_coords[file_idx] = nose_coords[file_idx]/cm_to_pixels[file_idx]

# %% Calibrate pixels to cm
cm_to_pixels = np.empty(len(paths_vid))
for video_idx, video_path in enumerate(paths_vid):
    cm_to_pixels[video_idx] = psy_beh.calibrate_pixels_to_cm(video_path,
                                                             real_world_cm=LENGTH_CM,
                                                             frame_number=0)

# %% Define all relevant ROIs if you have not already.
# Note: ROIs are saved for each video after defining them.
# If you want to re-define an ROI, you must delete the existing one.

# Preallocate lists.
open_field_base = [None]*(len(paths_vid))
edge_rois = [None]*(len(paths_vid))
center_rois = [None]*(len(paths_vid))
corner_rois = [None]*(len(paths_vid))

object_one = [None]*len(paths_vid)
object_two = [None]*len(paths_vid)
object_one_roi = [None]*len(paths_vid)
object_two_roi = [None]*len(paths_vid)


print("""
      DEFINE ROIs (Follow command-line instructions) or READ IN SAVED ROIs.
      """)
for video_idx, video_file in enumerate(paths_vid):
    open_field_base[video_idx] = psy_beh.get_roi_coords(video_file,
                                                        "open_field_base",
                                                        coordinates=body_coords[video_idx])

    object_one[video_idx] = psy_beh.get_roi_mask(video_file, "object_one", shape=OBJECT_ONE_SHAPE)
    object_two[video_idx] = psy_beh.get_roi_mask(video_file, "object_two", shape=OBJECT_TWO_SHAPE)
    object_one_roi[video_idx] = psy_beh.resize_mask_by_cm(object_one[video_idx],
                                                          cm_offset=3,
                                                          cm_per_pixel=cm_to_pixels[video_idx])
    object_two_roi[video_idx] = psy_beh.resize_mask_by_cm(object_two[video_idx],
                                                          cm_offset=3,
                                                          cm_per_pixel=cm_to_pixels[video_idx])

# %% Object exploration analysis
# Exploration is defined as:
# - nose within 3 cm of object, body is not in object mask (animal cannot be sitting on object)
# - animal must also be looking at object.

START_FRAME = np.round(START_SEC*FRAMERATE).astype(int)
END_FRAME = np.round(END_SEC*FRAMERATE).astype(int)


def get_exploration_metrics(object_roi, object_mask):
    """
    Function to get exploration metrics for a given object.
    Parameters:
    - object_roi: list of object roi masks for each video.
    - object_mask: list of object masks for each video.

    Returns:
    - object_time: array of time spent exploring object for each video.
    - object_latency: array of latency to first exploration for each video.
    - object_visits: array of number of visits to object for each video.
    """
    object_time = np.empty(len(paths_vid))
    object_latency = np.empty(len(paths_vid))
    object_visits = np.empty(len(paths_vid))
    for vid_i in range(len(paths_vid)):
        object_one_nose_frames = psy_beh.get_timepoints_in_mask(nose_coords[vid_i],
                                                                object_roi[vid_i])
        object_one_body_frames = psy_beh.get_timepoints_in_mask(body_coords[vid_i],
                                                                object_mask[vid_i])
        common_frames = np.intersect1d(object_one_nose_frames, object_one_body_frames)
        object_one_exploration = np.setdiff1d(object_one_nose_frames, common_frames)

        # Animal must also be looking towards the object.
        m = cv2.moments(object_mask[vid_i])
        object_coord = np.array([m["m10"] / m["m00"], m["m01"] / m["m00"]]).reshape(1, 2)
        time_looking_object_one = psy_beh.timepoints_looking_at_object(nose_coords[vid_i],
                                                                       body_coords[vid_i],
                                                                       object_coord,
                                                                       angle_thresh_deg=30)

        object_one_frames = np.intersect1d(object_one_exploration, time_looking_object_one)
        object_one_frames = object_one_frames[(object_one_frames >= START_FRAME[vid_i]) &
                                              (object_one_frames <= END_FRAME[vid_i])]

        object_time[vid_i] = len(object_one_frames)/FRAMERATE[vid_i]

        object_one_frames_visits = psy_beh.split_into_visits(object_one_frames,
                                                             min_length=int(FRAMERATE[vid_i]/6))
        object_visits[vid_i] = len(object_one_frames_visits)
        if object_visits[vid_i] > 0:
            object_latency[vid_i] = object_one_frames_visits[0][0]/FRAMERATE[vid_i]
        else:
            object_latency[vid_i] = END_SEC-START_SEC

    return object_time, object_latency, object_visits


object_one_time, object_one_latency, object_one_visits = get_exploration_metrics(object_one_roi,
                                                                                 object_one)
object_two_time, object_two_latency, object_two_visits = get_exploration_metrics(object_two_roi,
                                                                                 object_two)

total_time = object_one_time + object_two_time
average_time = total_time/2
total_visits = object_one_visits + object_two_visits
latency_exploration = np.min(np.vstack((object_one_latency, object_two_latency)).T, axis=1)

# %% Anxiety-related behaviors analysis
# divide the open field base into N_BOXES number of equal portions.
# From the boxes, define edge, center, and corner ROIs.
for video_idx, video_file in enumerate(paths_vid):
    open_field_width = open_field_base[video_idx][1] - open_field_base[video_idx][0]
    open_field_length = open_field_base[video_idx][3] - open_field_base[video_idx][2]
    cm_to_pixels_x = LENGTH_CM/open_field_width
    cm_to_pixels_y = LENGTH_CM/open_field_length

    n_rows = int(np.sqrt(N_BOXES))
    sub_width = open_field_width / n_rows
    sub_length = open_field_length / n_rows

    sub_rois = []
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


def timestamps_in_roi(roi_type, percent_time=True):
    """
    Function to get timestamps and time spent in a particular roi type.

    Parameters:
    - roi_type: list of lists of rois for each video.
    - percent_time: boolean, whether to return time as percent of total analysis time.

    Returns:
    - frames_in_roi: list of arrays of frames in roi for each video.
    - time_in_roi: array of time spent in roi for each video.
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
frames_in_corner, time_in_corners = timestamps_in_roi(corner_rois)

# %% Locomotor measures
# Preallocate variables.
distance_travelled = np.empty(len(paths_vid))
velocity = np.empty(len(paths_vid))

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

# %% Create data dictionary for summary data.
data_summary_dict = {
    'Animal_ID': paths_vid,
    'Sex': sex_key,
    'Treatment': treatment_key,
    'Institution': INSTITUTION,
    'Avg_time': average_time,
    'Latency_explore': latency_exploration,
    'Object_one_time': object_one_time,
    'Object_two_time': object_two_time,
    'Proportion_time': object_one_time/(object_two_time + object_one_time)*100,
    'Time_Edges': time_in_edges,
    'Time_Corners': time_in_corners,
    'Distance': distance_travelled,
    'Velocity': velocity
    }

data_summary_df = pd.DataFrame(data_summary_dict)
column_names = data_summary_df.columns[4:].tolist()
column_labels = ["Average exploration time (s)", "Latency first exploration (s)",
                 "Object 1 time (s)", "Object 2 time (s)", "Percent time (%)",
                 "Time in edges (%)", "Time in corners (%)",
                 "Distance travelled (cm)", "Velocity (cm/s)"]
data_all = data_summary_df

# %% Do statistical tests
stats_summary_results = {}
for col in column_names:
    stats_summary_results[col] = psy_stats.stats_treatment_sex(data_summary_df, col)

# %% Plot data
os.makedirs(os.path.join(FOLDER_PATH, "saved_data"), exist_ok=True)

fig, ax = plt.subplots(1, len(column_names), figsize=(1.2*len(column_names), 1.5))
for col_i, col in enumerate(column_names):
    psy_plot.plot_individual_lab(ax[col_i], data_summary_df, col, column_labels[col_i],
                                 INSTITUTION, stats_summary_results[col]["significance"])
plt.savefig(os.path.join(FOLDER_PATH, "saved_data", "summary_data_plots.png"))
plt.close()

# %% Plot animal traces for visualization
for video_idx, video_path in enumerate(paths_vid):
    fig, ax = plt.subplots()
    psy_plot.plot_traces(ax, video_path, body_coords[video_idx], nose_coords[video_idx])
    # show object rois.
    ax.imshow(object_one_roi[video_idx], cmap='Reds', alpha=0.25)
    ax.imshow(object_two_roi[video_idx], cmap='Reds', alpha=0.25)
    # show open field base rectangle
    psy_plot.plot_roi_coords(ax, open_field_base[video_idx])
    plt.savefig(os.path.join(FOLDER_PATH, "saved_data", f"{video_idx:02}_trace.png"))
    plt.close()

# %% Save raw data
data_all.to_csv(os.path.join(FOLDER_PATH, "saved_data", 'NOE_data.csv'), index=False)
