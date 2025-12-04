# -*- coding: utf-8 -*-
"""
Created on Sat Jul  5 18:42:30 2025

@author: olu

SIT Analysis:
Instructions:
    1. Put all video and coordinate files from a single experiment in a folder.
    2. Adjust the variables in the "Variables to change" section to match your settings.
    3. Run the script. Follow command-line instructions as prompted.
    4. Update Supplementary table 1 (meta data sheet) with data in the "data_all" variable.

The variables in data_all include summary data for the following measures:
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
FOLDER_PATH = r"Y:/PsyBAANC/paperExperiments/SIT_SM/SIT_all"  # path to folder with data
VIDEO_TYPE = "avi"  # options: "mp4", "avi", others also likely OK.
COORDINATE_FILE_TYPE = "csv"  # options: "csv", "xlsx"
START_SEC = 0  # the time in seconds that you wish to begin analysis.
END_SEC = 10*60  # the time in seconds that you wish to end analysis.

LENGTH_CM = 40  # true size of your open field box in cm

X_NOSE_COORDINATE_INDEX = 1  # Index of nose x-coord column in coordinates file (index starts at 0).
Y_NOSE_COORDINATE_INDEX = 2  # Index of nose y-coord column in coordinates file (index starts at 0).
X_BODY_COORDINATE_INDEX = 4  # Index of body x-coord column in coordinates file (index starts at 0).
Y_BODY_COORDINATE_INDEX = 5  # Index of body y-coord column in coordinates file (index starts at 0).
ROW_INDEX = 4  # Row number that position data starts in coordinate files
DATA_DLC = True  # Is your data from deeplabcut (DLC)? true or false.
COORDINATES_CM = False  # Are your coordinates in centimeters? (And not pixels)

OBJECT_ONE_SHAPE = "ellipse"  # options: "circle", "ellipse", "rectangle", "polygon"
OBJECT_TWO_SHAPE = "ellipse"  # options: "circle", "ellipse", "rectangle", "polygon"

sex_key = ["M"]*14 + ["F"]*14  # List of sex of the animals, i.e., ["M", "F", "M", etc.]
treatment_key = (["P"]*3 + ["S"]*7 + ["P"]*4 +
                 ["P"]*3 + ["S"]*6 + ["P"]*5)  # List of treatment of animals.

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
empty_chamber = [None]*(len(paths_vid))
social_chamber = [None]*(len(paths_vid))

empty_cup = [None]*len(paths_vid)
social_cup = [None]*len(paths_vid)
empty_cup_roi = [None]*len(paths_vid)
social_cup_roi = [None]*len(paths_vid)

print("""
      DEFINE ROIs (Follow command-line instructions) or READ IN SAVED ROIs.
      """)
for video_idx, video_file in enumerate(paths_vid):
    empty_chamber[video_idx] = psy_beh.get_roi_coords(video_file, "empty_chamber_base",
                                                      body_coords[video_idx])
    social_chamber[video_idx] = psy_beh.get_roi_coords(video_file, "social_chamber_base",
                                                       body_coords[video_idx])
    empty_cup[video_idx] = psy_beh.get_roi_mask(video_file,
                                                "empty_cup", shape=OBJECT_ONE_SHAPE)
    social_cup[video_idx] = psy_beh.get_roi_mask(video_file,
                                                 "social_cup", shape=OBJECT_TWO_SHAPE)
    empty_cup_roi[video_idx] = psy_beh.resize_mask_by_cm(empty_cup[video_idx], cm_offset=3,
                                                         cm_per_pixel=cm_to_pixels[video_idx])
    social_cup_roi[video_idx] = psy_beh.resize_mask_by_cm(social_cup[video_idx], cm_offset=3,
                                                          cm_per_pixel=cm_to_pixels[video_idx])

# %% Object exploration analysis
# Exploration is defined as:
# - nose within 3 cm of cup, body is not within cup mask. (animal cannot be sitting on object)
# - animal must also be looking at object.

START_FRAME = np.round(START_SEC*FRAMERATE).astype(int)
END_FRAME = np.round(END_SEC*FRAMERATE).astype(int)


def get_exploration_metrics(cup_roi, cup):
    """Get exploration metrics for cup object.

    Parameters:
    - cup_roi: list of cup ROI masks (expanded by 3 cm)
    - cup: list of cup ROI masks
    Returns:
    - cup_time: time spent exploring cup (s)
    - cup_latency: latency to first visit to cup (s)
    - cup_visits: number of visits to cup
    """
    cup_time = np.empty(len(paths_vid))
    cup_latency = np.empty(len(paths_vid))
    cup_visits = np.empty(len(paths_vid))
    for vid_i in range(len(paths_vid)):
        empty_cup_nose_frames = psy_beh.get_timepoints_in_mask(nose_coords[vid_i],
                                                               cup_roi[vid_i])
        empty_cup_body_frames = psy_beh.get_timepoints_in_mask(body_coords[vid_i],
                                                               cup[vid_i])
        common_frames = np.intersect1d(empty_cup_nose_frames, empty_cup_body_frames)
        empty_cup_exploration = np.setdiff1d(empty_cup_nose_frames, common_frames)

        # Animal must also be looking towards the object.
        m = cv2.moments(cup[vid_i])
        object_coord = np.array([m["m10"] / m["m00"], m["m01"] / m["m00"]]).reshape(1, 2)
        time_looking_empty_cup = psy_beh.timepoints_looking_at_object(nose_coords[vid_i],
                                                                      body_coords[vid_i],
                                                                      object_coord,
                                                                      angle_thresh_deg=30)

        empty_cup_frames = np.intersect1d(empty_cup_exploration, time_looking_empty_cup)
        empty_cup_frames = empty_cup_frames[(empty_cup_frames >= START_FRAME[vid_i]) &
                                            (empty_cup_frames <= END_FRAME[vid_i])]

        cup_time[vid_i] = len(empty_cup_frames)/FRAMERATE[vid_i]

        empty_cup_frames_visits = psy_beh.split_into_visits(empty_cup_frames,
                                                            min_length=FRAMERATE[vid_i]/6)
        cup_visits[vid_i] = len(empty_cup_frames_visits)
        if cup_visits[vid_i] > 0:
            cup_latency[vid_i] = empty_cup_frames_visits[0][0]/FRAMERATE[vid_i]
        else:
            cup_latency[vid_i] = END_SEC-START_SEC

    return cup_time, cup_latency, cup_visits


empty_cup_time, empty_cup_latency, empty_cup_visits = get_exploration_metrics(empty_cup_roi,
                                                                              empty_cup)
social_cup_time, social_cup_latency, social_cup_visits = get_exploration_metrics(social_cup_roi,
                                                                                 social_cup)

total_time = empty_cup_time + social_cup_time
social_index = (social_cup_time-empty_cup_time)/total_time

# %% Time spent in chambers.
LENGTH_EXPERIMENT = END_SEC-START_SEC

frames_empty = [None]*len(paths_vid)
frames_social = [None]*len(paths_vid)

time_empty = np.empty(len(paths_vid))
time_social = np.empty(len(paths_vid))
time_center = np.empty(len(paths_vid))

distance_empty = np.empty(len(paths_vid))
distance_social = np.empty(len(paths_vid))

for video_idx, video_path in enumerate(paths_vid):
    frames_empty[video_idx] = psy_beh.find_timepoints_in_coords(body_coords[video_idx],
                                                                empty_chamber[video_idx])
    frames_empty[video_idx] = frames_empty[video_idx][
        (frames_empty[video_idx] >= START_FRAME[video_idx]) &
        (frames_empty[video_idx] <= END_FRAME[video_idx])]
    time_empty[video_idx] = len(frames_empty[video_idx])/FRAMERATE[video_idx]
    frames_social[video_idx] = psy_beh.find_timepoints_in_coords(body_coords[video_idx],
                                                                 social_chamber[video_idx])
    frames_social[video_idx] = frames_social[video_idx][
        (frames_social[video_idx] >= START_FRAME[video_idx]) &
        (frames_social[video_idx] <= END_FRAME[video_idx])]
    time_social[video_idx] = len(frames_social[video_idx])/FRAMERATE[video_idx]
    time_center[video_idx] = LENGTH_EXPERIMENT - time_empty[video_idx] - time_social[video_idx]

time_center = time_center/LENGTH_EXPERIMENT*100
time_empty = time_empty/LENGTH_EXPERIMENT*100
time_social = time_social/LENGTH_EXPERIMENT*100

# %% Locomotion metrics
# Preallocate variables.
distance_travelled = np.empty(len(paths_vid))
velocity = np.empty(len(paths_vid))

distance_empty = np.empty(len(paths_vid))
distance_social = np.empty(len(paths_vid))

for video_idx, vid_path in enumerate(paths_vid):
    # get distance array: distance travelled per frame.
    diff_array = np.diff(body_coords[video_idx], axis=0)
    diff_array[:, 0] = diff_array[:, 0] * cm_to_pixels[video_idx]
    diff_array[:, 1] = diff_array[:, 1] * cm_to_pixels[video_idx]
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

    distance_array_empty = distance_array[(frames_empty[video_idx])]
    distance_empty[video_idx] = sum(distance_array_empty)

    distance_array_social = distance_array[(frames_social[video_idx])]
    distance_social[video_idx] = sum(distance_array_social)

# %% Create data dictionary for summary data.
data_summary_dict = {
    'Animal_ID': paths_vid,
    'Sex': sex_key,
    'Treatment': treatment_key,
    'Institution': INSTITUTION,
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

data_summary_df = pd.DataFrame(data_summary_dict)
column_names = data_summary_df.columns[4:].tolist()
column_labels = ["Empty cup time (s)", "Social cup time (s)", "Social preference index",
                 "Latency empty cup (s)", "Latency social cup (s)",
                 "Visits empty cup", "Visits social cup",
                 "Empty chamber time (%)", "Social chamber time (%)", "Center chamber time (%)",
                 "Distance empty chamber (cm)", "Distance social chamber (s)",
                 "Distance travelled (cm)", "Velocity (cm/s)"]
data_all = data_summary_df

# %% Do statistical tests.
stats_summary_results = {}
for col in column_names:
    stats_summary_results[col] = psy_stats.stats_treatment_sex(data_summary_df, col)

# %% Plot data.
os.makedirs(os.path.join(FOLDER_PATH, "saved_data"), exist_ok=True)

fig, ax = plt.subplots(1, len(column_names), figsize=(1.2*len(column_names), 1.5))
for col_i, col in enumerate(column_names):
    psy_plot.plot_individual_lab(ax[col_i], data_summary_df, col, column_labels[col_i],
                                 INSTITUTION, stats_summary_results[col]["significance"])
plt.savefig(os.path.join(FOLDER_PATH, "saved_data", "summary_data_plots.png"))
plt.close()

# %% Plot animal traces for visualization.
for video_idx, video_path in enumerate(paths_vid):
    fig, ax = plt.subplots()
    psy_plot.plot_traces(ax, video_path, body_coords[video_idx], nose_coords[video_idx])
    # show object rois.
    ax.imshow(empty_cup_roi[video_idx], cmap='Reds', alpha=0.25)
    ax.imshow(social_cup_roi[video_idx], cmap='Reds', alpha=0.25)
    # show open field base rectangle
    psy_plot.plot_roi_coords(ax, empty_chamber[video_idx])
    psy_plot.plot_roi_coords(ax, social_chamber[video_idx])
    plt.savefig(os.path.join(FOLDER_PATH, "saved_data", f"{video_idx:02}_trace.png"))
    plt.close()

# %% Save raw data if wanted.
data_all.to_csv(os.path.join(FOLDER_PATH, "saved_data", 'SIT_data.csv'), index=False)
