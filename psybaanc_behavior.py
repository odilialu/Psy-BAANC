# -*- coding: utf-8 -*-
"""
Created on Mon Jun 30 10:10:56 2025

@author: olu
Functions used in psy-baanc behavioral video analyses.

"""
import os
import cv2
import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#%%
def get_file_paths(directory, file_type):
    """
    Parameters:
    - directory: string path to files
    - file_type : string of file type
    
    Returns:
    - file_paths : List of all the file paths in a directory with specified file type.
    """
    
    file_paths = []
    filenames = os.listdir(directory)
    for file in filenames:
        if file.endswith(file_type):
            joined = os.path.join(directory, file)
            file_paths.append(joined)
    return file_paths

#%%
def get_roi_from_frame(video_file_path, roi_name, coordinates = None):
    """
    Opens a video and allows the user to select a rectangular ROI on a frame. 
    If the ROI already exists, load the ROI.

    Parameters:
    - video_file_path: Path to the video file.

    Returns:
    - roi_coords: A tuple (x_min, x_max, y_min, y_max) representing the selected region.
    """

    # Folder path where rois are saved 
    video_folder = os.path.dirname(video_file_path)
    filename = os.path.split(video_file_path)[-1]
    filename_cut = os.path.splitext(filename)[0]
    path_roi_folder = video_folder + "/" + roi_name
    pathROI = path_roi_folder + "/" + roi_name + "_" + filename_cut + ".pickle"

    if os.path.exists(pathROI):
        with open(pathROI, "rb") as file:
            roi_coords = pickle.load(file)

    else:
        cap = cv2.VideoCapture(video_file_path)
        ret, frame = cap.read()  # Read the first frame
        cap.release()  # Close the video file

        if not ret:
            print("Error: Could not read a frame.")
            return None
        
        if coordinates is not None:
            coordinates= np.round(coordinates)
            for i in range(1, len(coordinates)):
                pt1 = tuple(coordinates[i - 1].astype(int))
                pt2 = tuple(coordinates[i].astype(int))
                cv2.line(frame, pt1, pt2, color=(0, 255, 0), thickness=1)  # green line

        # Select ROI manually
        roi = cv2.selectROI("Select ROI", frame, fromCenter=False, showCrosshair=True)
        cv2.destroyAllWindows()  # Close the window after selection

        # Convert ROI to (x_min, x_max, y_min, y_max)
        x_min, y_min, width, height = roi
        x_max, y_max = x_min + width, y_min + height
        roi_coords = (x_min, x_max, y_min, y_max) 

        # save the roi coordinates 
        if not os.path.exists(path_roi_folder):
            os.makedirs(path_roi_folder)

        with open(path_roi_folder + "/" + roi_name + "_" + filename_cut + ".pickle", "wb") as file:
            pickle.dump(roi_coords, file)

    print("ROI Coordinates: ", roi_coords)
    return roi_coords

#%%
def find_time_in_location(xy_coords, timestamps, roi):
    """
    Find timestamps when the mouse is in a specified location.

    Parameters:
    - xy_coords: np.array of shape (N,2) with x and y coordinates.
    - timestamps: np.array of shape (N,) with corresponding timestamps.
    - roi: tuple (x_min, x_max, y_min, y_max) defining the region of interest.

    Returns:
    - times_in_roi: List of timestamps when the mouse is inside the ROI.
    """

    x_min, x_max, y_min, y_max = roi

    # Find indices where the mouse is inside the ROI
    in_roi = (xy_coords[:, 0] >= x_min) & (xy_coords[:, 0] <= x_max) & \
             (xy_coords[:, 1] >= y_min) & (xy_coords[:, 1] <= y_max)

    # Get timestamps corresponding to these indices
    times_in_roi = timestamps[in_roi==True]

    return times_in_roi

#%%
def get_body_parts(data_dlc, body_part_idx_x = 13, likelihood=0.6):
    """
    Return a data frame that includes deeplabcut coordinates over time.
    Coordinates with likelihood below threshold cutoff are interpolated.
    
    Parameters:
    - data_dlc: data frame with deeplabcut data
    - body_part_idx_x: the index of the x-coordinates for the body part of interest
    - likelihood: likelihood cutoff threshold
    """

    data = {'x': data_dlc.iloc[:, body_part_idx_x],
            'y': data_dlc.iloc[:, body_part_idx_x+1],
            'likelihood': data_dlc.iloc[:, body_part_idx_x+2]}
    df = pd.DataFrame(data)

    df['x'] = np.where(df['likelihood']<=likelihood, np.nan, df['x'])
    df['y'] = np.where(df['likelihood']<=likelihood, np.nan, df['y'])

    df = df.interpolate().to_numpy()[:, 0:2]

    return df

#%%
def zscore_dataframe(data_final):
    m_sal = data_final[(data_final["Sex"] == "M") & (data_final["Treatment"] == "S")]
    f_sal = data_final[(data_final["Sex"] == "F") & (data_final["Treatment"] == "S")]
    m_psi = data_final[(data_final["Sex"] == "M") & (data_final["Treatment"] == "P")]
    f_psi = data_final[(data_final["Sex"] == "F") & (data_final["Treatment"] == "P")]
    
    m_sal_mean = m_sal.iloc[:, 3:].mean()
    f_sal_mean = f_sal.iloc[:, 3:].mean()
    m_sal_std = m_sal.iloc[:, 3:].std()
    f_sal_std = f_sal.iloc[:, 3:].std()
    
    
    def zscore(original_data, sal_mean, sal_std):
        zscored_data = original_data[["Animal_ID", "Sex", "Treatment"]].copy()
        for col in original_data.columns[3:]:
            zscore = np.empty(len(original_data))
            for animal in range(len(original_data)):
                zscore[animal] = (original_data[col].iloc[animal] - sal_mean[col])/sal_std[col]
            zscored_data[col] = zscore
        return zscored_data
    
    m_sal_zscore = zscore(m_sal, m_sal_mean, m_sal_std)
    f_sal_zscore = zscore(f_sal, f_sal_mean, f_sal_std)
    m_psi_zscore = zscore(m_psi, m_sal_mean, m_sal_std)
    f_psi_zscore = zscore(f_psi, f_sal_mean, f_sal_std)
    
    data_final_zscore = pd.concat([m_sal_zscore, f_sal_zscore, m_psi_zscore, f_psi_zscore], axis=0)

    return data_final_zscore

#%%
def plot_traces(video, coordinates):
    cap  = cv2.VideoCapture(video)
    ret, image = cap.read()
    
    plt.figure()
    plt.imshow(image)
    plt.plot(coordinates[:, 0], coordinates[:, 1], lw=1)
    plt.show()
    