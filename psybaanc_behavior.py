# -*- coding: utf-8 -*-
"""
Created on Mon Jun 30 10:10:56 2025

@author: olu
Functions used in psy-baanc behavioral video analyses.

"""
#%%
# Import packages
import os
import cv2
import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Global drawing variables
drawing = False
roi_type = None
ix, iy = -1, -1
points = []
temp_img = None
mask = None

#%%
def get_file_paths(directory, file_type):
    """
    From the directory, get a list of all the files of specified file type

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
    - timestamps: np.array of shape (N,) with corresponding timestamps. Note times start at 1, not 0.
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
    
    Returns:
    - df: data frame with x, y coordinates over time.
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
    """
    Z-scores the data in the data_final data frame, excluding the first three columns
    which are identity columns
    
    Parameters:
    - data_final : pd.DataFrame containing the data to be z-scored. The first three columns are identity columns and not data columns.
    
    Returns 
    - data_final_zscore : pd.DataFrame containing the z-scored data.
    """
    m_sal = data_final[(data_final["Sex"] == "M") & (data_final["Treatment"] == "S")]
    f_sal = data_final[(data_final["Sex"] == "F") & (data_final["Treatment"] == "S")]
    m_psi = data_final[(data_final["Sex"] == "M") & (data_final["Treatment"] == "P")]
    f_psi = data_final[(data_final["Sex"] == "F") & (data_final["Treatment"] == "P")]
    
    m_sal_mean = m_sal.iloc[:, 3:].mean()
    f_sal_mean = f_sal.iloc[:, 3:].mean()
    m_sal_std = m_sal.iloc[:, 3:].std()
    f_sal_std = f_sal.iloc[:, 3:].std()
    
    
    def zscore(original_data, sal_mean, sal_std):
        """
        Z-scores the data in a DataFrame, excluding the first three columns
        which are assumed to be identity columns.
        
        Parameters:
        original_data : pd.DataFrame containing the data to be z-scored. The first three columns are identity columns and not data columns.
        - sal_mean : pd.Series containing the mean values for the saline control treatment.
        - sal_std : pd.Series containing the standard deviation values for the saline control treatment.
        
        Returns:
        - zscored_data : pd.DataFrame containing the z-scored data.
        """
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
    """
    Plots the coordinates on the first frame of a video.
    Parameters:
    - video: Path to the video file.
    - coordinates: np.array of shape (N, 2) with x and y coordinates.
    """
    cap  = cv2.VideoCapture(video)
    ret, image = cap.read()
    
    plt.figure()
    plt.imshow(image)
    plt.plot(coordinates[:, 0], coordinates[:, 1], lw=1)
    plt.show()

#%%
def draw_callback(event, x, y, flags, param):
    """ Handles mouse events for drawing ROIs on the image.
    Parameters:
    - event: Mouse event type (e.g., cv2.EVENT_LBUTTONDOWN, cv2.EVENT_MOUSEMOVE, cv2.EVENT_LBUTTONUP).
    - x: x-coordinate of the mouse event.
    - y: y-coordinate of the mouse event.
    - flags: Additional flags (not used here).
    - param: Additional parameters (the image to draw on).
    
    This function allows the user to draw different types of ROIs (rectangle, circle, ellipse, polygon) on the image.
    It uses global variables to keep track of the drawing state and the temporary image being modified.
    """
    # Global drawing variables
    global ix, iy, drawing, temp_img, mask, points

    if roi_type in ['rectangle', 'circle', 'ellipse']:
        if event == cv2.EVENT_LBUTTONDOWN:
            drawing = True
            ix, iy = x, y
        elif event == cv2.EVENT_MOUSEMOVE and drawing:
            temp_img = param.copy()
            if roi_type == 'rectangle':
                cv2.rectangle(temp_img, (ix, iy), (x, y), (0, 255, 0), 2)
            elif roi_type == 'circle':
                radius = int(np.hypot(x - ix, y - iy))
                cv2.circle(temp_img, (ix, iy), radius, (0, 255, 0), 2)
            elif roi_type == 'ellipse':
                axes = (abs(x - ix), abs(y - iy))
                center = (ix, iy)
                cv2.ellipse(temp_img, center, axes, 0, 0, 360, (0, 255, 0), 2)
            cv2.imshow("Select ROI", temp_img)
        elif event == cv2.EVENT_LBUTTONUP:
            drawing = False
            if roi_type == 'rectangle':
                cv2.rectangle(mask, (ix, iy), (x, y), 1, -1)
            elif roi_type == 'circle':
                radius = int(np.hypot(x - ix, y - iy))
                cv2.circle(mask, (ix, iy), radius, 1, -1)
            elif roi_type == 'ellipse':
                axes = (abs(x - ix), abs(y - iy))
                center = (ix, iy)
                cv2.ellipse(mask, center, axes, 0, 0, 360, 1, -1)
            cv2.destroyWindow("Select ROI")

    elif roi_type == 'polygon':
        if event == cv2.EVENT_LBUTTONDOWN:
            points.append((x, y))
            if len(points) > 1:
                cv2.line(temp_img, points[-2], points[-1], (0, 255, 0), 2)
            cv2.circle(temp_img, (x, y), 3, (0, 255, 0), -1)
            cv2.imshow("Select ROI", temp_img)
        elif event == cv2.EVENT_RBUTTONDOWN:
            if len(points) >= 3:
                cv2.polylines(temp_img, [np.array(points)], isClosed=True, color=(0, 255, 0), thickness=2)
                cv2.fillPoly(mask, [np.array(points)], 1)
                cv2.destroyWindow("Select ROI")
            else:
                print("Need at least 3 points for polygon")
                
                
def get_roi_flexible(video_file_path, roi_name, shape='rectangle', coordinates=None):
    """
    Opens a video file and allows the user to draw a flexible ROI on the first frame.
    If the ROI already exists, it loads the ROI from a pickle file.

    Parameters:
    - video_file_path: Path to the video file.
    - roi_name: Name of the ROI to be saved.
    - shape: Shape of the ROI to be drawn ('rectangle', 'circle', 'ellipse', 'polygon').
    - coordinates: Optional coordinates to draw lines on the frame before selecting the ROI.

    Returns:
    - mask: A binary mask of the selected ROI.
    """
    global roi_type, temp_img, mask, points
    
    # Folder path where rois are saved 
    video_folder = os.path.dirname(video_file_path)
    filename = os.path.split(video_file_path)[-1]
    filename_cut = os.path.splitext(filename)[0]
    path_roi_folder = video_folder + "/" + roi_name
    pathROI = path_roi_folder + "/" + roi_name + "_" + filename_cut + ".pickle"

    if os.path.exists(pathROI):
        with open(pathROI, "rb") as file:
            mask = pickle.load(file)

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

        print(f"Draw {shape.upper()} ROI for {roi_name}. Left-click to define points of polygon, right click to close polygon and move onto next video.")
    
        roi_type = shape
        temp_img = frame.copy()
        mask = np.zeros(frame.shape[:2], dtype=np.uint8)
        points = []
    
        cv2.namedWindow("Select ROI")
        cv2.setMouseCallback("Select ROI", draw_callback, frame)
    
        cv2.imshow("Select ROI", frame)
        cv2.waitKey(0)
        cv2.destroyAllWindows()
        
        # save the mask coordinates 
        if not os.path.exists(path_roi_folder):
            os.makedirs(path_roi_folder)

        with open(path_roi_folder + "/" + roi_name + "_" + filename_cut + ".pickle", "wb") as file:
            pickle.dump(mask, file)

    return mask

#%%
def get_timepoints_in_mask(coordinates, mask):
    """    
    Given coordinates and a mask, returns the timepoints where the coordinates are inside the mask.
    
    Parameters:
    - coordinates: np.array of shape (N, 2) with x and y coordinates.
    - mask: binary mask of shape (height, width) where 1 indicates the area of interest.
    
    Returns:
    - timepoints_in_mask: np.array of indices where the coordinates are inside the mask.
    """
    x = np.round(coordinates[:, 0]).astype(int)
    y = np.round(coordinates[:, 1]).astype(int)

    # Ensure x and y are within bounds of the mask
    in_bounds = (x >= 0) & (x < mask.shape[1]) & (y >= 0) & (y < mask.shape[0])
    
    # Check if each in-bounds coordinate is inside the mask
    inside_mask = np.zeros_like(x, dtype=bool)
    inside_mask[in_bounds] = mask[y[in_bounds], x[in_bounds]].astype(bool)

    # Get the timepoints (indices) where coordinate is inside the mask
    timepoints_in_mask = np.where(inside_mask)[0]
    
    return timepoints_in_mask

#%%
def split_into_visits(frames, min_length = 1):
    """
    Splits a list of frames into visits based on gaps in the frame numbers.
    
    Parameters:
    - frames: List or array of frame numbers.
    - min_length: Minimum length in frames of a visit to be considered valid.

    Returns:
    - List of visits, where each visit is a list of frame numbers.
    """
    frames = np.sort(np.array(frames))
    if len(frames) == 0:
        return []

    # Find indices where the difference between consecutive frames is > 1
    split_indices = np.where(np.diff(frames) > 1)[0] + 1

    # Split at those points
    visits = np.split(frames, split_indices)

    # Filter visits by minimum length
    filtered_visits = [visit.tolist() for visit in visits if len(visit) >= min_length]

    return filtered_visits

#%%
def calibrate_pixels_to_cm(video_path, real_world_cm, frame_number=0):
    """
    Opens a frame from a video and lets user click two points to define a known distance.
    Returns pixel-to-cm conversion factor. 

    Parameters:
    - video_path: Path to the video file.
    - real_world_cm: Known distance in centimeters between the two clicked points.
    - frame_number: Frame number to read from the video (default is 0).

    Returns:
    - cm_per_pixel: Conversion factor from pixels to centimeters.
    
    """
    # Folder path where rois are saved 
    video_folder = os.path.dirname(video_path)
    filename = os.path.split(video_path)[-1]
    filename_cut = os.path.splitext(filename)[0]
    path_calibration_folder = video_folder + "/" + "calibration"
    path_calibration = path_calibration_folder + "/" + "calibration" + "_" + filename_cut + ".pickle"

    if os.path.exists(path_calibration):
        with open(path_calibration, "rb") as file:
            cm_per_pixel = pickle.load(file)

    else:
        
        # Load video
        cap = cv2.VideoCapture(video_path)
        cap.set(cv2.CAP_PROP_POS_FRAMES, frame_number)
        ret, frame = cap.read()
        cap.release()
        
        if not ret:
            print("Failed to read frame.")
            return None
    
        points = []
    
        def click_event(event, x, y, flags, param):
            if event == cv2.EVENT_LBUTTONDOWN and len(points) < 2:
                points.append((x, y))
                cv2.circle(display_frame, (x, y), 5, (0, 255, 0), -1)
                cv2.imshow("Click two points", display_frame)
    
        display_frame = frame.copy()
        cv2.imshow("Click two points", display_frame)
        cv2.setMouseCallback("Click two points", click_event)
    
        print("Click two points to define a known real-world distance.")
        while len(points) < 2:
            if cv2.waitKey(1) & 0xFF == 27:  # Escape key to quit
                print("Canceled.")
                cv2.destroyAllWindows()
                return None
    
        cv2.destroyAllWindows()
    
        # Calculate pixel distance
        pt1, pt2 = points
        pixel_dist = np.linalg.norm(np.array(pt1) - np.array(pt2))
    
        # Ask user for real-world distance
        cm_per_pixel = real_world_cm / pixel_dist
        
        # save the mask coordinates 
        if not os.path.exists(path_calibration_folder):
            os.makedirs(path_calibration_folder)

        with open(path_calibration_folder + "/" + "calibration" + "_" + filename_cut + ".pickle", "wb") as file:
            pickle.dump(cm_per_pixel, file)

    return cm_per_pixel