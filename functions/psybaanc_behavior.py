# -*- coding: utf-8 -*-
"""
Created on Mon Jun 30 10:10:56 2025

@author: olu
Functions used in psy-baanc behavioral video analyses.

"""
# %% Import packages
import os
import cv2
import pickle
import pandas as pd
import numpy as np
from shapely.geometry import Polygon


# %%
def get_file_paths(directory, file_type):
    """
    Get sorted list of file paths in a directory with a specific file type.

    Parameters:
    - directory: Path to the directory to search.
    - file_type: File extension to filter by (e.g., '.csv', '.xlsx').

    Returns:
    - file_paths: Sorted list of full file paths matching the file type.
    """

    file_paths = []
    filenames = os.listdir(directory)
    for file in filenames:
        if file.endswith(file_type) and (file.startswith('.') is False
                                         and file.startswith('~') is False):
            joined = os.path.join(directory, file)
            file_paths.append(joined)
    file_paths.sort()
    return file_paths


# %% Coordinate related functions
def dlc_body_parts(DATA_DLC, body_part_idx_x, likelihood=0.6, nan_begin=120):
    """
    Extracts x,y coordinates for a specific body part from DeepLabCut data,
    applying a likelihood threshold to filter out low-confidence points.

    Parameters:
    - DATA_DLC: DataFrame containing DeepLabCut output data.
    - body_part_idx_x: Index of the x-coordinate column for the body part.
    - likelihood: Likelihood threshold below which coordinates are set to NaN.

    Returns:
    - df: DataFrame with filtered x and y coordinates for the specified body part.
    """

    data = {'x': DATA_DLC.iloc[:, body_part_idx_x],
            'y': DATA_DLC.iloc[:, body_part_idx_x+1],
            'likelihood': DATA_DLC.iloc[:, body_part_idx_x+2]}
    df = pd.DataFrame(data)

    df['x'] = np.where(df['likelihood'] <= likelihood, np.nan, df['x'])
    df['y'] = np.where(df['likelihood'] <= likelihood, np.nan, df['y'])

    df = df.drop(columns=['likelihood'])

    if nan_begin is not None:
        diff_array = df.diff()
        distance_array = np.sqrt((diff_array["x"]**2) + (diff_array["y"]**2))
        jump_thresh = 5*np.std(distance_array, ddof=1) + np.mean(distance_array)
        indices_jump = np.where(distance_array > jump_thresh)[0]
        if len(indices_jump) > 0 and indices_jump[0] < nan_begin:
            index_max = indices_jump[indices_jump < nan_begin][-1] + 1
        else:
            index_max = 0
        df.iloc[:index_max, :] = np.nan

    return df


def read_coordinates(paths_coordinates, x_coordinate_i, y_coordinate_i, row_index,
                     DATA_DLC=False, likelihood=0.6):
    """
    Reads coordinate files and returns a list of numpy arrays with x,y coordinates.

    Parameters:
    - paths_coordinates: list of file paths to coordinate files
    - x_coordinate_i: index of x-coordinate column
    - y_coordinate_i: index of y-coordinate column
    - row_index: row index to start reading data from (to skip headers)
    - DATA_DLC: boolean indicating if data is from DeepLabCut
    - likelihood: likelihood threshold for DeepLabCut data

    Returns:
    - body_coords: list of numpy arrays with x,y coordinates for each file
    """

    coordinate_file_type = os.path.splitext(paths_coordinates[0])[-1]
    body_coords = [None]*len(paths_coordinates)

    for file_idx, coordinate_file in enumerate(paths_coordinates):
        if DATA_DLC is True:
            coordinates = pd.read_csv(coordinate_file, skiprows=list(range(0, row_index-2)))
            body_coords[file_idx] = dlc_body_parts(coordinates,
                                                   body_part_idx_x=x_coordinate_i,
                                                   likelihood=likelihood)
        else:
            if "xlsx" in coordinate_file_type:
                body_coords[file_idx] = pd.read_excel(coordinate_file,
                                                      usecols=[x_coordinate_i, y_coordinate_i],
                                                      skiprows=list(range(0, (row_index-2))))
            elif "csv" in coordinate_file_type:
                body_coords[file_idx] = pd.read_csv(coordinate_file,
                                                    usecols=[x_coordinate_i, y_coordinate_i],
                                                    skiprows=list(range(0, (row_index-2))))

            # Replace various non-numeric values with NaN
            body_coords[file_idx] = body_coords[file_idx].replace(
                ['-', '', ' ', 'NaN', 'nan', 'NULL', 'null'], np.nan)
            body_coords[file_idx] = body_coords[file_idx].astype(float)

        # Print to check data shapes.
        print(
            f"File {file_idx}: "
            f"Length of coordinates: {body_coords[file_idx].shape[0]}"
            )
        print(
            f"File {file_idx}: "
            f"NaN counts in coordinates colX, colY: {body_coords[file_idx].isnull().sum().tolist()}"
            )

        # Interpolate nans.
        body_coords[file_idx] = body_coords[file_idx].interpolate()
        body_coords[file_idx] = body_coords[file_idx].bfill().to_numpy()

    return body_coords


# %% Define ROIs
def get_roi_coords(video_file_path, roi_name, coordinates=None):
    """
    Opens a video file and allows the user to draw a rectangular ROI on the first frame.
    If the ROI already exists, it loads the ROI from a pickle file.

    Parameters:
    - video_file_path: Path to the video file.
    - roi_name: Name of the ROI to be saved.
    - coordinates: Optional coordinates to draw lines on the frame before selecting the ROI.

    Returns:
    - roi_coords: Tuple (x_min, x_max, y_min, y_max) defining the selected ROI.
    """

    # Folder path where rois are saved
    video_folder = os.path.dirname(video_file_path)
    filename = os.path.split(video_file_path)[-1]
    filename_cut = os.path.splitext(filename)[0]
    path_roi_folder = video_folder + "/" + roi_name
    path_roi = path_roi_folder + "/" + roi_name + "_" + filename_cut + ".pickle"

    if os.path.exists(path_roi):
        with open(path_roi, "rb") as file:
            roi_coords = pickle.load(file)
            print("Reading saved " + roi_name + " for " + video_file_path)

    else:
        cap = cv2.VideoCapture(video_file_path)
        cap.set(cv2.CAP_PROP_POS_FRAMES, 100)
        ret, frame = cap.read()
        cap.release()  # Close the video file

        if not ret:
            print("Error: Could not read a frame.")
            return None

        if coordinates is not None:
            coordinates = np.round(coordinates)
            for i in range(1, len(coordinates)):
                pt1 = tuple(coordinates[i - 1].astype(int))
                pt2 = tuple(coordinates[i].astype(int))
                cv2.line(frame, pt1, pt2, color=(0, 255, 0), thickness=1)  # green line

        # Select ROI manually
        print(f"Draw {roi_name} for {video_file_path}. Press ENTER when done.\n")
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


class ROISelector:
    def __init__(self, frame, roi_type="rectangle"):
        self.frame = frame
        self.roi_type = roi_type

        self.temp_img = frame.copy()
        self.mask = np.zeros(frame.shape[:2], dtype=np.uint8)     # 0/1 mask

        self.ix = None
        self.iy = None
        self.drawing = False
        self.points = []

    def clamp(self, x, y):
        h, w = self.frame.shape[:2]
        return max(0, min(x, w - 1)), max(0, min(y, h - 1))

    # ---------------- Mouse Callback ---------------- #
    def callback(self, event, x, y, flags, param):
        x, y = self.clamp(x, y)

        # ---- SHAPE ROIs ----
        if self.roi_type in ("rectangle", "circle", "ellipse"):

            if event == cv2.EVENT_LBUTTONDOWN:
                self.drawing = True
                self.ix, self.iy = x, y

            elif event == cv2.EVENT_MOUSEMOVE and self.drawing:
                self.temp_img = self.frame.copy()

                if self.roi_type == "rectangle":
                    cv2.rectangle(self.temp_img, (self.ix, self.iy), (x, y), (0, 255, 0), 2)

                elif self.roi_type == "circle":
                    r = int(np.hypot(x - self.ix, y - self.iy))
                    cv2.circle(self.temp_img, (self.ix, self.iy), r, (0, 255, 0), 2)

                elif self.roi_type == "ellipse":
                    axes = (abs(x - self.ix), abs(y - self.iy))
                    cv2.ellipse(self.temp_img, (self.ix, self.iy), axes, 0, 0, 360, (0, 255, 0), 2)

            elif event == cv2.EVENT_LBUTTONUP and self.drawing:
                self.drawing = False

                # Draw on mask using value 1 instead of 255
                if self.roi_type == "rectangle":
                    cv2.rectangle(self.mask, (self.ix, self.iy), (x, y), 1, -1)

                elif self.roi_type == "circle":
                    r = int(np.hypot(x - self.ix, y - self.iy))
                    cv2.circle(self.mask, (self.ix, self.iy), r, 1, -1)

                elif self.roi_type == "ellipse":
                    axes = (abs(x - self.ix), abs(y - self.iy))
                    cv2.ellipse(self.mask, (self.ix, self.iy), axes, 0, 0, 360, 1, -1)

        # ---- POLYGON ROI ----
        elif self.roi_type == "polygon":
            if event == cv2.EVENT_LBUTTONDOWN:
                self.points.append((x, y))

                self.temp_img = self.frame.copy()
                # draw lines between polygon points
                for i in range(1, len(self.points)):
                    cv2.line(self.temp_img, self.points[i - 1], self.points[i],
                             (0, 255, 0), 2)
                cv2.circle(self.temp_img, (x, y), 3, (0, 255, 0), -1)

            elif event == cv2.EVENT_RBUTTONDOWN:
                # close polygon
                if len(self.points) >= 3:
                    pts = np.array(self.points, dtype=np.int32)
                    cv2.fillPoly(self.mask, [pts], 1)
                self.points = []

    # ---------------- Main Loop ---------------- #
    def select(self):
        cv2.namedWindow("ROI Selector", cv2.WINDOW_NORMAL)
        cv2.setMouseCallback("ROI Selector", self.callback)
        cv2.imshow("ROI Selector", self.temp_img)
        cv2.waitKey(200)

        while True:
            k = cv2.waitKey(20) & 0xFF

            # ENTER closes polygon automatically
            if k == 13:
                if self.roi_type == "polygon" and len(self.points) >= 3:
                    pts = np.array(self.points, dtype=np.int32)
                    cv2.fillPoly(self.mask, [pts], 1)
                break

            if k in (27, ord('q')):  # ESC or q to cancel
                break

            cv2.imshow("ROI Selector", self.temp_img)

        cv2.destroyAllWindows()
        return self.mask


def get_roi_mask(video_file_path, roi_name, shape="rectangle", coordinates=None):
    video_folder = os.path.dirname(video_file_path)
    filename = os.path.basename(video_file_path)
    base = os.path.splitext(filename)[0]

    save_dir = os.path.join(video_folder, roi_name)
    save_path = os.path.join(save_dir, f"{roi_name}_{base}.pickle")

    # Load if exists
    if os.path.exists(save_path):
        with open(save_path, "rb") as f:
            return pickle.load(f)

    # Read first frame
    cap = cv2.VideoCapture(video_file_path)
    cap.set(cv2.CAP_PROP_POS_FRAMES, 100)
    ret, frame = cap.read()
    cap.release()

    if not ret:
        print("Error: Cannot read video frame.")
        return None

    # Optional guided lines
    if coordinates is not None:
        coords = np.round(coordinates).astype(int)
        for i in range(1, len(coords)):
            cv2.line(frame, tuple(coords[i - 1]), tuple(coords[i]), (0, 255, 0), 1)

    print(f"\n Defining ROI for {video_file_path}")
    if "polygon" in shape:
        print(f"Draw a tight {shape.upper()} ROI boundary for {roi_name.upper()}."
              "Left-click to define points of polygon,"
              "right click to close polygon. Press ENTER when done.")
    else:
        print(f"Draw a tight {shape.upper()} ROI boundary for {roi_name.upper()}."
              "Press and hold to draw shape, release to confirm. Press ENTER when done.")

    selector = ROISelector(frame, roi_type=shape)
    mask = selector.select()

    os.makedirs(save_dir, exist_ok=True)
    with open(save_path, "wb") as f:
        pickle.dump(mask, f)

    return mask


def resize_mask_by_cm(MASK, cm_offset, cm_per_pixel):
    """
    Resize a binary mask by a specified offset in centimeters.

    Parameters:
    - mask: binary mask of shape (height, width)
    - cm_offset: offset in centimeters to resize the mask (positive to expand, negative to shrink)
    - cm_per_pixel: conversion factor from pixels to centimeters

    Returns:
    - resized_mask: binary mask resized by the specified offset in centimeters.
    """
    pixel_offset = cm_offset / cm_per_pixel
    contours, _ = cv2.findContours(MASK.astype(np.uint8),
                                   cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    resized_mask = np.zeros_like(MASK, dtype=np.uint8)

    for cnt in contours:
        if cnt.shape[0] < 3:
            continue  # skip degenerate contours

        poly = Polygon(cnt.squeeze())
        if not poly.is_valid or poly.is_empty:
            continue

        poly_resized = poly.buffer(pixel_offset)
        if poly_resized.is_empty:
            continue

        # handle both Polygon and MultiPolygon
        polys = [poly_resized] if poly_resized.geom_type == 'Polygon' else list(poly_resized)

        for p in polys:
            if not p.is_valid or p.is_empty:
                continue
            coords = np.array(p.exterior.coords).round().astype(np.int32)
            cv2.fillPoly(resized_mask, [coords], 1)

    return resized_mask.astype(bool)


# %% ROI analyses
def find_timepoints_in_coords(xy_coords, roi):
    """
    Given x,y coordinates and a rectangular ROI, returns the timepoints where the
    coordinates are inside the ROI.

    Parameters:
    - xy_coords: np.array of shape (N, 2) with x and y coordinates.
    - roi: tuple (x_min, x_max, y_min, y_max) defining the rectangular ROI.

    Returns:
    - times_in_roi: np.array of indices where the coordinates are inside the ROI.
    """

    x_min, x_max, y_min, y_max = roi

    # Find indices where the mouse is inside the ROI
    in_roi = (xy_coords[:, 0] >= x_min) & (xy_coords[:, 0] <= x_max) & \
             (xy_coords[:, 1] >= y_min) & (xy_coords[:, 1] <= y_max)

    # Get timestamps corresponding to these indices
    times_in_roi = np.where(in_roi)[0]

    return times_in_roi


def get_timepoints_in_mask(coordinates, MASK):
    """
    Given x,y coordinates and a binary mask, returns the timepoints where the
    coordinates are inside the mask.

    Parameters:
    - coordinates: np.array of shape (N, 2) with x and y coordinates
    - mask: binary mask of shape (height, width)

    Returns:
    - timepoints_in_mask: np.array of indices where the coordinates are inside the mask.
    """
    x = np.round(coordinates[:, 0]).astype(int)
    y = np.round(coordinates[:, 1]).astype(int)

    # Ensure x and y are within bounds of the MASK
    in_bounds = (x >= 0) & (x < MASK.shape[1]) & (y >= 0) & (y < MASK.shape[0])

    # Check if each in-bounds coordinate is inside the MASK
    inside_mask = np.zeros_like(x, dtype=bool)
    inside_mask[in_bounds] = MASK[y[in_bounds], x[in_bounds]].astype(bool)

    # Get the timepoints (indices) where coordinate is inside the MASK
    timepoints_in_mask = np.where(inside_mask)[0]

    return timepoints_in_mask


# %% Pixel CM conversion-related functions
def calibrate_pixels_to_cm(video_path, real_world_cm, frame_number=100):
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
    path_calibration = f"{path_calibration_folder}/calibration_{filename_cut}.pickle"

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

        cv2.namedWindow("Click two points", cv2.WINDOW_NORMAL)
        cv2.resizeWindow("Click two points", frame.shape[1], frame.shape[0])
        cv2.setWindowProperty("Click two points", cv2.WND_PROP_AUTOSIZE, 1)

        display_frame = frame.copy()
        cv2.namedWindow("Click two points", cv2.WINDOW_AUTOSIZE)
        cv2.imshow("Click two points", display_frame)
        cv2.setMouseCallback("Click two points", click_event)

        print(f"Calibration for {video_path}: "
              "Click two points to define a known real-world distance.")
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

        # save the MASK coordinates
        if not os.path.exists(path_calibration_folder):
            os.makedirs(path_calibration_folder)

        with open(f"{path_calibration_folder}/calibration_{filename_cut}.pickle", "wb") as file:
            pickle.dump(cm_per_pixel, file)

    return cm_per_pixel


# %% Others
def split_into_visits(frames, min_length=1):
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


def timepoints_looking_at_object(nose_coordinates, body_coordinates, object_coordinates,
                                 angle_thresh_deg=30):
    """
    Determines timepoints when the subject is looking at an object based on head direction.

    Parameters:
    - nose_coordinates: np.array of shape (N, 2) with x,y coordinates of the nose.
    - body_coordinates: np.array of shape (N, 2) with x,y coordinates of the body.
    - object_coordinates: np.array of shape (N, 2) with x,y coordinates of the object.
    - angle_thresh_deg: Angle threshold in degrees to consider "looking at" the object.

    Returns:
    - looking_timepoints: np.array of indices where the subject is looking at the object.
    """
    x_body = body_coordinates[:, 0]
    y_body = body_coordinates[:, 1]

    x_nose = nose_coordinates[:, 0]
    y_nose = nose_coordinates[:, 1]

    x_obj = object_coordinates[:, 0]
    y_obj = object_coordinates[:, 1]

    # Vectors: head direction and object direction
    head_vec = np.stack([x_nose - x_body, y_nose - y_body], axis=1)
    obj_vec = np.stack([x_obj - x_body, y_obj - y_body], axis=1)

    # Normalize vectors (avoid div by zero)
    head_norm = np.linalg.norm(head_vec, axis=1, keepdims=True)
    obj_norm = np.linalg.norm(obj_vec, axis=1, keepdims=True)

    valid = (head_norm[:, 0] > 0) & (obj_norm[:, 0] > 0)
    head_unit = np.zeros_like(head_vec)
    obj_unit = np.zeros_like(obj_vec)
    head_unit[valid] = head_vec[valid] / head_norm[valid]
    obj_unit[valid] = obj_vec[valid] / obj_norm[valid]

    # Angle between head direction and object vector
    cos_sim = np.sum(head_unit * obj_unit, axis=1)
    cos_sim = np.clip(cos_sim, -1.0, 1.0)
    angles_deg = np.degrees(np.arccos(cos_sim))

    # Proximity condition: nose closer to object than body
    dist_nose = np.sqrt((x_nose - x_obj)**2 + (y_nose - y_obj)**2)
    dist_body = np.sqrt((x_body - x_obj)**2 + (y_body - y_obj)**2)
    proximity_ok = dist_nose < dist_body

    # Final condition: angle < threshold AND proximity OK AND vectors valid
    looking_mask = (angles_deg < angle_thresh_deg) & proximity_ok & valid

    return np.where(looking_mask)[0]
