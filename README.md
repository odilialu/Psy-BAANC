### Usage

1. **Install requirements in a new conda environment:**

```
conda create --name psybaanc python=3.10
```
```
pip install numpy pandas matplotlib opencv-python shapely statsmodels scipy scikit-posthocs openpyxl
````

2. **Format videos:**

    Crop videos using `videoCrop.m` or `videoCropGUI.m` so that only one animal is in the video. Crop out any background activity that may interfere with accurate coordinate tracking. 

3. **Export coordinates:**

    Use DeepLabCut, Biobserve, or Ethovision to export mouse body and nose coordinates as csv from video files.

* Each row is a timepoint.
* For each body part tracked, have one column with X coordinates, another column with Y coordinates.
* The coordinates must be relative to an origin (0,0) that is located at the top left corner of the video. This is the default for DLC but may need to be set for Biobserve / Ethovision. 

4. **Format data:**
   
    Create a folder for each experiment (`OFT_acute`, `OFT_persistent`, `EPM_acute`, `EPM_persistent`, `NOE_acute`, `NOE_persistent`, `SIT_acute_habituation`, `SIT_acute_test`, `SIT_persistent_habituation`, `SIT_persistent_test`).

    Place the following files in each experiment's folder, where `animal_id` is the unique identifier for that animal:

* `{animal_id}.csv` (coordinate file for each animal)
* `{animal_id}.mp4` (video file for each animal)

5. **Run scripts:**
`psybaanc_oft.py`, `psybaanc_epm.py`, `psybaanc_noe.py`, `psybaanc_sit.py`

* Adjust all variables in the "Variables to change" section.
* Run the script and follow all command-line instructions.

6. **Update data:**

* Update supplementary table 1 (raw data) with data found in data_final and data_final_zscored dataframes: https://docs.google.com/spreadsheets/d/1jZRU8h_juDAiZj3U_H3D_dsHmteQ5KVj/edit?usp=sharing&ouid=106513721789887588563&rtpof=true&sd=true. Use Lammel lab's inputs as reference. 
* Update supplementary table 2  (stats) with stats found in stats_results data frame, and the levene's test results and sample sizes printed on the command line: https://docs.google.com/spreadsheets/d/1R5cXtwkwpU6tYqsGRJM9rcBsQPA0_p8K/edit?usp=sharing&ouid=106513721789887588563&rtpof=true&sd=true. Use Lammel lab's inputs as reference. 
* Update Prism files. Make a copy of the Lammel prism file of interest (https://drive.google.com/drive/folders/1S6dTklMii7Y0SjlsDXLQR8yxjmKqcHv-?usp=drive_link), and replace it with your data. Confirm formatting looks the same, and then upload your Prism file to the google drive, in your lab's folder. 
* Notify Odilia when completed. 
