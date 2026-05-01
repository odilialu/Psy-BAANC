# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 13:25:25 2026
cohens kappa
@author: olu
"""
# %% Import packages
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import cohen_kappa_score
# matplotlib plotting parameters
plt.rcParams['xtick.bottom'] = True
plt.rcParams['ytick.left'] = True
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 8
plt.rcParams['axes.titlesize'] = 8
plt.rcParams['axes.labelsize'] = 8
plt.rcParams['axes.titleweight'] = 'regular'
plt.rcParams['figure.dpi'] = 3000
plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['xtick.major.width'] = 0.5
plt.rcParams['ytick.major.width'] = 0.5
plt.rcParams['xtick.major.pad'] = 1.5
plt.rcParams['ytick.major.pad'] = 1.5
plt.rcParams['axes.labelpad'] = 1
plt.rcParams['svg.fonttype'] = 'none'


# %% Variables to change
csvs_path = r"C:/Users/olu/Downloads/FST_cohens_kappa_scoring/Scored Data"
save_path = r"Y:/PsyBAANC/figures/final/kappas"

labs = ["Stanford", "Berkeley1", "Berkeley2", "UCSF1", "UCSF2"]
scorers = ["StanfordSetareh", "StanfordKendall",
           "Berkeley1Odilia",
           "Berkeley2Katrina", "Berkeley2Sam",
           "UCSF1Kelsey", "UCSF1Esther",
           "UCSF2Austin", "UCSF2Ellie", "UCSF2Kira", "UCSF2Lena", "UCSF2Alex"]
vid_names = ["Stanford_Malenka_Setareh", "Berkeley1_Lammel_Odilia", "Berkeley2_Gomez_Katrina",
             "UCSF1_Kheirbek_Kelsey", "UCSF2_Sohal_Lena"]
scorers_short = ["S_S", "S_K",
                 "B1_O",
                 "B2_K", "B2_S",
                 "U1_K", "U1_E",
                 "U2_Au", "U2_E", "U2_K", "U2_L", "U2_Al"]

# %% Get file information and read in csvs.
file_ids = os.listdir(csvs_path)
full_paths = []
for file in file_ids:
    full_paths.append(os.path.join(csvs_path, file))

csvs = {}
for path, file in zip(full_paths, file_ids):
    csvs[file] = pd.read_csv(path)

# %% Load the files and get kappas
kappas_dict = {}
for vid in vid_names:
    kappas_df = pd.DataFrame(index=scorers, columns=scorers)
    files_of_interest_i = [i for i in range(len(file_ids)) if vid in file_ids[i]]
    n_files = len(files_of_interest_i)

    for i in range(n_files):
        for j in range(n_files):
            file1 = file_ids[files_of_interest_i[i]]
            file2 = file_ids[files_of_interest_i[j]]

            scorer1 = file1.split("_")[3]
            scorer2 = file2.split("_")[3]

            csv1 = csvs[file1].to_numpy()
            csv2 = csvs[file2].to_numpy()

            print(f"----Comparing {scorer1} vs. {scorer2} for video: {vid}.")
            print(f"File lengths = {len(csv1)}, {len(csv2)}")

            csv_length = np.min((len(csv1), len(csv2)))
            csv1 = csv1[:csv_length, :]
            csv2 = csv2[:csv_length, :]

            df = pd.DataFrame(
                data={"frame": csv1[:, 0], f"{scorer1}": csv1[:, 1], f"{scorer2}": csv2[:, 1]})
            df_clean = df.dropna()

            kappa_score = cohen_kappa_score(df_clean[scorer1], df_clean[scorer2])

            print(f"Cohen's Kappa: {kappa_score}")

            kappas_df.loc[scorer1, scorer2] = kappa_score
    kappas_dict[vid] = kappas_df

# %%
fig, ax = plt.subplots(nrows=1, ncols=len(vid_names), figsize=(7, 2))
for i, vid in enumerate(vid_names):
    sns.heatmap(kappas_dict[vid].to_numpy().astype(float), cmap='coolwarm',
                ax=ax[i], vmin=0, vmax=1, cbar=True,  xticklabels=scorers_short)
    ax[i].set_yticks([])
    ax[i].set_title(vid)
    plt.tight_layout()

plt.savefig(os.path.join(save_path, "kappas.svg"), transparent=True)

# %% visualization of the raw data.
for vid in vid_names:
    files_of_interest_i = [i for i in range(len(file_ids)) if vid in file_ids[i]]
    file_length = min([len(csvs[file_ids[files_of_interest_i[i]]])
                       for i in range(len(files_of_interest_i))])
    vid_data = pd.DataFrame({"frame": list(range(file_length))})

    for i, file_i in enumerate(files_of_interest_i):
        file = file_ids[file_i]
        csv = csvs[file].to_numpy()[:file_length, 1:].flatten()
        scorer = file.split("_")[3]
        vid_data[scorer] = csv

    vid_data = vid_data.dropna().reset_index()
    vid_data = vid_data[scorers]

    vid_data_resample = pd.DataFrame(index=list(range(240)), columns=scorers)
    indices = np.linspace(0, len(vid_data) - 1, 240).astype(int)
    for scorer in scorers:
        for frame in range(len(indices)-1):
            frame_start = indices[frame]
            frame_end = indices[frame+1]
            vid_data_resample.loc[frame, scorer] = np.sum(
                vid_data.loc[frame_start:frame_end, scorer])
    vid_data_resample = vid_data_resample.astype(float).T

    fig, ax = plt.subplots(figsize=(7, 1.75))
    sns.heatmap(vid_data_resample, cmap='Reds',  # cbar_kws={'ticks': [0, 1]},
                xticklabels=20, yticklabels=True)
    ax.set_title(vid)

    plt.tight_layout()
    plt.savefig(os.path.join(save_path, f"{vid}.svg"), transparent=True)

# %%
for vid in vid_names:
    triu_idx = np.triu_indices_from(kappas_dict[vid], k=1)
    kappas_lab = []
    for x, y in zip(triu_idx[0], triu_idx[1]):
        kappa = kappas_dict[vid].iloc[x, y]
        kappas_lab.append(kappa)
    print(
        f"Cohen's Kappa for {vid}: {round(np.mean(kappas_lab), 2)} ± {round(np.std(kappas_lab, ddof=1), 2)}")
