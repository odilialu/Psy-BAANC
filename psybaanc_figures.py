# -*- coding: utf-8 -*-
"""
Created on Wed Oct 29 09:49:02 2025
FINAL. 
psybaanc figures
@author: olu
"""

# %% Import packages
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
import functions.psybaanc_stats as psy_stats
import functions.psybaanc_plot as psy_plot

# Dataframe print settings (do not change)
pd.set_option('display.max_rows', None)  # Display all rows
pd.set_option('display.max_columns', None)  # Display all columns
pd.set_option('display.width', 1000)  # Adjust display width to prevent line breaks
pd.set_option('display.max_colwidth', None)  # Display full content of each column

# %% Variables to change
data_path = r"C:\Users\olu\Documents\Psy-BAANC\Paper Drafts\revision materials\Supplementary Table 1_031726.xlsx"
save_path = r"Y:/PsyBAANC/figures/final"
labs = ["Stanford", "Berkeley 1", "Berkeley 2", "UCSF 1", "UCSF 2"]
custom_palette = ["darkorange", "orangered", "maroon", "purple", "slateblue"]
markers = ["o", "s", "v", "D", "p"]
id_cols = ["mouse_ID", "institution", "experiment", "treatment", "sex"]

sex_order = ["M", "F"]
treatment_order = ["S", "P"]
treatment_colors_bars = ["darkgray", "limegreen"]
treatment_colors_points = ["gray", "green"]

pval_limits = [0.05, 0.01, 0.001, 0.0001]
pval_stars = ["ns", "*", "**", "***"]

# matplotlib plotting parameters
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 8
plt.rcParams['axes.titlesize'] = 8  # Set title font size
plt.rcParams['axes.labelsize'] = 8  # Adjust the value as needed
plt.rcParams['axes.titleweight'] = 'regular'
plt.rcParams['figure.dpi'] = 3000
plt.rcParams['axes.linewidth'] = 0.5  # Example: set linewidth to 2.0
plt.rcParams['xtick.major.width'] = 0.5  # Major x-tick width
plt.rcParams['ytick.major.width'] = 0.5  # Major y-tick width
plt.rcParams['xtick.major.size'] = 2  # Length of major ticks on x-axis
plt.rcParams['ytick.major.size'] = 2  # Length of major ticks on y-axis
plt.rcParams['xtick.major.pad'] = 1.5
plt.rcParams['ytick.major.pad'] = 1.5
plt.rcParams['svg.fonttype'] = 'none'

# %% Read in raw data.
data_all = pd.read_excel(data_path, sheet_name=None)

# %% Initialize dicts
stats_results = {}
stats_results_print = {}


# %% Pipelines
def summary_lab_individual(data_expt, variable, stats_results, stats_results_print,
                           ymin=None, ymax=None, ybin=None):
    # First get per-lab figures and stats
    var_id = f"{experiment_name}_{variable}"
    stats_results_labs = {}
    stats_results_labs_print = {}
    fig, ax = plt.subplots(1, len(labs), figsize=(4.2, 1.6), sharey=True)
    for lab_i, lab in enumerate(labs):
        data_expt_lab = data_expt[data_expt["institution"] == lab]
        stats_results_labs[lab], stats_results_labs_print[lab] = psy_stats.stats_treatment_sex(
            data_expt_lab, variable)
        stars = stats_results_labs[lab]["significance"]
        psy_plot.plot_individual_lab(ax[lab_i], data_expt_lab, variable, variable_names[var_i], lab,
                                     stars=stars,
                                     ymin=ymin, ymax=ymax, ybin=ybin)

    stats_results[var_id] = stats_results_labs
    stats_results_print[var_id] = stats_results_labs_print
    plt.savefig(os.path.join(
        save_path, f"{var_id}.svg"), transparent=True)
    plt.savefig(os.path.join(save_path, f"{var_id}.png"))
    plt.close()
    return stats_results, stats_results_print


def summary_lab_all(data_expt, variable, stats_results,
                    ymin=None, ymax=None, ybin=None):
    # Then, Get all lab figure and stats
    var_id = f"{experiment_name}_{variable}_pooled"

    stats_results[var_id], stats_results_print[var_id] = psy_stats.stats_treatment_sex_lab(
        data_expt, variable)
    stars = stats_results[var_id]["significance"]
    fig, ax = plt.subplots(1, 1, figsize=(1.18, 1.7), sharey=True)
    psy_plot.plot_all_labs(ax, data_expt, variable, variable_names[var_i],
                           stars=stars,
                           ymin=ymin, ymax=ymax, ybin=ybin)
    fig.subplots_adjust(wspace=-1)
    plt.savefig(os.path.join(
        save_path, f"{var_id}.svg"), transparent=True)
    plt.savefig(os.path.join(
        save_path, f"{var_id}.png"))
    plt.close()
    return stats_results, stats_results_print


def time_lab_individual(data_expt, variable, stats_results, xticks, xname, xlabel,
                        ymin=None, ymax=None, ybin=None, stats=False):
    # First get per-lab figures and stats
    var_id = f"{experiment_name}_{variable}_time"
    stats_results_labs = {}
    stats_results_labs_print = {}
    fig, ax = plt.subplots(1, len(labs), figsize=(4.2, 1.6), sharey=True)
    for lab_i, lab in enumerate(labs):
        data_expt_lab = data_expt[data_expt["institution"] == lab]
        data_expt_lab = data_expt_lab.dropna(subset=[variable])
        if stats:
            stats_results_labs[lab], stats_results_labs_print[lab] = psy_stats.stats_treatment_sex_third(
                data_expt_lab, variable, xname)
        psy_plot.plot_over_time(ax[lab_i], data_expt_lab, variable, variable_names[var_i],
                                xticks, xname, xlabel,
                                lab, ymin=ymin, ymax=ymax, ybin=ybin)
    plt.savefig(os.path.join(
        save_path, f"{var_id}.svg"), transparent=True)
    plt.savefig(os.path.join(save_path, f"{var_id}.png"))
    plt.close()

    if stats:
        stats_results[var_id] = stats_results_labs
        stats_results_print[var_id] = stats_results_labs_print
        return stats_results, stats_results_print


def time_lab_all(data_expt, variable, stats_results, xticks, xname, xlabel,
                 ymin=None, ymax=None, ybin=None, stats=False):
    # First get per-lab figures and stats
    var_id = f"{experiment_name}_{variable}_time_pooled"
    fig, ax = plt.subplots(1, 1, figsize=(1.18, 1.7), sharey=True)
    psy_plot.plot_over_time(ax, data_expt, variable, variable_names[var_i],
                            xticks, xname, xlabel,
                            ymin=ymin, ymax=ymax, ybin=ybin)
    fig.subplots_adjust(wspace=-1)
    plt.savefig(os.path.join(
        save_path, f"{var_id}.svg"), transparent=True)
    plt.savefig(os.path.join(save_path, f"{var_id}.png"))
    plt.close()

    if stats:
        stats_results[var_id], stats_results_print[var_id] = psy_stats.stats_treatment_sex_third(
            data_expt, variable, xname)
        return stats_results, stats_results_print


def pairs_lab_individual(data_expt, variable, stats_results, paired_var,
                         ymin=None, ymax=None, ybin=None, third_sex=False):
    # First get per-lab figures and stats
    var_id = f"{experiment_name}_{variable}_paired"
    stats_results_labs = {}
    stats_results_labs_print = {}
    fig, ax = plt.subplots(1, len(labs), figsize=(4.2, 1.6), sharey=True)
    for lab_i, lab in enumerate(labs):
        ylabel = lab_i == 0
        data_expt_lab = data_expt[data_expt["institution"] == lab]

        stats_results_labs[lab], stats_results_labs_print[lab] = psy_stats.stats_treatment_sex_third(
            data_expt_lab, variable, paired_var, third_factor_compare=third_sex, third_paired=True)
        stars = stats_results_labs[lab]["significance"]
        psy_plot.plot_pairs(ax[lab_i], data_expt_lab, variable, variable_names[var_i], lab, paired_var,
                            stars=stars,
                            ymin=ymin, ymax=ymax, ybin=ybin, ylabel=ylabel)
    stats_results[var_id] = stats_results_labs
    stats_results_print[var_id] = stats_results_labs_print
    plt.savefig(os.path.join(
        save_path, f"{var_id}.svg"), transparent=True)
    plt.savefig(os.path.join(save_path, f"{var_id}.png"))
    plt.close()
    return stats_results, stats_results_print


def pairs_lab_all(data_expt, variable, stats_results, paired_var,
                  ymin=None, ymax=None, ybin=None, third_sex=False, lab="All labs"):
    # First get per-lab figures and stats
    var_id = f"{experiment_name}_{variable}_paired_pooled"
    fig, ax = plt.subplots(1, 1, figsize=(1.18, 1.7), sharey=True)
    stats_results[var_id], stats_results_print[var_id] = psy_stats.stats_treatment_sex_third(
        data_expt, variable, paired_var, third_factor_compare=third_sex, third_paired=True)
    stars = stats_results[var_id]["significance"]
    psy_plot.plot_pairs(ax, data_expt, variable, variable_names[var_i], lab, paired_var,
                        stars=stars,
                        ymin=ymin, ymax=ymax, ybin=ybin, ylabel=True)
    plt.savefig(os.path.join(save_path, f"{var_id}.svg"), transparent=True)
    plt.savefig(os.path.join(save_path, f"{var_id}.png"))
    plt.close()
    return stats_results, stats_results_print


def stress_lab_individual(data_expt, variable, stats_results, stats_results_print,
                          ymin=None, ymax=None, ybin=None):
    # First get per-lab figures and stats
    var_id = f"{experiment_name}_{variable}"
    stats_results_labs = {}
    stats_results_labs_print = {}
    fig, ax = plt.subplots(1, len(labs), figsize=(4.2, 1.6), sharey=True)
    for lab_i, lab in enumerate(labs):
        data_expt_lab = data_expt[data_expt["institution"] == lab]
        stats_results_labs[lab], stats_results_labs_print[lab] = psy_stats.stats_treatment_sex_stress(
            data_expt_lab, variable)
        stars = stats_results_labs[lab]["significance"]
        psy_plot.plot_bars_thirdfactor(ax[lab_i], data_expt_lab, variable, variable_names[var_i], lab,
                                       stars=stars,
                                       ymin=ymin, ymax=ymax, ybin=ybin)
    stats_results[var_id] = stats_results_labs
    stats_results_print[var_id] = stats_results_labs_print
    plt.savefig(os.path.join(
        save_path, f"{var_id}.svg"), transparent=True)
    plt.savefig(os.path.join(save_path, f"{var_id}.png"))
    plt.close()
    return stats_results, stats_results_print


def stress_lab_all(data_expt, variable, stats_results,
                   ymin=None, ymax=None, ybin=None, lab="All labs"):
    # Then, Get all lab figure and stats
    var_id = f"{experiment_name}_{variable}_pooled"

    stats_results[var_id], stats_results_print[var_id] = psy_stats.stats_treatment_sex_stress(
        data_expt, variable)
    stars = stats_results[var_id]["significance"]
    fig, ax = plt.subplots(1, 1, figsize=(1.18, 1.6), sharey=True)
    for lab_idx, lab in enumerate(labs):
        data_expt_lab = data_expt[data_expt["institution"] == lab]
        psy_plot.plot_bars_thirdfactor(ax, data_expt_lab, variable, variable_names[var_i], lab,
                                       stars=stars,
                                       ymin=ymin, ymax=ymax, ybin=ybin, data_expt=data_expt)
    ax.set_title("All labs", color='black')
    title_rect = patches.Rectangle((0, 1.05),
                                   1, 0.17, facecolor='white', edgecolor='black', linewidth=0.5,
                                   zorder=0, transform=ax.transAxes, clip_on=False)

    ax.add_patch(title_rect)
    fig.subplots_adjust(wspace=-1)
    plt.savefig(os.path.join(
        save_path, f"{var_id}.svg"), transparent=True)
    plt.savefig(os.path.join(
        save_path, f"{var_id}.png"))
    plt.close()
    return stats_results, stats_results_print


def zscore_all(data_expt, variable,
               ymin=-2, ymax=2, ybin=2):
    psi_zscore_data = psy_stats.get_zscores(data_expt, variable)
    fig, ax = plt.subplots(figsize=(1.18, 1.45))
    psy_plot.plot_zscores(ax, psi_zscore_data, variable,
                          ymin=ymin, ymax=ymax, ybin=ybin)
    fig.subplots_adjust(wspace=-1)
    plt.savefig(os.path.join(
        save_path, f"{experiment_name}_{variable}_zscores.svg"), transparent=True)
    plt.savefig(os.path.join(
        save_path, f"{experiment_name}_{variable}_zscores.png"))
    plt.close()
    return psi_zscore_data


def zscore_all_stress(data_expt, variable,
                      ymin=-2, ymax=2, ybin=2):
    psi_zscore_data = psy_stats.get_zscores(data_expt, variable, stress=True)
    psi_zscore_data['sex_stress'] = (psi_zscore_data['stress'] + '\n' + psi_zscore_data['sex'] + ',' +
                                     psi_zscore_data['treatment'].astype(str))
    x_order = list(np.unique(psi_zscore_data['sex_stress']))
    x_order.sort()
    order = [1, 0, 5, 3, 4, 2]
    x_order = [x_order[i] for i in order]
    fig, ax = plt.subplots(figsize=(2, 1.6))
    psy_plot.plot_zscores(ax, psi_zscore_data, variable,
                          ymin=ymin, ymax=ymax, ybin=ybin, x="sex_stress", x_order=x_order)
    fig.subplots_adjust(wspace=-1)
    plt.savefig(os.path.join(
        save_path, f"{experiment_name}_{variable}_zscores.svg"), transparent=True)
    plt.savefig(os.path.join(
        save_path, f"{experiment_name}_{variable}_zscores.png"))
    plt.close()
    return psi_zscore_data


# %% Get plot for head twitch response
experiment_name = "HTR"
variables = ["htr_total"]
variable_names = ["Head twitches (#)"]
data_expt = psy_stats.organize_categories(data_all["HTR"])

data_expt_time = data_expt[id_cols + list(data_expt.columns[-10:])]
data_expt_time = pd.melt(data_expt_time, id_vars=id_cols, value_vars=list(data_expt.columns[-10:]),
                         var_name='Time', value_name='htr_total')
data_expt_time['Time'] = data_expt_time['Time'].str.replace("htr_", '', regex=False).astype(float)


for var_i, variable in enumerate(variables):
    stats_results, stats_results_print = summary_lab_individual(
        data_expt, variable, stats_results, stats_results_print, ymin=0, ymax=90, ybin=30)
    stats_results, stats_results_print = summary_lab_all(data_expt, variable, stats_results,
                                                         ymin=0, ymax=90, ybin=30)

    time_lab_individual(data_expt_time, variable, stats_results,
                        xticks=[0, 10, 20], xname='Time', xlabel="Time (min)",
                        ymin=0, ymax=8, ybin=4)
    time_lab_all(data_expt_time, variable, stats_results,
                 xticks=[0, 10, 20], xname='Time', xlabel="Time (min)",
                 ymin=0, ymax=8, ybin=4)

    psi_zscore_data = zscore_all(data_expt, variable, ymin=0, ymax=40, ybin=20)

# %% Open field test
experiment_names = ["aOFT", "pOFT"]
variables = ["time_center", "time_corners", "distance_in_center", "velocity"]
variable_names = ["Time in center (%)", "Time in corners (%)",
                  "Distance in center (m)", "Velocity (cm/s)"]
data_OFT = data_all["OFT"]
data_OFT = psy_stats.organize_categories(data_OFT)

distance_converted_OFT = False
if distance_converted_OFT is False:
    for col in data_OFT.columns:
        if "distance_in_center" in col:
            data_OFT[col] = data_OFT[col]/100
            distance_converted_OFT = True

ymin = {"aOFT": [0, 10, 0, 0], "pOFT": [0, 10, 0, 0]}
ymax = {"aOFT": [60, 70, 80, 15], "pOFT": [60, 70, 80, 15], "time": [60, 70, 12, 15]}
ybin = {"aOFT": [30, 30, 40, 5], "pOFT": [30, 30, 40, 5], "time": [30, 30, 4, 5]}
ymin_zscore = {"aOFT": [-2, -1, -3, -3], "pOFT": [-2, -2, -3, -3]}
ymax_zscore = {"aOFT": [2, 3, 1, 3], "pOFT": [2, 2, 3, 3]}
ybin_zscore = {"aOFT": [2, 3, 3, 3], "pOFT": [2, 2, 3, 3]}

for experiment_name in experiment_names:
    data_expt = data_OFT[data_OFT["experiment"] == experiment_name]
    for var_i, variable in enumerate(variables):
        cols_time = [col for col in data_expt.columns if f"{variable}_" in col and len(
            col) == len(variable)+2]
        data_expt_time = data_expt[id_cols + cols_time]
        data_expt_time = pd.melt(data_expt_time, id_vars=id_cols, value_vars=cols_time,
                                 var_name='Time', value_name=variable)
        data_expt_time['Time'] = data_expt_time['Time'].str.replace(
            f"{variable}_", '', regex=False).astype(float)
        data_expt_time['Time'] = data_expt_time['Time']*5

        stats_results, stats_results_print = summary_lab_individual(
            data_expt, variable, stats_results, stats_results_print,
            ymin=ymin[experiment_name][var_i], ymax=ymax[experiment_name][var_i], ybin=ybin[experiment_name][var_i])

        stats_results, stats_results_print = summary_lab_all(
            data_expt, variable, stats_results,
            ymin=ymin[experiment_name][var_i], ymax=ymax[experiment_name][var_i], ybin=ybin[experiment_name][var_i])

        time_lab_individual(
            data_expt_time, variable, stats_results,
            xticks=[0, 10, 20, 30], xname='Time', xlabel="Time (min)",
            ymin=ymin[experiment_name][var_i], ymax=ymax["time"][var_i], ybin=ybin["time"][var_i])

        time_lab_all(
            data_expt_time, variable, stats_results,
            xticks=[0, 10, 20, 30], xname='Time', xlabel="Time (min)",
            ymin=ymin[experiment_name][var_i], ymax=ymax["time"][var_i], ybin=ybin["time"][var_i])

        psi_zscore_data = zscore_all(
            data_expt, variable,
            ymin=ymin_zscore[experiment_name][var_i], ymax=ymax_zscore[experiment_name][var_i], ybin=ybin_zscore[experiment_name][var_i])

# %% Elevated Plus Maze
experiment_names = ["aEPM", "pEPM"]
variables = ["time_open", "time_center", "latency_open", "distance"]
variable_names = ["Time in open arms (%)", "Time in center (%)", "Open arm latency (s)", "Distance (m)"
                  ]
data_EPM = data_all["EPM"]
data_EPM = psy_stats.organize_categories(data_EPM)

distance_converted_EPM = False
if distance_converted_EPM is False:
    data_EPM["distance"] = data_EPM["distance"]/100
    distance_converted_EPM = True

ymin = {"aEPM": [0, 0, 0, 0], "pEPM": [0, 0, 0, 0]}
ymax = {"aEPM": [40, 50, 120, 45], "pEPM": [40, 50, 150, 40]}
ybin = {"aEPM": [20, 25, 60, 20], "pEPM": [20, 25, 50, 20]}
ymin_zscore = {"aEPM": [-3, -3, -2, -3], "pEPM": [-2, -3, -3, -3]}
ymax_zscore = {"aEPM": [3, 3, 8, 1], "pEPM": [2, 3, 5, 3]}
ybin_zscore = {"aEPM": [3, 3, 2, 3], "pEPM": [2, 3, 3, 3]}

for experiment_name in experiment_names:
    data_expt = data_EPM[data_EPM["experiment"] == experiment_name]
    for var_i, variable in enumerate(variables):
        stats_results, stats_results_print = summary_lab_individual(
            data_expt, variable, stats_results, stats_results_print,
            ymin=ymin[experiment_name][var_i], ymax=ymax[experiment_name][var_i], ybin=ybin[experiment_name][var_i])

        stats_results, stats_results_print = summary_lab_all(
            data_expt, variable, stats_results,
            ymin=ymin[experiment_name][var_i], ymax=ymax[experiment_name][var_i], ybin=ybin[experiment_name][var_i])

        psi_zscore_data = zscore_all(
            data_expt, variable,
            ymin=ymin_zscore[experiment_name][var_i], ymax=ymax_zscore[experiment_name][var_i], ybin=ybin_zscore[experiment_name][var_i])

# %% Novel object exploration
experiment_names = ["aNOE", "pNOE"]
variables = ["average_time", "latency_explore",
             "time_edges", "distance"]
variable_names = ["Average time (s)", "Latency (s)",
                  "Time in edges (%)", "Distance (m)"]
data_NOE = data_all["NOE"]
data_NOE = psy_stats.organize_categories(data_NOE)

distance_converted_NOE = False
if distance_converted_NOE is False:
    data_NOE["distance"] = data_NOE["distance"]/100
    distance_converted_NOE = True

ymin = {"aNOE": [0, 0, 50, 0], "pNOE": [0, 0, 50, 0]}
ymax = {"aNOE": [45, 50, 110, 80], "pNOE": [45, 50, 110, 80]}
ybin = {"aNOE": [20, 25, 50, 40], "pNOE": [20, 25, 50, 40]}
ymin_zscore = {"aNOE": [-3, -2, -2, -3], "pNOE": [-2, -2, -4, -2]}
ymax_zscore = {"aNOE": [3, 5, 8, 3], "pNOE": [2, 2, 4, 2]}
ybin_zscore = {"aNOE": [3, 2, 2, 3], "pNOE": [2, 2, 4, 2]}

for experiment_name in experiment_names:
    data_expt = data_NOE[data_NOE["experiment"] == experiment_name]
    for var_i, variable in enumerate(variables):
        stats_results, stats_results_print = summary_lab_individual(
            data_expt, variable, stats_results, stats_results_print,
            ymin=ymin[experiment_name][var_i], ymax=ymax[experiment_name][var_i], ybin=ybin[experiment_name][var_i])

        stats_results, stats_results_print = summary_lab_all(
            data_expt, variable, stats_results,
            ymin=ymin[experiment_name][var_i], ymax=ymax[experiment_name][var_i], ybin=ybin[experiment_name][var_i])

        psi_zscore_data = zscore_all(
            data_expt, variable,
            ymin=ymin_zscore[experiment_name][var_i], ymax=ymax_zscore[experiment_name][var_i], ybin=ybin_zscore[experiment_name][var_i])

variable_names = ["Exploration (s)"]
var_i = 0
variables = ["object_one_time", "object_two_time"]
data_expt_wide = data_NOE[id_cols + variables]
data_expt_long = pd.melt(data_expt_wide, id_vars=id_cols, value_vars=variables,
                         var_name='object', value_name="exploration")
data_expt_long.loc[data_expt_long['object'] == 'object_one_time', 'object'] = '1'
data_expt_long.loc[data_expt_long['object'] == 'object_two_time', 'object'] = '2'
for experiment_name in experiment_names:
    data_expt = data_expt_long[data_expt_long["experiment"] == experiment_name]
    stats_results, stats_results_print = pairs_lab_individual(data_expt, "exploration", stats_results,
                                                              "object", ymin=0, ymax=50, ybin=25)
    stats_results, stats_results_print = pairs_lab_all(data_expt, "exploration", stats_results,
                                                       "object", ymin=0, ymax=50, ybin=25)


# %% Social interaction test
experiment_names = ["aSIT", "pSIT"]
variables = ["social_index", "distance"]
variable_names = ["Social index", "Distance (m)"]
data_SIT = data_all["SIT"]
data_SIT = psy_stats.organize_categories(data_SIT)

distance_converted_SIT = False
if distance_converted_SIT is False:
    data_SIT["distance"] = data_SIT["distance"]/100
    data_SIT["distance_chamber_empty"] = data_SIT["distance_chamber_empty"]/100
    data_SIT["distance_chamber_social"] = data_SIT["distance_chamber_social"]/100

    distance_converted_SIT = True

ymin = {"aSIT": [-0.5, 0], "pSIT": [-0.5, 0]}
ymax = {"aSIT": [1.05, 68], "pSIT": [1.05, 68]}
ybin = {"aSIT": [0.5, 30], "pSIT": [0.5, 30]}
ymin_zscore = {"aSIT": [-3, -4], "pSIT": [-2, -3]}
ymax_zscore = {"aSIT": [3, 2], "pSIT": [2, 3]}
ybin_zscore = {"aSIT": [3, 2], "pSIT": [2, 3]}

for experiment_name in experiment_names:
    data_expt = data_SIT[data_SIT["experiment"] == experiment_name]
    for var_i, variable in enumerate(variables):
        stats_results, stats_results_print = summary_lab_individual(
            data_expt, variable, stats_results, stats_results_print,
            ymin=ymin[experiment_name][var_i], ymax=ymax[experiment_name][var_i], ybin=ybin[experiment_name][var_i]
        )

        stats_results, stats_results_print = summary_lab_all(
            data_expt, variable, stats_results,
            ymin=ymin[experiment_name][var_i], ymax=ymax[experiment_name][var_i], ybin=ybin[experiment_name][var_i]
        )

        psi_zscore_data = zscore_all(
            data_expt, variable,
            ymin=ymin_zscore[experiment_name][var_i], ymax=ymax_zscore[experiment_name][var_i], ybin=ybin_zscore[experiment_name][var_i])

variables = [["cup_time_empty", "cup_time_social"],
             ["latency_empty", "latency_social"],
             ["visits_empty", "visits_social"],
             ["chamber_time_empty", "chamber_time_social"]]
variable_names = ["Cup time (s)", "Exploration latency (s)",
                  "Visits (#)", "Chamber time (%)"]

ymin = [0, 0, 0, 0]
ymax = [150, 150, 120, 80]
ybin = [50, 50, 50, 40]

for experiment_name in experiment_names:
    data_expt = data_SIT[data_SIT["experiment"] == experiment_name]
    for var_i, variable_set in enumerate(variables):
        variable = variable_set[0].partition("_e")[0]
        data_expt_wide = data_expt[id_cols + variable_set]
        data_expt_long = pd.melt(data_expt_wide, id_vars=id_cols, value_vars=variable_set,
                                 var_name='zone', value_name=variable)
        data_expt_long['zone'] = data_expt_long['zone'].str.replace(f"{variable}_", '', regex=False)

        stats_results, stats_results_print = pairs_lab_individual(data_expt_long, variable, stats_results,
                                                                  "zone", ymin=ymin[var_i], ymax=ymax[var_i], ybin=ybin[var_i],
                                                                  third_sex=True)
        stats_results, stats_results_print = pairs_lab_all(data_expt_long, variable, stats_results,
                                                           "zone", ymin=ymin[var_i], ymax=ymax[var_i], ybin=ybin[var_i], third_sex=True)

# %% TST in naive mice
experiment_names = ["pTST"]
variables = ["time_immobile"]
variable_names = ["Immobility (s)"]
data_TST = data_all["pTST"]
data_TST = psy_stats.organize_categories(data_TST)
ymin = [0]
ymax = [400]
ybin = [100]

for experiment_name in experiment_names:
    data_expt = data_TST[data_TST["experiment"] == experiment_name]
    for var_i, variable in enumerate(variables):
        stats_results, stats_results_print = summary_lab_individual(
            data_expt, variable, stats_results, stats_results_print,
            ymin=ymin[var_i], ymax=ymax[var_i], ybin=ybin[var_i])

        stats_results, stats_results_print = summary_lab_all(
            data_expt, variable, stats_results,
            ymin=ymin[var_i], ymax=ymax[var_i], ybin=ybin[var_i])

        psi_zscore_data = zscore_all(data_expt, variable, ymin=-4, ymax=4, ybin=4)

# %% Repeated Forced swim test and tail suspension test expts.
experiment_names = ["rFST"]
variables = ["day_3_immobility"]
variable_names = ["Day 3 immobility (s)"]
data_FST = data_all["rFST"]
data_FST = psy_stats.organize_categories(data_FST)
ymin = [0]
ymax = [300]
ybin = [100]

for experiment_name in experiment_names:
    data_expt = data_FST[data_FST["experiment"] == experiment_name]
    for var_i, variable in enumerate(variables):
        stats_results, stats_results_print = summary_lab_individual(
            data_expt, variable, stats_results, stats_results_print,
            ymin=ymin[var_i], ymax=ymax[var_i], ybin=ybin[var_i])

        stats_results, stats_results_print = summary_lab_all(
            data_expt, variable, stats_results,
            ymin=ymin[var_i], ymax=ymax[var_i], ybin=ybin[var_i])

        psi_zscore_data = zscore_all(data_expt, variable, ymin=-4, ymax=4, ybin=4)

variable_names = ["Immobility (s)"]
var_i = 0
variables = ["day_1_immobility", "day_3_immobility"]
data_expt_wide = data_FST[id_cols + variables]
data_expt_long = pd.melt(data_expt_wide, id_vars=id_cols, value_vars=variables,
                         var_name='day', value_name="immobility")
data_expt_long.loc[data_expt_long['day'] == 'day_1_immobility', 'day'] = '1'
data_expt_long.loc[data_expt_long['day'] == 'day_3_immobility', 'day'] = '3'
for experiment_name in experiment_names:
    data_expt = data_expt_long[data_expt_long["experiment"] == experiment_name]
    stats_results, stats_results_print = pairs_lab_individual(data_expt, "immobility", stats_results,
                                                              "day", ymin=0, ymax=300, ybin=100, third_sex=True)
    stats_results, stats_results_print = pairs_lab_all(data_expt, "immobility", stats_results,
                                                       "day", ymin=0, ymax=300, ybin=100, third_sex=True)

# %%  Tail suspension test Repeat expt
experiment_names = ["rTST"]
variables = ["day_3_immobility"]
variable_names = ["Day 3 immobility (s)"]
data_TST = data_all["rTST"]
data_TST = psy_stats.organize_categories(data_TST)
ymin = [0]
ymax = [400]
ybin = [100]

for experiment_name in experiment_names:
    data_expt = data_TST[data_TST["experiment"] == experiment_name]
    for var_i, variable in enumerate(variables):

        stats_results, stats_results_print = summary_lab_individual(
            data_expt, variable, stats_results, stats_results_print,
            ymin=ymin[var_i], ymax=ymax[var_i], ybin=ybin[var_i])

        stats_results, stats_results_print = summary_lab_all(
            data_expt, variable, stats_results,
            ymin=ymin[var_i], ymax=ymax[var_i], ybin=ybin[var_i])

        psi_zscore_data = zscore_all(data_expt, variable, ymin=-3, ymax=3, ybin=3)


variable_names = ["Immobility (s)"]
var_i = 0
variables = ["day_1_immobility", "day_3_immobility"]
data_expt_wide = data_TST[id_cols + variables]
data_expt_long = pd.melt(data_expt_wide, id_vars=id_cols, value_vars=variables,
                         var_name='day', value_name="immobility")
data_expt_long.loc[data_expt_long['day'] == 'day_1_immobility', 'day'] = '1'
data_expt_long.loc[data_expt_long['day'] == 'day_3_immobility', 'day'] = '3'
for experiment_name in experiment_names:
    data_expt = data_expt_long[data_expt_long["experiment"] == experiment_name]
    stats_results, stats_results_print = pairs_lab_individual(data_expt, "immobility", stats_results,
                                                              "day", ymin=0, ymax=400, ybin=100, third_sex=True)
    stats_results, stats_results_print = pairs_lab_all(data_expt, "immobility", stats_results,
                                                       "day", ymin=0, ymax=400, ybin=100, third_sex=True)


# %% cort EPM
experiment_name = "cortEPM"
variables = ["time_open", "time_center", "latency_open", "distance"]
variable_names = ["Time in open arms (%)", "Time in center (%)", "Open arm latency (s)", "Distance (m)"
                  ]
data_cortEPM = data_all["CORT (EPM)"]
data_cortEPM = psy_stats.organize_categories(data_cortEPM)
data_expt = data_cortEPM

distance_converted_cortEPM = False
if distance_converted_cortEPM is False:
    data_cortEPM["distance"] = data_cortEPM["distance"]/100
    distance_converted_cortEPM = True

ymin = [0]*len(variables)
ymax = [35, 35, 150, 40]
ybin = [10, 10, 50, 20]
ymin_zscore = [-2, -2.5, -2, -2.5]
ymax_zscore = [3, 2.5, 13, 2.5]
ybin_zscore = [2, 2.5, 2, 2.5]

for var_i, variable in enumerate(variables):
    stats_results, stats_results_print = stress_lab_individual(
        data_expt, variable, stats_results, stats_results_print,
        ymin=ymin[var_i], ymax=ymax[var_i], ybin=ybin[var_i])

    stats_results, stats_results_print = stress_lab_all(
        data_expt, variable, stats_results,
        ymin=ymin[var_i], ymax=ymax[var_i], ybin=ybin[var_i])

    psi_zscore_data = zscore_all_stress(data_expt, variable,
                                        ymin=ymin_zscore[var_i], ymax=ymax_zscore[var_i], ybin=ybin_zscore[var_i])


# %% cort SPT
experiment_name = "cortSPT"
variables = ["sucrose", "water", "sucrose_pref_index", "total"]
variable_names = ["Sucrose (g)", "Water (g)", "Sucrose pref", "Total (g)"]
data_cortSPT = data_all["CORT (SPT)"]
data_cortSPT = psy_stats.organize_categories(data_cortSPT)
data_expt = data_cortSPT

ymin = [0]*len(variables)
ymax = [8, 3, 1.35, 8]
ybin = [4, 1, 0.5, 4]

ymin_zscore = [-4, -3, -6, -5]
ymax_zscore = [2, 4, 2, 2]
ybin_zscore = [2, 3, 2, 5]

for var_i, variable in enumerate(variables):
    stats_results, stats_results_print = stress_lab_individual(
        data_expt, variable, stats_results, stats_results_print,
        ymin=ymin[var_i], ymax=ymax[var_i], ybin=ybin[var_i])

    stats_results, stats_results_print = stress_lab_all(
        data_expt, variable, stats_results,
        ymin=ymin[var_i], ymax=ymax[var_i], ybin=ybin[var_i])

    psi_zscore_data = zscore_all_stress(
        data_expt, variable,
        ymin=ymin_zscore[var_i], ymax=ymax_zscore[var_i], ybin=ybin_zscore[var_i])


# %% cort FST
experiment_name = "cortFST"
variables = ["time_immobile", "struggling_bouts"]
variable_names = ["Immobility (s)", "Struggling bouts (#)"]
data_cortFST = data_all["CORT (FST)"]
data_cortFST = psy_stats.organize_categories(data_cortFST)
data_expt = data_cortFST

ymin = [0]*len(variables)
ymax = [300, 40]
ybin = [100, 25]
ymin_zscore = [-4, -3]
ymax_zscore = [2, 3]
ybin_zscore = [2, 3]

for var_i, variable in enumerate(variables):
    stats_results, stats_results_print = stress_lab_individual(
        data_expt, variable, stats_results, stats_results_print,
        ymin=ymin[var_i], ymax=ymax[var_i], ybin=ybin[var_i])

    stats_results, stats_results_print = stress_lab_all(
        data_expt, variable, stats_results,
        ymin=ymin[var_i], ymax=ymax[var_i], ybin=ybin[var_i])

    psi_zscore_data = zscore_all_stress(
        data_expt, variable,
        ymin=ymin_zscore[var_i], ymax=ymax_zscore[var_i], ybin=ybin_zscore[var_i])


# %% Fear conditioning
experiment_names = ["pre_retrieval", "post_retrieval"]
session_names = ['encoding_freezing', 'retrieval_freezing',
                 'extinction1_freezing', 'extinction2_freezing']
data_FC = data_all["FC"]
data_FC = psy_stats.organize_categories(data_FC)


for experiment_name in experiment_names:
    var_id = f"FC_{experiment_name}"
    data_expt = data_FC[data_FC["experiment"] == experiment_name]
    variable_names = data_expt.columns[-31:].tolist()

    data_final = data_expt[id_cols + variable_names]
    df_long = pd.wide_to_long(data_final,
                              stubnames=session_names,
                              i=id_cols,
                              j='cue',
                              sep='_').reset_index()

    # individual labs figure
    fig, ax = plt.subplots(4, 5, figsize=(5, 5.5), sharey=True)
    for session_i, session in enumerate(session_names):
        stats_results_labs = {}
        stats_results_print_labs = {}
        stats_results_labs_b = {}
        stats_results_print_labs_b = {}
        for lab_i, lab in enumerate(labs):
            df_lab = df_long[df_long["institution"] == lab]
            df_session = df_lab.dropna(subset=[session])
            df_baseline_stats = df_session[df_session["cue"] == 0]
            df_baseline_stats.rename(columns={f'{session}': f'{session}_b'}, inplace=True)
            df_cue_stats = df_session[df_session["cue"] != 0]
            if session == "retrieval_freezing":
                df_cue_stats["cue"] = df_cue_stats["cue"].astype(str)
            stats_results_labs_b[lab], stats_results_print_labs_b[lab] = psy_stats.stats_treatment_sex(
                df_baseline_stats, f"{session}_b")
            stats_results_labs[lab], stats_results_print_labs[lab] = psy_stats.stats_treatment_sex_third(
                df_cue_stats, session, third_factor="cue", third_factor_compare=False, third_paired=True)
            max_cue = df_session["cue"].max()
            psy_plot.plot_over_time(ax[session_i, lab_i],
                                    df_session, session, "Freezing (%)",
                                    [1, max_cue], x_name="cue", xlabel="Cue (#)",
                                    ymin=0, ymax=100, ybin=50, lab=lab)
        stats_results[f"{var_id}_{session}"] = stats_results_labs
        stats_results[f"{var_id}_{session}_b"] = stats_results_labs_b
        stats_results_print[f"{var_id}_{session}"] = stats_results_print_labs
        stats_results_print[f"{var_id}_{session}_b"] = stats_results_print_labs_b
    for a in ax.flat:
        a.label_outer()
    plt.tight_layout()
    plt.savefig(os.path.join(save_path, f"{var_id}_time.svg"), transparent=True)
    plt.savefig(os.path.join(save_path, f"{var_id}_time.png"), transparent=True)
    plt.close()

    # all labs figure
    fig, ax = plt.subplots(4, 1, figsize=(1.3, 5.5), sharey=True)
    for session_i, session in enumerate(session_names):
        df_session = df_long.dropna(subset=[session])
        df_baseline_stats = df_session[df_session["cue"] == 0]
        df_baseline_stats.rename(columns={f'{session}': f'{session}_b'}, inplace=True)
        df_cue_stats = df_session[df_session["cue"] != 0]
        df_cue_stats["cue"] = pd.to_numeric(df_cue_stats["cue"])
        if session == "retrieval_freezing":
            df_cue_stats["cue"] = df_cue_stats["cue"].astype(str)
        stats_results[f"{var_id}_{session}_b_pooled"], stats_results_print[f"{var_id}_{session}_b_pooled"] = psy_stats.stats_treatment_sex_lab(
            df_baseline_stats, f"{session}_b")
        try:
            stats_results[f"{var_id}_{session}_pooled"], stats_results_print[f"{var_id}_{session}_pooled"] = psy_stats.stats_treatment_sex_third(
                df_cue_stats, session, third_factor="cue", third_factor_compare=False, third_paired=True)
        except:
            stats_results[f"{var_id}_{session}_pooled"], stats_results_print[f"{var_id}_{session}_pooled"] = [
                "Error"], ["Error"]
        max_cue = df_session["cue"].max()
        psy_plot.plot_over_time(ax[session_i],
                                df_session, session, "Freezing (%)",
                                [1, max_cue], x_name="cue", xlabel="Cue (#)",
                                ymin=0, ymax=100, ybin=50)
    for a in ax.flat:
        a.label_outer()
    plt.tight_layout()
    plt.savefig(os.path.join(save_path, f"{var_id}_time_pooled.svg"), transparent=True)
    plt.savefig(os.path.join(save_path, f"{var_id}_time_pooled.png"), transparent=True)

    plt.close()

# %% fear conditioning summary
experiment_names = ["pre_retrieval", "post_retrieval"]
variables = ["encoding_freezing_avg",
             "retrieval_freezing_avg", "extinction1_freezing_avg", "extinction2_freezing_avg"]
variable_names = ["Avg freezing (%)"]*4
data_FC = data_all["FC"]
data_FC = psy_stats.organize_categories(data_FC)

ymin = [0]*4
ymax = [120]*4
ybin = [50]*4
ymin_zscore = {"pre_retrieval": [-2, -4, -2, -2], "post_retrieval": [-2, -2, -3, -3]}
ymax_zscore = {"pre_retrieval": [2, 1, 2, 2], "post_retrieval": [2, 2, 3, 3]}
ybin_zscore = {"pre_retrieval": [2, 2, 2, 2], "post_retrieval": [2, 2, 3, 3]}

for experiment_name in experiment_names:
    data_expt = data_FC[data_FC["experiment"] == experiment_name]
    for var_i, variable in enumerate(variables):
        data_expt = data_expt.dropna(subset=variable)
        stats_results, stats_results_print = summary_lab_individual(
            data_expt, variable, stats_results, stats_results_print,
            ymin=ymin[var_i], ymax=ymax[var_i], ybin=ybin[var_i])

        stats_results, stats_results_print = summary_lab_all(
            data_expt, variable, stats_results,
            ymin=ymin[var_i], ymax=ymax[var_i], ybin=ybin[var_i])

        psi_zscore_data = zscore_all(
            data_expt, variable,
            ymin=ymin_zscore[experiment_name][var_i], ymax=ymax_zscore[experiment_name][var_i], ybin=ybin_zscore[experiment_name][var_i])


# %% sCPP
# juveniles vs. adults.
experiment_name = "sCPP_JA"
data_sCPP = data_all["sCPP"]

data_sCPP_JA = data_sCPP[data_sCPP["treatment"] != "P"]
data_sCPP_JA = data_sCPP_JA.rename(
    columns={"treatment": "treatment_NAN", "experiment": "treatment"})
data_sCPP_JA["treatment"] = data_sCPP_JA["treatment"].map({"juv": "S", "adult": "P"})
data_sCPP_JA = psy_stats.organize_categories(data_sCPP_JA)

variable_names = ["Time spent (%)"]
var_i = 0
variables = ["time_social_side_pre", "time_social_side_post"]
id_cols_sCPP = ['mouse_ID', 'institution', 'treatment', 'sex']
data_expt_wide = data_sCPP_JA[id_cols_sCPP + variables]
data_expt_long = pd.melt(data_expt_wide, id_vars=id_cols_sCPP, value_vars=variables,
                         var_name='day_scpp', value_name="time_social")
data_expt_long['day_scpp'] = data_expt_long['day_scpp'].str.replace(
    "time_social_side_", '', regex=False)
data_expt_long["day_scpp"] = data_expt_long["day_scpp"].map({"pre": "1", "post": "5"})

for lab in labs:
    data_temp = data_expt_long[data_expt_long["institution"] == lab]
    stats_results[f"sCPP_JvA_prepost_{lab}"], stats_results_print[f"sCPP_JvA_prepost_{lab}"] = psy_stats.paired_t(
        data_temp, "time_social", "day_scpp", group_names=["Juv", "Adult"])
stats_results["sCPP_JvA_prepost_pooled"], stats_results_print["sCPP_JvA_prepost_pooled"] = psy_stats.paired_t(
    data_expt_long, "time_social", "day_scpp", group_names=["Juv", "Adult"])
stats_results, stats_results_print = pairs_lab_individual(data_expt_long, "time_social", stats_results,
                                                          "day_scpp", ymin=0, ymax=100, ybin=50, third_sex=True)

stats_results, stats_results_print = pairs_lab_all(data_expt_long, "time_social", stats_results,
                                                   "day_scpp", ymin=0, ymax=100, ybin=50, third_sex=True)

variable_names = ["Social preference"]
# pref score
stats_results, stats_results_print = summary_lab_individual(
    data_sCPP_JA, "social_preference_score", stats_results, stats_results_print,
    ymin=0, ymax=2.5, ybin=1)

stats_results, stats_results_print = summary_lab_all(
    data_sCPP_JA, "social_preference_score", stats_results,
    ymin=0, ymax=2.5, ybin=1)

psi_zscore_data = zscore_all(
    data_sCPP_JA, "social_preference_score",
    ymin=-3, ymax=1.5, ybin=3)

# %%
# Adults sal v. psi.
experiment_name = "sCPP_adults"
data_sCPP_adults = data_sCPP[data_sCPP["experiment"] == "adult"]
data_sCPP_adults = psy_stats.organize_categories(data_sCPP_adults)
variable_names = ["Time spent (%)"]
var_i = 0
variables = ["time_social_side_pre", "time_social_side_post"]
data_expt_wide = data_sCPP_adults[id_cols + variables]
data_expt_long = pd.melt(data_expt_wide, id_vars=id_cols, value_vars=variables,
                         var_name='day_scpp', value_name="time_social")
data_expt_long['day_scpp'] = data_expt_long['day_scpp'].str.replace(
    "time_social_side_", '', regex=False)
data_expt_long["day_scpp"] = data_expt_long["day_scpp"].map({"pre": "1", "post": "5"})

for lab in labs:
    data_temp = data_expt_long[data_expt_long["institution"] == lab]
    stats_results[f"sCPP_SvP_prepost_{lab}"], stats_results_print[f"sCPP_SvP_prepost_{lab}"] = psy_stats.paired_t(
        data_temp, "time_social", "day_scpp")
stats_results["sCPP_SvP_prepost_pooled"], stats_results_print["sCPP_SvP_prepost_pooled"] = psy_stats.paired_t(
    data_expt_long, "time_social", "day_scpp")
stats_results, stats_results_print = pairs_lab_individual(data_expt_long, "time_social", stats_results,
                                                          "day_scpp", ymin=0, ymax=100, ybin=50, third_sex=True)

stats_results, stats_results_print = pairs_lab_all(data_expt_long, "time_social", stats_results,
                                                   "day_scpp", ymin=0, ymax=100, ybin=50, third_sex=True)
# pref score
variable_names = ["Social preference"]
stats_results, stats_results_print = summary_lab_individual(
    data_sCPP_adults, "social_preference_score", stats_results, stats_results_print,
    ymin=0, ymax=2.5, ybin=1)

stats_results, stats_results_print = summary_lab_all(
    data_sCPP_adults, "social_preference_score", stats_results,
    ymin=0, ymax=2.5, ybin=1)

psi_zscore_data = zscore_all(
    data_sCPP_adults, "social_preference_score",
    ymin=-2, ymax=2, ybin=2)


# %% Write out statistics table:

def write_content(stats_results, keys, ws_name, startrow=0, empty_rows_between_dfs=3):
    output_file = r"C:\Users\olu\Documents\output_dataframes.xlsx"
    with pd.ExcelWriter(output_file, engine='openpyxl', mode='a', if_sheet_exists='overlay') as writer:
        for key in keys:
            results_to_write = stats_results[key]
            if type(results_to_write) is dict:
                for lab in labs:
                    labels_to_write = stats_results_print[key][lab]
                    results_to_write = stats_results[key][lab]
                    for label in labels_to_write:
                        pd.DataFrame({label: [label]}).to_excel(writer, sheet_name=ws_name,
                                                                index=False, header=False,
                                                                startrow=startrow)
                        startrow += 1

                    # Write dataframe
                    results_to_write.to_excel(writer, sheet_name=ws_name,
                                              startrow=startrow, index=True)
                    startrow += len(results_to_write) + empty_rows_between_dfs

            else:
                labels_to_write = stats_results_print[key]
                for label in labels_to_write:
                    pd.DataFrame({label: [label]}).to_excel(writer, sheet_name=ws_name,
                                                            index=False, header=False, startrow=startrow)
                    startrow += 1

                # Write dataframe
                results_to_write.to_excel(writer, sheet_name=ws_name, startrow=startrow, index=True)
                startrow += len(results_to_write) + empty_rows_between_dfs


# %% Main figs
write_content(stats_results, ["aOFT_time_center", "aOFT_time_center_pooled",
                              "aEPM_time_open", "aEPM_time_open_pooled",
                              "pOFT_time_center", "pOFT_time_center_pooled",
                              "pEPM_time_open", "pEPM_time_open_pooled"], "Fig 1")

write_content(stats_results, ["aNOE_average_time", "aNOE_average_time_pooled",
                              "pNOE_average_time", "pNOE_average_time_pooled"], "Fig 2")

write_content(stats_results, ["aSIT_cup_time_paired", "aSIT_cup_time_paired_pooled",
                              "aSIT_social_index", "aSIT_social_index_pooled",
                              "pSIT_cup_time_paired", "pSIT_cup_time_paired_pooled",
                              "pSIT_social_index", "pSIT_social_index_pooled"], "Fig 3")

write_content(stats_results, ["pTST_time_immobile", "pTST_time_immobile_pooled",
                              "rFST_immobility_paired", "rFST_immobility_paired_pooled",
                              "rFST_day_3_immobility", "rFST_day_3_immobility_pooled",
                              "rTST_immobility_paired", "rTST_immobility_paired_pooled",
                              "rTST_day_3_immobility", "rTST_day_3_immobility_pooled"], "Fig 4")

write_content(stats_results, ["cortEPM_time_center", "cortEPM_time_center_pooled",
                              "cortEPM_distance", "cortEPM_distance_pooled",
                              "cortSPT_sucrose_pref_index", "cortSPT_sucrose_pref_index_pooled",
                              "cortFST_time_immobile", "cortFST_time_immobile_pooled"], "Fig 5")

write_content(stats_results, ["FC_pre_retrieval_encoding_freezing_b", "FC_pre_retrieval_encoding_freezing",
                              "FC_pre_retrieval_retrieval_freezing_b", "FC_pre_retrieval_retrieval_freezing",
                              "FC_pre_retrieval_extinction1_freezing_b", "FC_pre_retrieval_extinction1_freezing",
                              "FC_pre_retrieval_extinction2_freezing_b", "FC_pre_retrieval_extinction2_freezing",
                              "FC_pre_retrieval_encoding_freezing_b_pooled", "FC_pre_retrieval_encoding_freezing_pooled",
                              "FC_pre_retrieval_retrieval_freezing_b_pooled", "FC_pre_retrieval_retrieval_freezing_pooled",
                              "FC_pre_retrieval_extinction1_freezing_b_pooled", "FC_pre_retrieval_extinction1_freezing_pooled",
                              "FC_pre_retrieval_extinction2_freezing_b_pooled", "FC_pre_retrieval_extinction2_freezing_pooled",

                              "FC_post_retrieval_encoding_freezing_b", "FC_post_retrieval_encoding_freezing",
                              "FC_post_retrieval_retrieval_freezing_b", "FC_post_retrieval_retrieval_freezing",
                              "FC_post_retrieval_extinction1_freezing_b", "FC_post_retrieval_extinction1_freezing",
                              "FC_post_retrieval_extinction2_freezing_b", "FC_post_retrieval_extinction2_freezing",
                              "FC_post_retrieval_encoding_freezing_b_pooled", "post_retrieval_encoding_freezing_avg_pooled",
                              "FC_post_retrieval_retrieval_freezing_b_pooled", "FC_post_retrieval_retrieval_freezing_pooled",
                              "FC_post_retrieval_extinction1_freezing_b_pooled", "FC_post_retrieval_extinction1_freezing_pooled",
                              "FC_post_retrieval_extinction2_freezing_b_pooled", "FC_post_retrieval_extinction2_freezing_pooled",
                              ], "Fig 6")

write_content(stats_results, ["sCPP_JvA_prepost_Stanford", "sCPP_JvA_prepost_Berkeley 1",
                              "sCPP_JvA_prepost_Berkeley 2", "sCPP_JvA_prepost_UCSF 1",
                              "sCPP_JvA_prepost_UCSF 2", "sCPP_JvA_prepost_pooled",
                              "sCPP_JA_social_preference_score", "sCPP_JA_social_preference_score_pooled",

                              "sCPP_SvP_prepost_Stanford", "sCPP_SvP_prepost_Berkeley 1",
                              "sCPP_SvP_prepost_Berkeley 2", "sCPP_SvP_prepost_UCSF 1",
                              "sCPP_SvP_prepost_UCSF 2", "sCPP_SvP_prepost_pooled",
                              "sCPP_adults_social_preference_score", "sCPP_adults_social_preference_score_pooled",
                              ], "Fig 7")

# %% ED
write_content(stats_results, ["HTR_htr_total", "HTR_htr_total_pooled"], "ED1")

write_content(stats_results, ["aOFT_time_corners_pooled", "aOFT_distance_in_center_pooled",
                              "pOFT_time_corners_pooled", "pOFT_distance_in_center_pooled"], "ED2")

write_content(stats_results, ["aOFT_velocity", "aOFT_velocity_pooled",
                              "pOFT_velocity", "pOFT_velocity_pooled"], "ED3")

write_content(stats_results, ["aEPM_time_center", "aEPM_time_center_pooled",
                              "aEPM_distance", "aEPM_distance_pooled",
                              "aEPM_latency_open", "aEPM_latency_open_pooled",
                              "pEPM_time_center", "pEPM_time_center_pooled",
                              "pEPM_distance", "pEPM_distance_pooled",
                              "pEPM_latency_open", "pEPM_latency_open_pooled"], "ED4")

write_content(stats_results, ["aNOE_exploration_paired", "aNOE_exploration_paired_pooled",
                              "aNOE_latency_explore", "aNOE_latency_explore_pooled",
                              "aNOE_distance", "aNOE_distance_pooled",
                              "aNOE_time_edges", "aNOE_time_edges_pooled",
                              "pNOE_exploration_paired", "pNOE_exploration_paired_pooled",
                              "pNOE_latency_explore", "pNOE_latency_explore_pooled",
                              "pNOE_distance", "pNOE_distance_pooled",
                              "pNOE_time_edges", "pNOE_time_edges_pooled"], "ED5")

write_content(stats_results, ["aSIT_visits_paired", "aSIT_visits_paired_pooled",
                              "aSIT_latency_paired", "aSIT_latency_paired_pooled",
                              "aSIT_chamber_time_paired", "aSIT_chamber_time_paired_pooled",
                              "aSIT_distance", "aSIT_distance_pooled",
                              "pSIT_visits_paired", "pSIT_visits_paired_pooled",
                              "pSIT_latency_paired", "pSIT_latency_paired_pooled",
                              "pSIT_chamber_time_paired", "pSIT_chamber_time_paired_pooled",
                              "pSIT_distance", "pSIT_distance_pooled",], "ED6")

write_content(stats_results, ["cortEPM_time_open", "cortEPM_time_open_pooled",
                              "cortEPM_latency_open", "cortEPM_latency_open_pooled",
                              "cortSPT_total", "cortSPT_total_pooled",
                              "cortSPT_sucrose", "cortSPT_sucrose_pooled",
                              "cortSPT_water", "cortSPT_water_pooled",
                              "cortFST_struggling_bouts", "cortFST_struggling_bouts_pooled"], "ED8")
