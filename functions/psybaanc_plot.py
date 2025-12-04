# -*- coding: utf-8 -*-
"""
Created on Mon Nov  3 14:45:54 2025

Functions used to create plots of psybaanc data
@author: olu
"""
import cv2
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import patches
import seaborn as sns

matplotlib.use('Agg')

# matplotlib plotting parameters
plt.rcParams['xtick.bottom'] = True
plt.rcParams['ytick.left'] = True
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 8
plt.rcParams['axes.titlesize'] = 8
plt.rcParams['axes.labelsize'] = 8
plt.rcParams['axes.titleweight'] = 'regular'
plt.rcParams['figure.dpi'] = 500
plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['xtick.major.width'] = 0.5
plt.rcParams['ytick.major.width'] = 0.5
plt.rcParams['xtick.major.pad'] = 1.5
plt.rcParams['ytick.major.pad'] = 1.5
plt.rcParams['axes.labelpad'] = 1
plt.rcParams['svg.fonttype'] = 'none'

labs = ["Stanford", "Berkeley 1", "Berkeley 2", "UCSF 1", "UCSF 2"]
custom_palette = dict(zip(labs, ["darkorange", "orangered", "maroon", "purple", "slateblue"]))
markers = dict(zip(labs, ["o", "s", "v", "D", "p"]))

sex_order = ["M", "F"]
treatment_order = ["S", "P"]
stress_order = ["Ctrl", "Stress"]
treatment_colors_bars = ["darkgray", "limegreen"]
treatment_colors_points = ["gray", "green"]


# %%
def add_legend(ax):
    """
    Adds a legend to the given axis.
    Parameters:
    - ax: matplotlib axis to add the legend to.
    """
    text_y = 0.9
    text_x_sal = 0.02
    text_x_psi = 0.30
    rect_y_offset = 0.01
    rect_x_offset = 0.01
    rect_w = 0.28
    rect_h = 0.2
    ax.text(x=text_x_sal, y=text_y, s='Sal', fontsize=8, transform=ax.transAxes)
    ax.text(x=text_x_psi, y=text_y, s='Psi', fontsize=8, transform=ax.transAxes)

    rect_sal = patches.Rectangle((text_x_sal-rect_x_offset, text_y-rect_y_offset),
                                 rect_w, rect_h, lw=0.5, facecolor=treatment_colors_bars[0],
                                 alpha=0.5, zorder=0, transform=ax.transAxes)
    rect_psi = patches.Rectangle((text_x_psi - rect_x_offset, text_y - rect_y_offset),
                                 rect_w, rect_h, facecolor=treatment_colors_bars[1], alpha=0.5,
                                 zorder=0, transform=ax.transAxes)
    ax.add_patch(rect_sal)
    ax.add_patch(rect_psi)


def plot_stars(stars, ax, data):
    """
    Plots significance stars on the given axis.
    Parameters:
    - stars: pd.Series with index as comparison names and values as significance stars.
    - ax: matplotlib axis to plot on.
    - data: pd.Series or np.array with the data used to determine y-axis limits.
    - y_max: Optional maximum y-value for plotting stars. If None, uses data.max().
    """
    # Get the x-coordinates for bars and stars
    bar_positions = []
    star_positions = []
    for patch in ax.patches:
        bar_positions.append(patch.get_x() + patch.get_width() / 2)
        star_positions.append(patch.get_x())
    sorted_bar_positions = sorted(bar_positions)
    sorted_star_positions = sorted(star_positions)

    # Get y-positions for bars and stars.
    ylim_upper = ax.get_ylim()[1]
    star_spacing = (ylim_upper - data.max())/5
    y_star_locs = np.arange(data.max()*1.1, ylim_upper, star_spacing)

    # First check if there are any post-hoc comparisons. If there are, plot those.
    comparisons = ["m,sal v. m,psi", "f,sal v. f,psi", "m,sal v. f,sal", "m,psi v. f,psi"]
    comparisons_plotted = [False]*len(comparisons)
    if comparisons[0] in stars.index.tolist():
        y_coord = y_star_locs[0]  # default y_coord.
        positions = [[0, 1], [2, 3], [0, 2], [1, 3]]
        star_shift = {"*": 0.13, "**": 0.065, "***": 0}
        for i, comparison in enumerate(comparisons):
            if i == 2 and comparisons_plotted[0] is True and comparisons_plotted[1] is True:
                y_coord = y_star_locs[1]
            if i == 3 and comparisons_plotted[1] is True:
                y_coord = y_star_locs[1]
                if comparisons_plotted[2] is True:
                    y_coord = y_star_locs[2]

            if "*" in stars[comparison]:
                comparisons_plotted[i] = True
                x_coords_bar = [sorted_bar_positions[positions[i][0]] + 0.05,
                                sorted_bar_positions[positions[i][1]] - 0.05]
                x_coords_star = [sorted_star_positions[positions[i][0]],
                                 sorted_star_positions[positions[i][1]]]
                y_coords = [y_coord, y_coord]
                # Plot the line over the specified bars
                ax.plot(x_coords_bar, y_coords, color='black', lw=0.5, markersize=0)
                ax.text(
                        x=np.mean(x_coords_star) + star_shift[stars[comparison]],
                        y=y_coords[0]*0.95,
                        s=stars[comparison],
                        fontsize=8,
                        )

    # If there were no posthocs, then plot main effects.
    if np.sum(comparisons_plotted) == 0:  # no comparisons were plotted. Then write main effects.
        if "Kruskal" in stars.index.tolist():
            if "*" in stars["Kruskal"]:
                plot_text = "kruskal" + stars["Kruskal"]
                ax.text(x=0.5, y=0.8, s=plot_text, ha='center', va='top', transform=ax.transAxes)

        elif "Treatment:Sex" in stars.index.tolist():
            if "*" in stars["Treatment:Sex"]:
                plot_text = "drug:sex" + stars["Treatment:Sex"]
                ax.text(x=0.5, y=0.8, s=plot_text, ha='center', va='top', transform=ax.transAxes)

            elif "*" in stars["Treatment"] or "*" in stars["Sex"]:
                if "*" in stars["Treatment"]:
                    plot_text = stars["Treatment"]
                    ax.text(x=0.6, y=0.88, s=stars["Treatment"], transform=ax.transAxes)
                if "*" in stars["Sex"]:
                    ax.set_xlabel(stars["Sex"])


def plot_individual_lab(ax, data_expt_lab, variable, variable_name, lab, stars=None,
                        ymin=None, ymax=None, ybin=None):
    """
    Plots data for an individual lab.
    Parameters:
    - ax: matplotlib axis to plot on.
    - data_expt_lab: pd.DataFrame with data for the specific lab.
    - variable: str, column name of the variable to plot.
    - variable_name: str, name of the variable for labeling.
    - lab: str, name of the lab.
    - stars: pd.Series with significance stars for comparisons.
    - ymin: Optional minimum y-value for the plot. If None, defaults based on variable.
    - ymax: Optional maximum y-value for the plot. If None, defaults based on variable.
    - ybin: Optional y-axis bin size. If None, defaults to half the range.
    """
    ymin, ymax, ybin = y_scale_auto(ymin, ymax, ybin, data_expt_lab[variable], variable_name)
    sns.barplot(data=data_expt_lab,
                x="Sex", order=sex_order,
                y=variable,
                hue="Treatment", hue_order=treatment_order,
                errorbar='se', capsize=0.2, err_kws={'linewidth': 0.5},
                palette=treatment_colors_bars, alpha=0.5,
                legend=False, ax=ax)

    sns.stripplot(data=data_expt_lab,
                  x="Sex", order=sex_order,
                  y=variable,
                  hue="Treatment", hue_order=treatment_order,
                  palette=treatment_colors_points,
                  size=2, linewidth=0.2, alpha=0.75, marker=markers[lab],
                  dodge=True, jitter=0.2,
                  legend=False, ax=ax)

    ax.set_ylim([ymin, ymax])
    ax.set_yticks(list(range(ymin, ymax+1, ybin)))
    ax.set_xlabel("")
    ax.set_xticks(sex_order, ["M", "F"])
    ax.set_title(lab, color=custom_palette[lab])
    ax.set_ylabel(variable_name)

    # then do stats
    if stars is not None:
        plot_stars(stars, ax, data_expt_lab[variable])
    add_legend(ax)

    sns.despine()
    plt.tight_layout()


def plot_over_time(ax, data_binned_dict, col, col_name, time_range):
    """
    Plots line plots over time for males and females.
    Parameters:
    - ax: matplotlib axis to plot on.
    - data_binned_dict: dict of pd.DataFrames with binned data.
    - col: str, column name of the variable to plot.
    - col_name: str, name of the variable for labeling.
    - time_range: tuple, (start_time, end_time) for x-axis labeling.
    """
    data_of_interest = data_binned_dict[col]
    data_males = data_of_interest[data_of_interest["Sex"] == "M"]
    data_females = data_of_interest[data_of_interest["Sex"] == "F"]
    sns.lineplot(data=data_males, x="Interval", y=col, hue="Treatment",
                 palette=treatment_colors_bars, ax=ax, lw=0.5,
                 errorbar='se', err_style='band', err_kws={'lw': 0, 'alpha': 0.2},
                 # err_kws={'linewidth': 0.5, 'capsize':2, 'capthick': 0.5},
                 legend=False)
    sns.lineplot(data=data_females, x="Interval", y=col, hue="Treatment",
                 palette=treatment_colors_bars, ax=ax, ls="--", lw=0.5,
                 errorbar='se', err_style='band', err_kws={'lw': 0, 'alpha': 0.2},
                 # err_kws={'linewidth': 0.5, 'capsize':2, 'capthick': 0.5},
                 legend=False)

    ax.set_xlabel("Time (min)")
    n_intervals = int(len(np.unique(data_of_interest["Interval"])))
    time_bin = int((time_range[1] - time_range[0]) / n_intervals)
    ax.set_xticks(list(range(n_intervals)),
                  list(range(int(time_range[0]), int(time_range[1]), time_bin)))
    ax.set_ylabel(col_name)
    sns.despine()
    plt.tight_layout()


def plot_all_labs(data_expt, variable, variable_name, stars=None,
                  ymin=None, ymax=None, ybin=None, save_path=None):
    """
    Plots data for all labs combined.
    Parameters:
    - data_expt: pd.DataFrame with data for all labs.
    - variable: str, column name of the variable to plot.
    - variable_name: str, name of the variable for labeling.
    - ymin: Optional minimum y-value for the plot. If None, defaults based on variable.
    - ymax: Optional maximum y-value for the plot. If None, defaults based on variable.
    - ybin: Optional y-axis bin size. If None, defaults to 20.
    """
    ymin, ymax, ybin = y_scale_auto(ymin, ymax, ybin, data_expt[variable], variable_name)
    fig, ax = plt.subplots(1, 1, figsize=(1.18, 1.6), sharey=True)
    sns.barplot(data=data_expt,
                x="Sex",
                y=variable,
                hue="Treatment", hue_order=treatment_order,
                errorbar='se', capsize=0.2, err_kws={'linewidth': 0.5},
                palette=treatment_colors_bars, alpha=0.5,
                legend=False, ax=ax)
    for lab_idx, lab in enumerate(labs):
        sns.stripplot(data=data_expt[data_expt["Institution"] == lab],
                      x="Sex", order=sex_order,
                      y=variable,
                      hue="Treatment", hue_order=treatment_order,
                      palette=treatment_colors_points,
                      size=1.35, linewidth=0, alpha=0.75, jitter=0.2, dodge=True,
                      marker=markers[lab],
                      legend=False, ax=ax)

    ax.set_ylim([ymin, ymax])
    ax.set_yticks(list(range(ymin, ymax+1, ybin)))
    ax.set_xlabel("")
    ax.set_xticks(sex_order, ["M", "F"])
    ax.set_title("All labs")
    ax.set_ylabel(variable_name)

    # plot stars
    if stars is not None:
        plot_stars(stars, ax, data_expt[variable])

    # add title and legend
    title_rect = patches.Rectangle((0, 1.05),
                                   1, 0.17, facecolor='white', edgecolor='black', linewidth=0.5,
                                   zorder=0, transform=ax.transAxes, clip_on=False)

    ax.add_patch(title_rect)
    add_legend(ax)

    fig.subplots_adjust(wspace=-1)
    sns.despine()
    plt.tight_layout()
    if save_path is not None:
        plt.savefig(save_path + "summary_results.svg")
        plt.close()


def plot_zscores(psi_data, variable, ymin=-5, ymax=5, ybin=5, save_path=None):
    """
    Plots standardized psi effect data.
    Parameters:
    - psi_data: pd.DataFrame with standardized psi data.
    - variable: str, column name of the standardized variable to plot.
    - ymin: minimum y-value for the plot.
    - ymax: maximum y-value for the plot.
    - ybin: y-axis bin size.
    """
    fig, ax = plt.subplots(figsize=(1.18, 1.5))
    sns.pointplot(
                data=psi_data, x="Sex", y=variable, hue="Institution", hue_order=labs,
                errorbar="se", capsize=0.1, err_kws={'linewidth': 0.5},
                markersize=2, linewidth=0.5,
                palette=custom_palette, linestyle="none", dodge=0.5,
                legend=False
                )
    ax.set_ylabel("z-score")
    ax.set_xlabel("")
    ax.set_ylim([ymin, ymax])
    ax.set_yticks(list(range(0, ymax+1, ybin)))
    ax.axhline(0, lw=0.5, ls='--', color='black')
    sns.despine()
    ax.set_title("Psi Effect")
    plt.tight_layout()
    if save_path is not None:
        plt.savefig(save_path + "/zscores.svg")
        plt.close()


def plot_bars_thirdfactor(ax, data_expt_lab, variable, variable_name, lab,
                          stars=None, ymin=None, ymax=None, ybin=None, save_path=None):
    ymin, ymax, ybin = y_scale_auto(ymin, ymax, ybin, data_expt_lab[variable], variable_name)
    sns.barplot(data=data_expt_lab,
                x="Stress", order=stress_order,
                y=variable,
                hue="Treatment", hue_order=treatment_order,
                errorbar='se', capsize=0.2, err_kws={'linewidth': 0.5},
                palette=treatment_colors_bars, alpha=0.5,
                legend=False, ax=ax)

    sns.stripplot(data=data_expt_lab[data_expt_lab["Sex"] == 'M'],
                  x="Stress", order=stress_order,
                  y=variable,
                  hue="Treatment", hue_order=treatment_order,
                  palette=treatment_colors_points,
                  size=2, linewidth=0.2, alpha=0.75, marker=markers[lab],
                  facecolor='none', edgecolor='auto',
                  dodge=True, jitter=0.2,
                  legend=False, ax=ax)
    # sns.stripplot(data=data_expt_lab[data_expt_lab["Sex"] == 'F'],
    #               x="Stress", order=stress_order,
    #               y=variable,
    #               hue="Treatment", hue_order=treatment_order,
    #               palette=treatment_colors_points,
    #               size=2, linewidth=0.2, alpha=0.75, marker=markers[lab],
    #               dodge=True, jitter=0.2,
    #               legend=False, ax=ax)

    ax.set_ylim([ymin, ymax])
    ax.set_yticks(list(range(ymin, ymax+1, ybin)))
    ax.set_xlabel("")
    ax.set_xticks(stress_order, ["Ctrl", "Stress"])
    ax.set_title(lab, color=custom_palette[lab])
    ax.set_ylabel(variable_name)

    # then do stats
    if stars is not None:
        plot_stars(stars, ax, data_expt_lab[variable])
    add_legend(ax)

    sns.despine()
    plt.tight_layout()


# %% Make plots of animal behavior and rois for sanity check
def plot_traces(ax, video, coordinates, coordinates2=None):
    """
    Plots movement traces on top of the first frame of the video.
    Parameters:
    - ax: matplotlib axis to plot on.
    - video: path to the video file.
    - coordinates: np.array of shape (n, 2) with x and y coordinates
        for the first trace.
    - coordinates2: Optional np.array of shape (n, 2) with x and y coordinates
        for the second trace.
    """
    cap = cv2.VideoCapture(video)
    _, image = cap.read()

    ax.imshow(image)
    ax.plot(coordinates[:, 0], coordinates[:, 1], lw=1, alpha=0.5)
    if coordinates2 is not None:
        ax.plot(coordinates2[:, 0], coordinates2[:, 1], lw=1, alpha=0.5)


def plot_roi_coords(ax, roi):
    """
    Plots the region of interest (ROI) rectangle on the given axis.
    Parameters:
    - ax: matplotlib axis to plot on.
    - roi: tuple of (xmin, xmax, ymin, ymax) defining the ROI.
    """

    # roi = (xmin, xmax, ymin, ymax)
    # show open field base rectangle
    width_roi = roi[1] - roi[0]
    height_roi = roi[3] - roi[2]
    rect = patches.Rectangle((roi[0], roi[2]),
                             width_roi, height_roi,
                             linewidth=1, edgecolor='r', facecolor='none')
    ax.add_patch(rect)


# %% Misc
def y_scale_auto(ymin, ymax, ybin, data, variable_name):
    if ymin is None:
        if "index" in variable_name:
            ymin = int(-1)
        else:
            ymin = int(0)
    if ymax is None:
        if "index" in variable_name:
            ymax = int(1)
        else:
            ymax = int(round(np.max(data)*1.5+5.1, -1))
    if ybin is None:
        ybin = int((ymax-ymin)/2)
    return ymin, ymax, ybin
