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
import functions.psybaanc_stats as psy_stats

matplotlib.use('Agg')

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

labs = ["Stanford", "Berkeley 1", "Berkeley 2", "UCSF 1", "UCSF 2"]
custom_palette = dict(
    zip(labs, ["darkorange", "orangered", "maroon", "purple", "slateblue"]))
markers = dict(zip(labs, ["o", "s", "v", "D", "p"]))

sex_order = ["M", "F"]
treatment_order = ["S", "P"]
cup_order = ["empty", "social"]
stress_order = ["Ctrl", "Stress"]
treatment_colors_bars = ["darkgray", "limegreen"]  # ["bisque", "darkgray"]
treatment_colors_points = ["gray", "green"]  # ["orange", "gray"]


# %%
def add_legend(ax, alpha=0.5, scpp=False):
    """
    Adds a legend to the given axis.
    Parameters:
    - ax: matplotlib axis to add the legend to.
    """
    text_y = 0.9
    text_x_sal = 0.02
    text_x_psi = 0.30  # 0.30
    if scpp:
        text_x_psi = 0.32
    rect_y_offset = 0.01
    rect_x_offset = 0.01
    rect_w = 0.28
    rect_w_psi = 0.28
    if scpp:
        rect_w = 0.3
        rect_w_psi = 0.41
    rect_h = 0.2
    ax.text(x=text_x_sal, y=text_y, s='Sal',  # s='Juv'
            fontsize=8, transform=ax.transAxes)
    ax.text(x=text_x_psi, y=text_y, s='Psi',  # ,s='Adult'
            fontsize=8, transform=ax.transAxes)

    rect_sal = patches.Rectangle((text_x_sal-rect_x_offset, text_y-rect_y_offset),
                                 rect_w, rect_h, lw=0.5, facecolor=treatment_colors_bars[0],
                                 alpha=alpha, zorder=0, transform=ax.transAxes)
    rect_psi = patches.Rectangle((text_x_psi - rect_x_offset, text_y - rect_y_offset),
                                 rect_w_psi, rect_h, facecolor=treatment_colors_bars[1], alpha=alpha,
                                 zorder=0, transform=ax.transAxes)
    ax.add_patch(rect_sal)
    ax.add_patch(rect_psi)


def plot_stars(stars, ax, data, third=None):
    """
    Plots significance stars on the given axis.
    Parameters:
    - stars: pd.Series with index as comparison names and values as significance stars.
    - ax: matplotlib axis to plot on.
    - data: pd.Series or np.array with the data used to determine y-axis limits.
    """
    ax.set_xlabel(" ")
    comparisons_default = ["M,S v. M,P", "F,S v. F,P",
                           "M,S v. F,S", "M,P v. F,P"]
    positions = [[0, 1], [2, 3], [0, 2], [1, 3]]

    third_groups = {"stress": ["Ctrl", "Stress"],
                    "zone": ["empty", "social"],
                    "day": ["1", "3"],
                    "day_scpp": ["1", "5"],
                    "object": ["1", "2"],
                    "timepoint": ["pre", "post"]}

    if third is not None:
        comparisons = []
        third_options = third_groups[third]
        for treatment in ["S", "P"]:
            comparisons.append(f"{third_options[0]},{treatment} v. {third_options[1]},{treatment}")
        for option in third_options:
            comparisons.append(f"{option},S v. {option},P")
        third_text = third.lower()
    else:
        comparisons = comparisons_default

    if third == "stress":
        comparisons = [comparisons[2], comparisons[3], comparisons[0], comparisons[1]]

    stars_sig = stars.loc[stars.isin(['*', '**', '***'])]
    try:
        stars_sig = stars_sig.drop('Intercept')
    except KeyError:
        stars_sig = stars_sig

    posthocs_present = any("v." in s for s in stars_sig.index.tolist())

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

    comparisons_plotted = [False]*len(comparisons)
    if posthocs_present:
        y_coord = y_star_locs[0]  # default y_coord.
        star_shift = {"*": 0.13, "**": 0.065, "***": 0}
        for i, comparison in enumerate(comparisons):
            if i == 2 and comparisons_plotted[0] is True and comparisons_plotted[1] is True:
                y_coord = y_star_locs[1]
            if i == 3 and comparisons_plotted[1] is True and comparisons_plotted[2] is False:
                y_coord = y_star_locs[1]
            if i == 3 and comparisons_plotted[1] is True and comparisons_plotted[2] is True:
                y_coord = y_star_locs[2]
            if i == 3 and comparisons_plotted[1] is False and comparisons_plotted[2] is True:
                y_coord = y_star_locs[1]

            matching_keys = [
                key for key in stars_sig.keys() if comparison in key]
            if len(matching_keys) > 0:
                comparisons_plotted[i] = True
                x_coords_bar = [sorted_bar_positions[positions[i][0]] + 0.05,
                                sorted_bar_positions[positions[i][1]] - 0.05]
                x_coords_star = [sorted_star_positions[positions[i][0]],
                                 sorted_star_positions[positions[i][1]]]
                y_coords = [y_coord, y_coord]
                # Plot the line over the specified bars
                ax.plot(x_coords_bar, y_coords,
                        color='black', lw=0.5, markersize=0)

                # Plot the star with appropriate text.
                star_text = ""
                for key in matching_keys:
                    if key[-1] == "]":
                        star_text = (
                            star_text + f"{(key[key.find('[')+1:-1]).lower()}{stars_sig[key]}")
                        y_loc = y_coords[0]*1.05
                        x_loc = np.mean(x_coords_star)
                    else:
                        star_text = star_text + stars_sig[key]
                        y_loc = y_coords[0]*0.95
                        x_loc = np.mean(x_coords_star) + \
                            star_shift[stars_sig[key]]

                ax.text(
                    x=x_loc,
                    y=y_loc,
                    s=star_text,
                    fontsize=8,
                )

    # If there were no posthocs, then plot main effects.
    # no comparisons were plotted. Then write main effects.
    if np.sum(comparisons_plotted) == 0:
        if third is not None:
            if "*" in stars[f"treatment:sex:{third}"]:
                plot_text = f"drug:sex:{third_text}" + \
                    stars[f"treatment:sex:{third}"]
                ax.text(x=0.5, y=0.8, s=plot_text, ha='center',
                        va='top', transform=ax.transAxes)
            elif ("*" in stars["treatment:sex"] or
                  "*" in stars[f"treatment:{third}"] or
                  "*" in stars[f"sex:{third}"]):
                if "*" in stars["treatment:sex"]:
                    plot_text = "drug:sex" + stars["treatment:sex"]
                    ax.text(x=0.5, y=0.8, s=plot_text,
                            ha='center', va='top', transform=ax.transAxes)
                    if "*" in stars[f"{third}"]:
                        if third == "stress":
                            ax.set_xlabel(stars[f"{third}"])
                        elif third == "zone":
                            ax.text(x=0.5, y=0, s=stars[f"{third}"], ha='center', va='bottom',
                                    transform=ax.transAxes)
                if "*" in stars[f"treatment:{third}"]:
                    plot_text = f"drug:{third_text}" + \
                        stars[f"treatment:{third}"]
                    ax.text(x=0.5, y=0.8, s=plot_text,
                            ha='center', va='top', transform=ax.transAxes)
                if "*" in stars[f"sex:{third}"]:
                    plot_text = f"sex:{third_text}" + stars[f"sex:{third}"]
                    ax.text(x=0.5, y=0.8, s=plot_text,
                            ha='center', va='top', transform=ax.transAxes)
                    if "*" in stars["treatment"]:
                        ax.text(x=0.6, y=0.88,
                                s=stars["treatment"], transform=ax.transAxes)
            elif "*" in stars["treatment"] or "*" in stars["sex"] or "*" in stars[f"{third}"]:
                if "*" in stars["treatment"]:
                    plot_text = stars["treatment"]
                    ax.text(x=0.6, y=0.88,
                            s=stars["treatment"], transform=ax.transAxes)
                if "*" in stars["sex"]:
                    plot_text = "sex" + stars["sex"]
                    ax.text(x=0.5, y=0.8, s=plot_text,
                            ha='center', va='top', transform=ax.transAxes)
                if "*" in stars[f"{third}"]:
                    if third == "stress":
                        ax.set_xlabel(stars[f"{third}"])
                    elif third == "zone":
                        ax.text(x=0.5, y=0, s=stars[f"{third}"], ha='center', va='bottom',
                                transform=ax.transAxes)
                    else:
                        ax.text(x=0.5, y=0, s=stars[f"{third}"], ha='center', va='bottom',
                                transform=ax.transAxes)

        else:
            if "*" in stars["treatment:sex"]:
                plot_text = "drug:sex" + stars["treatment:sex"]
                ax.text(x=0.5, y=0.8, s=plot_text, ha='center',
                        va='top', transform=ax.transAxes)

            elif "*" in stars["treatment"] or "*" in stars["sex"]:
                if "*" in stars["treatment"]:
                    plot_text = stars["treatment"]
                    ax.text(x=0.6, y=0.88,
                            s=stars["treatment"], transform=ax.transAxes)
                if "*" in stars["sex"]:
                    ax.set_xlabel(stars["sex"])


# %%
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
    data_expt_lab = psy_stats.organize_categories(data_expt_lab, plot=True)
    ymin, ymax, ybin = y_scale_auto(
        ymin, ymax, ybin, data_expt_lab[variable], variable_name)
    sns.barplot(data=data_expt_lab,
                x="sex", order=sex_order,
                y=variable,
                hue="treatment", hue_order=treatment_order,
                errorbar='se', capsize=0.2, err_kws={'linewidth': 0.5},
                palette=treatment_colors_bars, alpha=0.5,
                legend=False, ax=ax)

    sns.stripplot(data=data_expt_lab,
                  x="sex", order=sex_order,
                  y=variable,
                  hue="treatment", hue_order=treatment_order,
                  palette=treatment_colors_points,
                  size=2, linewidth=0.2, alpha=0.75, marker=markers[lab],
                  dodge=True, jitter=0.2,
                  legend=False, ax=ax)

    ax.set_ylim([ymin, ymax])
    ax.set_yticks(np.arange(ymin, ymax+0.1, ybin))
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


def plot_over_time(ax, data_of_interest, col, col_name, xticks, x_name, xlabel,
                   lab=None, ymin=None, ymax=None, ybin=None):
    """
    Plots line plots over time for males and females.
    Parameters:
    - ax: matplotlib axis to plot on.
    - data_binned_dict: dict of pd.DataFrames with binned data.
    - col: str, column name of the variable to plot.
    - col_name: str, name of the variable for labeling.
    - time_range: tuple, (start_time, end_time) for x-axis labeling.
    """
    data_of_interest = psy_stats.organize_categories(data_of_interest, plot=True)
    ymin, ymax, ybin = y_scale_auto(
        ymin, ymax, ybin, data_of_interest[col], col_name)
    data_males = data_of_interest[data_of_interest["sex"] == "M"]
    data_females = data_of_interest[data_of_interest["sex"] == "F"]
    if lab is None:
        markersize = 0
        markertype = 'o'
    else:
        markersize = 2
        markertype = markers[lab]
    sns.lineplot(data=data_males, x=x_name, y=col,
                 hue="treatment", hue_order=treatment_order,
                 palette=treatment_colors_bars, ax=ax, lw=0.5,
                 # err_kws={'lw': 0.5, 'alpha': 0.2},
                 errorbar='se', err_style='bars',
                 marker=markertype, markersize=markersize, markeredgewidth=0,
                 err_kws={'linewidth': 0.5, 'capsize': 1, 'capthick': 0.5},
                 legend=False)
    sns.lineplot(data=data_females, x=x_name, y=col,
                 hue="treatment", hue_order=treatment_order,
                 palette=treatment_colors_bars, ax=ax, ls="--", lw=0.5,
                 # err_kws={'lw': 0.5, 'alpha': 0.2},
                 errorbar='se', err_style='bars',
                 marker=markertype, markersize=markersize, markeredgewidth=0,
                 err_kws={'linewidth': 0.5, 'capsize': 1, 'capthick': 0.5},
                 legend=False)

    ax.set_ylim([ymin, ymax])
    ax.set_yticks(np.arange(ymin, ymax+0.1, ybin))
    ax.set_ylabel(col_name)

    ax.set_xlim([-1, xticks[-1] + 1])
    ax.set_xlabel(xlabel)
    ax.set_xticks(xticks)

    if lab is None:
        ax.set_title("All labs")
        # add title and legend
        title_rect = patches.Rectangle((0, 1.05),
                                       1, 0.17, facecolor='white', edgecolor='black', linewidth=0.5,
                                       zorder=0, transform=ax.transAxes, clip_on=False)
        ax.add_patch(title_rect)
    else:
        ax.set_title(lab, color=custom_palette[lab])

    sns.despine()
    plt.tight_layout()


def plot_all_labs(ax, data_expt, variable, variable_name, stars=None,
                  ymin=None, ymax=None, ybin=None):
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
    data_expt = psy_stats.organize_categories(data_expt, plot=True)
    ymin, ymax, ybin = y_scale_auto(
        ymin, ymax, ybin, data_expt[variable], variable_name)
    sns.barplot(data=data_expt,
                x="sex", order=sex_order,
                y=variable,
                hue="treatment", hue_order=treatment_order,
                errorbar='se', capsize=0.2, err_kws={'linewidth': 0.5},
                palette=treatment_colors_bars, alpha=0.5,
                legend=False, ax=ax)
    for lab_idx, lab in enumerate(labs):
        sns.stripplot(data=data_expt[data_expt["institution"] == lab],
                      x="sex", order=sex_order,
                      y=variable,
                      hue="treatment", hue_order=treatment_order,
                      palette=treatment_colors_points,
                      size=1.35, linewidth=0, alpha=0.75, jitter=0.2, dodge=True,
                      marker=markers[lab],
                      legend=False, ax=ax)

    ax.set_ylim([ymin, ymax])
    ax.set_yticks(np.arange(ymin, ymax+0.1, ybin))
    ax.set_xlabel("")
    ax.set_xticks(sex_order, ["M", "F"])
    ax.set_ylabel(variable_name)

    # plot stars
    if stars is not None:
        plot_stars(stars, ax, data_expt[variable])

    # add title and legend
    ax.set_title("All labs")
    title_rect = patches.Rectangle((0, 1.05),
                                   1, 0.17, facecolor='white', edgecolor='black', linewidth=0.5,
                                   zorder=0, transform=ax.transAxes, clip_on=False)

    ax.add_patch(title_rect)
    add_legend(ax)
    sns.despine()
    plt.tight_layout()


def plot_zscores(ax, psi_data, variable, ymin=-2, ymax=2, ybin=2, x="sex", x_order=None):
    """
    Plots standardized psi effect data.
    Parameters:
    - psi_data: pd.DataFrame with standardized psi data.
    - variable: str, column name of the standardized variable to plot.
    - ymin: minimum y-value for the plot.
    - ymax: maximum y-value for the plot.
    - ybin: y-axis bin size.
    """
    if x_order is not None:
        order = x_order
    else:
        order = None
    sns.pointplot(
        data=psi_data, x=x, y=variable, order=order,
        hue="institution", hue_order=labs,
        errorbar="ci", capsize=0.1, err_kws={'linewidth': 0.5},
        markersize=2, linewidth=0.5,
        palette=custom_palette, linestyle="none", dodge=0.5,
        legend=False
    )
    ax.set_ylabel("z-score")
    ax.set_xlabel("")
    ax.set_ylim([ymin, ymax])
    ax.set_yticks(np.arange(ymin, ymax+1, ybin))
    ax.axhline(0, lw=0.5, ls='--', color='black')
    sns.despine()
    ax.set_title("Psi Effect")
    plt.tight_layout()


def plot_bars_thirdfactor(ax, data_expt_lab, variable, variable_name, lab,
                          stars=None, ymin=None, ymax=None, ybin=None, data_expt=None):
    data_expt_lab = psy_stats.organize_categories(data_expt_lab)
    ymin, ymax, ybin = y_scale_auto(
        ymin, ymax, ybin, data_expt_lab[variable], variable_name)
    if data_expt is not None:
        data_bars = data_expt
        alpha_bars = 0.5/5
    else:
        data_bars = data_expt_lab
        alpha_bars = 0.5
    sns.barplot(data=data_bars,
                x="stress", order=stress_order,
                y=variable,
                hue="treatment", hue_order=treatment_order,
                errorbar='se', capsize=0.2, err_kws={'linewidth': 0.5},
                palette=treatment_colors_bars, alpha=alpha_bars,
                legend=False, ax=ax)

    sns.stripplot(data=data_expt_lab[data_expt_lab["sex"] == 'M'],
                  x="stress", order=stress_order,
                  y=variable,
                  hue="treatment", hue_order=treatment_order,
                  palette=treatment_colors_points,
                  size=2, linewidth=0.2, alpha=0.75, marker=markers[lab],
                  dodge=True, jitter=0.2,
                  legend=False, ax=ax)
    sns.stripplot(data=data_expt_lab[data_expt_lab["sex"] == 'F'],
                  x="stress", order=stress_order,
                  y=variable,
                  hue="treatment", hue_order=treatment_order,
                  palette=treatment_colors_points,
                  size=2, linewidth=0.2, alpha=0.5, marker=markers[lab],
                  dodge=True, jitter=0.2,
                  legend=False, ax=ax)
    collections = ax.collections[-4:]
    colors = ["gray", "green", "gray", "green"]
    for i, artist in enumerate(collections):
        artist.set_facecolors("white")
        artist.set_edgecolors(colors[i])
    ax.set_ylim([ymin, ymax])
    ax.set_yticks(np.arange(ymin, ymax+0.1, ybin))
    ax.set_xlabel("")
    ax.set_xticks(stress_order, ["Ctrl", "Stress"])
    ax.set_title(lab, color=custom_palette[lab])
    ax.set_ylabel(variable_name)

    # then do stats
    if stars is not None:
        plot_stars(stars, ax, data_expt_lab[variable], third="stress")
    add_legend(ax, alpha_bars)

    sns.despine()
    plt.tight_layout()


def plot_pairs(ax, data, variable, variable_name, lab, other_var,
               stars=None, ymin=None, ymax=None, ybin=None, ylabel=False):

    ymin, ymax, ybin = y_scale_auto(
        ymin, ymax, ybin, data[variable], variable_name)

    other_var_values = np.unique(data[other_var])
    value1 = other_var_values[0]
    value2 = other_var_values[1]
    v1 = value1[0]
    v2 = value2[0]

    # Positioning for mini-bars
    x_positions = {"S": 0, "P": 0.5}
    color_map_bars = {
        "S": treatment_colors_bars[0], "P": treatment_colors_bars[1]}
    color_map_lines = {
        "S": treatment_colors_points[0], "P": treatment_colors_points[1]}
    ls_map = {"M": '-', "F": '--'}
    offset = 0.18   # separation between empty/social mini-bars
    # ---- BAR MEANS ----
    for t in treatment_order:
        other_var_order = sorted([value1, value2])
        for i, z in enumerate(other_var_order):
            df_sub = data[(data["treatment"] == t) & (data[other_var] == z)]
            mean = df_sub[variable].mean()
            sem = df_sub[variable].sem()

            xpos = x_positions[t] + (i - 0.5) * offset
            # Bar with mean
            ax.bar(xpos, mean, yerr=sem, color=color_map_bars[t],
                   alpha=0.65, width=0.16, linewidth=0.5,
                   error_kw=dict(lw=0.5, capsize=1, capthick=0.5))
            # Add label under bar
            label = v1 if z == value1 else v2
            ax.text(xpos, 0,      # position inside the bar (10% up from bottom)
                    label, ha="center", va="bottom", fontsize=8)

    # ---- INDIVIDUAL CONNECTING LINES ----
    for mouse in data["mouse_ID"].unique():
        df_mouse = data[data["mouse_ID"] == mouse]
        sex = df_mouse["sex"].iloc[0]
        for t in treatment_order:
            df_pair = df_mouse[df_mouse["treatment"] == t]
            if df_pair.shape[0] == 2:   # has both other_var_values
                x_value1 = x_positions[t] - offset/2
                x_value2 = x_positions[t] + offset/2

                y_value1 = df_pair[df_pair[other_var]
                                   == value1][variable].values[0]
                y_value2 = df_pair[df_pair[other_var]
                                   == value2][variable].values[0]
                ax.plot(
                    [x_value1, x_value2],
                    [y_value1, y_value2],
                    linewidth=0.5, ls=ls_map[sex], color=color_map_lines[t], alpha=0.5
                )
    ax.set_xticks([x_positions["S"], x_positions["P"]], ["S", "P"])
    ax.set_xlabel("")
    ax.set_ylim([ymin, ymax])
    ax.set_yticks(np.arange(ymin, ymax+1, ybin))
    if ylabel:
        ax.set_ylabel(variable_name)

    if lab == "All labs":
        # add title and legend
        ax.set_title(lab)
        title_rect = patches.Rectangle((0, 1.05),
                                       1, 0.17, facecolor='white', edgecolor='black', linewidth=0.5,
                                       zorder=0, transform=ax.transAxes, clip_on=False)
        ax.add_patch(title_rect)

    else:
        ax.set_title(lab, color=custom_palette[lab])

    # then do stats
    if stars is not None:
        plot_stars(stars, ax, data[variable], third=other_var)
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
        ybin = ((ymax-ymin)/2)
    return ymin, ymax, ybin
