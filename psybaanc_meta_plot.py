# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 14:59:43 2026

@author: olu
"""

# %% Import packages
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

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
data_path = r"C:\Users\olu\Documents\Psy-BAANC\Data\Gomez\Gomez Final\loose_model_statistics_FINAL.xlsx"
save_path = r"Y:/PsyBAANC/figures/final/meta"


data_og = pd.read_excel(data_path, sheet_name="fixed_effects_clean")
data_interact = pd.read_excel(data_path, sheet_name="interactions_alpha_0.05")

replacement_key = {"expt_time_binlate": "time_late", "expt_time_binmiddle": "time_middle",
                   "enrichmentnestlet_hut": "enrich_max", "enrichmentnestlet_paper": "enrich_mid2",
                   "enrichmentnone": "enrich_min",
                   "mouse_sourceJax": "source_Jax",
                   "injection_contexthome_cage": "inject_home",
                   "lightingdimly_lit": "light_dim", "lightingfully_lit": "light_bright",
                   "exptr_sexM": "exptr_M", "exptr_sexM + F": "exptr_M+F",
                   "stressStress": "stress",
                   "same_treatment1": "same_drug",
                   "analysis_coordinatesDeepLabCut": "DeepLabCut", "analysis_coordinatesethovision": "Ethovision",
                   "experiment_room_colony_distancenearby_room": "nearby_room",
                   "object_type250 mL glass bottle": "250_bottle", "object_type50 mL Erlenmeyer Flask": "50_flask",
                   "object_typeLego tower (Lego 10698)": "lego_tower", "object_typelego": "lego",
                   "object_typeTissue flask (Merck CLS353018)": "tissue_flask",
                   "social_cup_baseArea": "cup_area",
                   "treatmentP": "drugP",
                   "cagemates": "cage_size",
                   "tone_frequency": "tone_freq",
                   "handling_acc": "handling",


                   "treatmentP:exptr_sexM": "drugP:exptr_M",
                   "treatmentP:exptr_sexM + F": "drugP:exptr_M+F",
                   "treatmentP:enrichmentnestlet_paper": "drugP:enrich_mid2",
                   "treatmentP:enrichmentnestlet_hut": "drugP:enrich_max",
                   "treatmentP:enrichmentnone": "drugP:enrich_min",
                   "treatmentP:stressStress": "drugP:stress",
                   "treatmentP:injection_contexthome_cage": "drugP:inject_home",
                   "treatmentP:lightingfully_lit": "drugP:light_bright",
                   "treatmentP:mouse_sourceJax": "drugP:source_Jax",
                   "treatmentP:other_mice": "drugP:other_mice",
                   "treatmentP:age": "drugP:age",
                   "treatmentP:sexF": "drugP:sexF",
                   "treatmentP:cagemates": "drugP:cage_size",
                   "treatmentP:same_treatment": "drugP:same_drug",
                   "treatmentP:lab_acc": "drugP:lab_acc",
                   "treatmentP:handling_acc": "drugP:handling"
                   }

# %% Functions


def sort_df(df):
    df.replace(replacement_key, inplace=True)
    condition = (df['variable'] == 'drugP') | (df['variable'] == 'sexF')
    condition2 = df['variable'].str.contains(':', na=False)
    first_rows = df[condition].sort_values(by='variable', ascending=True)
    middle_rows = df[~condition & ~condition2]
    last_rows = df[condition2]
    sorted = middle_rows.sort_values(by="variable", ascending=True)
    df_new = pd.concat([first_rows, sorted, last_rows])
    return df_new


# %% All interactions
experiments_interact = np.unique(data_interact["experiment"])
experiments_og = np.unique(data_og["experiment"])

for experiment in experiments_og:

    data_expt_og = data_og[data_og["experiment"] == experiment].sort_values(by='variable')
    data_expt_og["interactions"] = "no"
    data_expt_og = sort_df(data_expt_og)
    x_positions_og = np.arange(len(data_expt_og))

    if experiment in experiments_interact:
        data_expt_interact = data_interact[data_interact["experiment"] == experiment]
        data_expt_interact["interactions"] = "yes"
        data_expt_interact = sort_df(data_expt_interact)
        x_positions_interact = np.arange(len(data_expt_interact))

        data_expt = pd.concat((data_expt_interact, data_expt_og))
    else:
        data_expt = data_expt_og

    data_expt.replace(replacement_key, inplace=True)
    x_positions = np.arange(len(data_expt))

    # 2. Create the base point plot with error bars disabled
    plt.figure(figsize=(0.7, 2.5))
    ax = sns.pointplot(
        data=data_expt, x="estimate", y="variable",
        hue="interactions", hue_order=['yes', 'no'], palette=['red', 'darkgray'],
        join=False, legend=False,
        markersize=3, marker='D', markerfacecolor='white', markeredgewidth=0.5, alpha=0.75
    )

    if experiment in experiments_interact:
        y_err_lower = data_expt_interact['estimate'] - data_expt_interact['ci_low_95']
        y_err_upper = data_expt_interact['ci_high_95'] - data_expt_interact['estimate']
        ax.errorbar(x=data_expt_interact['estimate'], y=x_positions_interact, xerr=[y_err_lower, y_err_upper],
                    fmt='none', c='red', capsize=2, capthick=0.5, zorder=1, elinewidth=0.5, alpha=0.5)
    y_err_lower = data_expt_og['estimate'] - data_expt_og['ci_low_95']
    y_err_upper = data_expt_og['ci_high_95'] - data_expt_og['estimate']
    ax.errorbar(x=data_expt_og['estimate'], y=x_positions_og, xerr=[y_err_lower, y_err_upper],
                fmt='none', c='darkgray', capsize=2, capthick=0.5, zorder=1, elinewidth=0.5,
                alpha=0.75)

    ax.axvline(x=0, lw=0.25, ls='--', color='gray')
    ax.set_title(experiment)
    ax.set_ylabel("")
    sns.despine()
    plt.savefig(os.path.join(save_path, f"{experiment}_ModelSummary.svg"), transparent=True)
    plt.show()


# %% does adding interactions increase variance explained by fixed effects?

dfs = []
for data in [data_og, data_interact]:
    experiments = []
    r2_marginal = []
    r2_conditional = []
    for experiment in experiments_interact:
        experiments.append(experiment)
        r2_marginal.append(data_interact["R2_marginal"]
                           [data_interact["experiment"] == experiment].iloc[0])
        r2_conditional.append(data_interact["R2_conditional"]
                              [data_interact["experiment"] == experiment].iloc[0])

    dfs.append(pd.DataFrame(data={"experiment": experiments,
                            "r2_marginal": r2_marginal, "r2_conditional": r2_conditional}))

df_final = pd.concat(dfs)
df_final["interactions"] = ["no"]*len(experiments) + ["yes"]*len(experiments)
df_long = pd.melt(df_final, id_vars=['experiment', 'interactions'], value_vars=['r2_conditional', 'r2_marginal'],
                  var_name='type', value_name='r2')

fig, ax = plt.subplots(figsize=(0.7, 1))
sns.stripplot(ax=ax, data=df_final, x="interactions", y="r2_marginal", hue='experiment', palette=['black'],
              jitter=False, legend=False, size=2)
sns.lineplot(ax=ax, data=df_final, x='interactions', y='r2_marginal', hue='experiment', palette=['black'],
             estimator=None, color='.7', alpha=0.5, lw=0.5, legend=False)
ax.set_ylim([0, 1])
sns.despine()
plt.savefig(os.path.join(save_path, "variance_w_interactions.svg"), transparent=True)

# %% how much of variance between labs can be explained by identifiable factors?
# R2 conditional - R2 marginal = residual between lab variance still unexplained after modeling fixed effects
# ICC = baseline between lab variance

df_lab_variance = pd.DataFrame(columns=["experiment", "methods included", "lab variance"])
methods_lab_variance = []
null_lab_variance = []
for experiment in experiments_og:
    null_lab_variance.append(data_og["null model icc"][data_og["experiment"] == experiment].iloc[0])
    methods_lab_variance.append(data_og["R2_conditional - R2_marginal"]
                                [data_og["experiment"] == experiment].iloc[0])
    print(f"{experiment}: {round(methods_lab_variance[-1] - null_lab_variance[-1], 3)}")

df_lab_variance["lab variance"] = null_lab_variance + methods_lab_variance
df_lab_variance["methods included"] = ["no"]*len(experiments_og) + ["yes"]*len(experiments_og)
df_lab_variance["experiment"] = list(experiments_og) + list(experiments_og)

fig, ax = plt.subplots(figsize=(0.7, 1))
sns.stripplot(ax=ax, data=df_lab_variance, x="methods included", y="lab variance", hue='experiment', palette=['black'],
              jitter=False, legend=False, size=2)
sns.lineplot(ax=ax, data=df_lab_variance, x='methods included', y='lab variance', hue='experiment', palette=['black'],
             estimator=None, color='.7', alpha=0.5, lw=0.5, legend=False)

ax.set_ylim([0, 1])
# ax.axhline(0.2, c='red', ls='--', lw=0.5, alpha=0.5)
sns.despine()
plt.savefig(os.path.join(save_path, "variance_lab_w_methods.svg"), transparent=True)
