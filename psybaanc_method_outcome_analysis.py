# -*- coding: utf-8 -*-
"""
Created on Mon Feb 16 09:49:18 2026

@author: olu
"""
# %% Import packages
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics.pairwise import cosine_similarity
from scipy.stats import pearsonr

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
data_path = r"C:\Users\olu\Documents\Psy-BAANC\Paper Drafts\revision materials\Supplementary Table 1_031726.xlsx"
save_path = r"Y:/PsyBAANC/figures/final/cossim"
labs = ["Stanford", "Berkeley 1", "Berkeley 2", "UCSF 1", "UCSF 2"]

id_cols = ["mouse_ID", "cage_ID", "institution", "lab", "experiment", "treatment", "sex", "stress"]
cols_exclude_all = ["prior_expts_notes",
                    ]
cols_methods_zscore_all = [
    "timepoint",
    "weight", "age",
    "handling_acc", "injection_acc", "lab_acc",
    "h_wall", "L", "w", "h_floor", "h", "epm_lip", "diameter", "water_depth", "water_temp"
    "distance_objects_wall", "object_base_area", "object_height",
    "cagemates",
    "other_mice", "prior_expts",
    "same_treatment", "expt_time",
    "tone_frequency", "tone_amplitude", "expt_time_extinction1"
]  #
cols_methods_custom_all = ["exptr_sex",
                           "lighting", "lighting_open_arms", "lighting_closed_arms",
                           "experiment_room_colony_distance", "transport", "acclimation_room",
                           "enrichment", "mouse_source",
                           "injection_context",
                           "analysis_coordinates", "analysis_variables",
                           "experiment_location", "counterbalanced_location",
                           "pip_pure",
                           "object_type"
                           ]

experiment = "CORT (FST)"
experiment_name = "cortFST"
variable = "time_immobile"


# %% Functions
def mantel_test(dist1, dist2, n_perm=10000):
    """
    Pearson correlation between upper triangles
    with permutation test.
    """

    # Extract upper triangle
    triu_idx = np.triu_indices_from(dist1, k=1)
    # Stanford = triu_idx[0] == 0
    # triu_idx = (triu_idx[0][Stanford], triu_idx[1][Stanford])
    x = dist1[triu_idx]
    y = dist2[triu_idx]

    observed_r, _ = pearsonr(x, y)

    permuted_rs = []

    for _ in range(n_perm):
        perm = np.random.permutation(dist1.shape[0])
        permuted = dist2[perm][:, perm]
        y_perm = permuted[triu_idx]
        r_perm, _ = pearsonr(x, y_perm)
        permuted_rs.append(r_perm)

    permuted_rs = np.array(permuted_rs)

    p_value = np.mean(np.abs(permuted_rs) >= np.abs(observed_r))

    return observed_r, p_value, x, y


# %% Main
# Read in raw data.
data_all = pd.read_excel(data_path, sheet_name=None)

# %% Get dataframes for methods and outcomes
data = data_all[experiment]

data_expt = data[data["experiment"] == experiment_name]
cols_exclude = list(set(data_expt.columns) & set(cols_exclude_all))
data_expt_clean = data_expt.drop(columns=cols_exclude).dropna()
id_cols_expt = list(set(data_expt.columns) & set(id_cols))
data_id = data_expt_clean[id_cols_expt]

cols_methods_custom = list(set(data_expt_clean.columns) & set(cols_methods_custom_all))
data_methods_categ = data_expt_clean[cols_methods_custom]
data_methods_categ = pd.get_dummies(data_methods_categ, dtype='int')

# Get methods dataframe
cols_methods_zscore = list(set(data_expt_clean.columns) & set(cols_methods_zscore_all))
data_methods_zscore = data_expt_clean[cols_methods_zscore]
data_methods_zscore = pd.concat((data_methods_zscore, data_methods_categ), axis=1)
data_methods_zscore = (data_methods_zscore - data_methods_zscore.mean()) / \
    data_methods_zscore.std(ddof=1)
for col in data_methods_zscore.columns:
    if data_methods_zscore[col].isna().all():
        data_methods_zscore[col] = 0

data_methods = pd.concat([data_id, data_methods_zscore], axis=1)

# Get outcomes dataframe
data_outcomes = data_expt[id_cols_expt + [variable]]
data_outcomes = data_outcomes.dropna()

# %% Get vectors and cosine similarity scores for methods and outcomes
# Methods
method_vectors = np.empty(
    (len(labs), (len(data_methods_zscore.columns))*2))
for lab_i, lab in enumerate(labs):
    psi_methods = []
    for sex in ["M", "F"]:
        psi = data_methods[(data_methods["institution"] == lab) & (data_methods["sex"] == sex)
                           & (data_methods["treatment"] == "P")]
        psi = psi.drop(columns=id_cols_expt)
        psi_methods.extend(psi.mean().values)
    method_vectors[lab_i] = psi_methods

mask = ~np.isnan(method_vectors).any(axis=0)
method_vectors = method_vectors[:, mask]

# Outcomes
if "stress" in data_outcomes.columns:
    outcome_vectors = np.empty((len(labs), 6))
else:
    outcome_vectors = np.empty((len(labs), 2))
for lab_i, lab in enumerate(labs):
    psi_zscore = []
    for sex in ["M", "F"]:
        sal = data_outcomes[(data_outcomes["institution"] == lab) & (data_outcomes["sex"] == sex)
                            & (data_outcomes["treatment"] == "S")]
        if "stress" in data_outcomes.columns:
            sal_stress = sal[sal["stress"] == "Stress"]
            sal = sal[sal["stress"] == "Ctrl"]
        sal_mean = sal[variable].mean()
        sal_std = sal[variable].std(ddof=1)
        psi = data_outcomes[(data_outcomes["institution"] == lab) & (data_outcomes["sex"] == sex)
                            & (data_outcomes["treatment"] == "P")]
        if "stress" in data_outcomes.columns:
            for stress_group in ["Ctrl", "Stress"]:
                psi_new = psi[psi["stress"] == stress_group]
                psi_zscore.append(((psi_new[variable]-sal_mean) / sal_std).mean())
            psi_zscore.append(((sal_stress[variable]-sal_mean) / sal_std).mean())
        else:
            psi_zscore.append(((psi[variable]-sal_mean) / sal_std).mean())
    outcome_vectors[lab_i] = psi_zscore


cossim_methods = cosine_similarity(method_vectors)
cossim_outcomes = cosine_similarity(outcome_vectors)
method_distance = 1 - cossim_methods
behavior_distance = 1 - cossim_outcomes


# %% Mantel test
mantel_r, mantel_p, x, y = mantel_test(
    method_distance,
    behavior_distance,
    n_perm=10000
)

print("Mantel correlation:", mantel_r)
print("Permutation p-value:", mantel_p)


# %% Plot
# Set ticks
vmin_methods = round(float(max(np.min(x) - 0.1, 0)), 1)
vmax_methods = round(float(min(np.max(x) + 0.1, 2)), 1)
vmin_behavior = round(float(max(np.min(y) - 0.1, 0)), 1)
vmax_behavior = round(float(min(np.max(y) + 0.1, 2)), 1)

fig, ax = plt.subplots(nrows=3, ncols=1, figsize=[1.5, 4])
im = ax[0].imshow(method_distance, cmap='Reds', vmin=vmin_methods, vmax=vmax_methods)
ax[0].set_title("Methods")
ax[0].set_yticks(np.arange(method_distance.shape[0]))
ax[0].set_yticklabels(["S", "B1", "B2", "U1", "U2"])
ax[0].set_xticks(np.arange(method_distance.shape[0]))
ax[0].set_xticklabels(["S", "B1", "B2", "U1", "U2"])

im2 = ax[1].imshow(behavior_distance, cmap='Reds', aspect='auto',
                   vmin=vmin_behavior, vmax=vmax_behavior)
ax[1].set_title("Behavior")
ax[1].set_xticks(np.arange(method_distance.shape[0]))
ax[1].set_xticklabels(["S", "B1", "B2", "U1", "U2"])
ax[1].set_yticks(np.arange(method_distance.shape[0]))
ax[1].set_yticklabels(["S", "B1", "B2", "U1", "U2"])

ax[2].scatter(x, y, s=5,
              color='black'
              )
ax[2].set_box_aspect(1)  # Makes the scatter plot square, matching the image
ax[2].set_xticks([vmin_methods, vmax_methods])
ax[2].set_yticks([vmin_behavior, vmax_behavior])
ax[2].set_xlabel("Methods")
ax[2].set_ylabel("Behavior")
ax[2].set_title("Correlation")

cb = fig.colorbar(im, ax=ax[0], location='right', pad=0.05)
cb.set_ticks([vmin_methods, vmax_methods])

cb = fig.colorbar(im2, ax=ax[1], location='right', pad=0.05)
cb.set_ticks([vmin_behavior, vmax_behavior])

plt.tight_layout()
plt.savefig(os.path.join(save_path, f"{experiment_name}_{variable}_cossim.svg"), transparent=True)
