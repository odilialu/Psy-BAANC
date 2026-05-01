# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 10:36:16 2026

@author: olu
"""

# %% Import packages
import os
import matplotlib.colors as mcolors
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random
import seaborn as sns
import functions.psybaanc_stats as psy_stats
from scipy import stats
import statsmodels.api as sm
import statsmodels.formula.api as smf

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
data_path = r"C:/Users/olu/Documents/Psy-BAANC/Paper Drafts/revision materials\Supplementary Table 1_final_250210.xlsx"
save_path = r"Y:/PsyBAANC/figures/final/simulations"
labs = ["Stanford", "Berkeley 1", "Berkeley 2", "UCSF 1", "UCSF 2"]
custom_palette = dict(
    zip(labs, ["darkorange", "orangered", "maroon", "purple", "slateblue"]))
n_sims = 10000
n_samples = [5, 10, 20, 40]

orig_cmap = plt.get_cmap('binary')
new_colors = orig_cmap(np.linspace(0.25, 1.0, 128))
new_cmap = mcolors.ListedColormap(new_colors)

experiments_dict = {
    "HTR": ["HTR"],
    "OFT": ["aOFT", "pOFT"],
    "EPM": ["aEPM", "pEPM"],
    "NOE": ["aNOE", "pNOE"],
    "SIT": ["aSIT", "pSIT"],
    "FST": ["rFST"],
    "TST": ["rTST"],
    "FC": ["pre_retrieval", "post_retrieval"],
    "sCPP": ["juv", "adult"],
    "CORT (EPM)": ["cortEPM"],
    "CORT (SPT)": ["cortSPT"],
    "CORT (FST)": ["cortFST"]
}

vars_dict = {
    "HTR": ["htr_total"],
    "OFT": ["time_center"],
    "EPM": ["time_open"],
    "NOE": ["Average_time"],
    "SIT": ["social_index"],
    "FST": ["day_3_immobility"],
    "TST": ["day_3_immobility"],
    "FC": ["retrieval_avg"],
    "sCPP": ["social_preference_score"],
    "CORT (EPM)": ["time_center"],
    "CORT (SPT)": ["pref_index"],
    "CORT (FST)": ["immobility"]
}


# %% Read in raw data.
data_all = pd.read_excel(data_path, sheet_name=None)


# %% Effect of small sample on effect size
for experiment in ["NOE"]:
    data = data_all[experiment]
    experiment_names = experiments_dict[experiment]
    for experiment_name in ["aNOE"]:
        data_expt = data[data["experiment"] == experiment_name]
        for variable in vars_dict[experiment]:
            data_expt_zscores = psy_stats.get_zscores(data_expt, variable)
            zscore_labs = []
            zscore_labs_M = []
            zscore_labs_F = []
            for lab in labs:
                zscore_labs.append(data_expt_zscores[variable]
                                   [(data_expt_zscores["institution"] == lab)].mean())
                zscore_labs_M.append(data_expt_zscores[variable]
                                     [(data_expt_zscores["institution"] == lab)
                                      & (data_expt_zscores["sex"] == "M")].mean())
                zscore_labs_F.append(data_expt_zscores[variable]
                                     [(data_expt_zscores["institution"] == lab)
                                      & (data_expt_zscores["sex"] == "F")].mean())

            sample_zscores = np.full((n_sims, len(n_samples)), np.nan)
            for i, sample_size in enumerate(n_samples):
                for j in range(n_sims):
                    sample = random.sample(data_expt_zscores[variable].tolist(), sample_size)
                    sample_zscores[j, i] = np.mean(sample)

            sample_zscores = pd.DataFrame(sample_zscores, columns=n_samples)
            sample_zscores_long = pd.melt(
                sample_zscores, value_name='z-score', var_name='sample_size')

            effect = sample_zscores_long['z-score'].mean()
            summary_df = pd.DataFrame(columns=["sample_size", "effect_type", "probability"])
            for sample_size in n_samples:
                df_temp = sample_zscores_long[sample_zscores_long["sample_size"] == sample_size]
                if effect > 0:
                    n_op = np.sum(df_temp["z-score"] < 0)/len(df_temp)*100
                    n_stronger = np.sum(df_temp["z-score"] > (1+effect))/len(df_temp)*100
                elif effect < 0:
                    n_op = np.sum(df_temp["z-score"] > 0)/len(df_temp)*100
                    n_stronger = np.sum(df_temp["z-score"] < (effect-1))/len(df_temp)*100

                summary_df = pd.concat((summary_df, pd.DataFrame([{"sample_size": sample_size,
                                                                   "effect_type": "opposite",
                                                                   "probability": n_op
                                                                   }])))
                summary_df = pd.concat((summary_df, pd.DataFrame([{"sample_size": sample_size,
                                                                   "effect_type": "z+1",
                                                                   "probability": n_stronger
                                                                   }])))

            fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(1.25, 3.5))
            sns.kdeplot(sample_zscores_long, x='z-score', hue='sample_size', palette=new_cmap,
                        lw=1, legend=False, ax=ax[0])
            ax[0].axvline(data_expt_zscores[variable].mean(), color='red', ls='--', lw=0.5)
            # for lab_i, lab in enumerate(labs):
            # ax.axvline(zscore_labs[lab_i], color=custom_palette[lab], lw=0.5)
            # ax.axvline(zscore_labs_F[lab_i], color=custom_palette[lab], lw=0.5, ls='--')

            ax[0].set_ylim([0, 1])
            ax[0].set_yticks([0, 0.5, 1])
            ax[0].set_title(experiment_name)
            # plt.tight_layout()
            # plt.savefig(os.path.join(
            #     save_path, f"{experiment_name}_{variable}_sampleXeffect_curves.svg"), transparent=True)
            # plt.show()
            # plt.close()

            # fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(1.25, 2.2), sharey=False)
            sns.barplot(data=summary_df[summary_df["effect_type"] == "opposite"],
                        x="sample_size", y="probability",
                        hue="sample_size", palette=new_cmap,
                        legend=False, ax=ax[1])
            sns.barplot(data=summary_df[summary_df["effect_type"] == "z+1"],
                        x="sample_size", y="probability",
                        hue="sample_size", palette=new_cmap,
                        legend=False, ax=ax[2])
            ax[1].set_ylabel("op. z chance %")
            ax[1].set_xlabel("Sample size")
            ax[1].set_yticks([0, 20, 40])
            ax[2].set_yticks([0, 3, 6])
            ax[2].set_xlabel("Sample size")
            ax[2].set_ylabel("1+z chance %")
            plt.tight_layout()
            plt.savefig(os.path.join(
                save_path, f"{experiment_name}_{variable}_sampleXeffect.svg"), transparent=True)
            plt.show()
            plt.close()


# %%
for experiment in ["FST"]:
    data = data_all[experiment]
    experiment_names = experiments_dict[experiment]
    for experiment_name in ["rFST"]:
        data_expt = data[data["experiment"] == experiment_name]
        for variable in vars_dict[experiment]:
            data_expt_labs = []
            for lab in labs:
                data_expt_lab = data_expt[data_expt["institution"] == lab]
                data_sal = data_expt_lab[data_expt_lab["treatment"] == "S"]
                data_sal_M = data_sal[variable][data_sal["sex"] == "M"]
                data_sal_F = data_sal[variable][data_sal["sex"] == "F"]

                data_psi = data_expt_lab[data_expt_lab["treatment"] == "P"]
                data_psi_M = data_psi[variable][data_psi["sex"] == "M"]
                data_psi_F = data_psi[variable][data_psi["sex"] == "F"]

                data_expt_zscores = data_expt_lab.copy()
                data_expt_zscores["zscores"] = np.nan

                data_expt_zscores.loc[data_expt_zscores["sex"] == "M", "zscores"] = (
                    data_expt_lab[variable] - data_sal_M.mean()) / data_sal_M.std(ddof=1)

                data_expt_zscores.loc[data_expt_zscores["sex"] == "F", "zscores"] = (
                    data_expt_lab[variable] - data_sal_F.mean()) / data_sal_F.std(ddof=1)

                data_expt_labs.append(data_expt_zscores)

            data_expt_zscores = pd.concat(data_expt_labs)

            fig, ax = plt.subplots(figsize=(1.35, 1.35))
            sns.kdeplot(data_expt_zscores[data_expt_zscores["sex"] == "M"], x="zscores",
                        hue="treatment", hue_order=["S", "P"], palette=['gray', 'green'],
                        common_norm=False, legend=False,
                        lw=0.5, ax=ax)
            sns.kdeplot(data_expt_zscores[data_expt_zscores["sex"] == "F"],  x="zscores",
                        hue="treatment", hue_order=["S", "P"], palette=['gray', 'green'],
                        common_norm=False, legend=False,
                        lw=0.5, ls='--', ax=ax)
            plt.xticks([-5, 0, 5])
            plt.yticks([0, 0.2, 0.4, 0.6])
            plt.ylim([0, 0.6])
            plt.xlim([-5, 5])
            plt.xlabel('z-score')
            plt.tight_layout()
            plt.savefig(os.path.join(
                save_path, f"{experiment_name}_{variable}_kdeDistr.svg"), transparent=True)
            plt.show()

            file_p = os.path.join(save_path, f"{experiment_name}_{variable}_simsdata.npy")
            if os.path.exists(file_p):
                p_drug = np.load(file_p)
            else:
                p_interaction = np.full((n_sims, len(n_samples)), np.nan)
                p_drug = np.full((n_sims, len(n_samples)), np.nan)
                for i, sample_size in enumerate(n_samples):
                    print(f"-------- Sample size: {sample_size} -----------")
                    for j in range(n_sims):
                        sample_dfs = []
                        for sex in ["M", "F"]:
                            for treatment in ["S", "P"]:
                                data_temp = data_expt_zscores[(data_expt_zscores["sex"] == sex) & (
                                    data_expt_zscores["treatment"] == treatment)]
                                data_kde = stats.gaussian_kde(data_temp["zscores"])
                                sample_from_kde = data_kde.resample(sample_size).ravel()
                                sample_dfs.append(pd.DataFrame(
                                    data={'y': sample_from_kde, 'drug': treatment, 'sex': sex}))

                        sample_df = pd.concat(sample_dfs)
                        formula = 'y ~ C(drug, Sum)*C(sex, Sum)'
                        ols_lm = smf.ols(formula, data=sample_df)
                        fit = ols_lm.fit()
                        results = sm.stats.anova_lm(fit, typ=3)
                        p_interaction[j, i] = results.loc["C(drug, Sum):C(sex, Sum)", "PR(>F)"]
                        p_drug[j, i] = results.loc["C(drug, Sum)", "PR(>F)"]
                        if j % 1000 == 0:
                            print(f"Sim no. {j}")
                np.save(file_p, p_drug)

        fig, ax = plt.subplots(figsize=(1.45, 1.35))
        sns.histplot(data=p_drug[:, 1], bins=20, alpha=0.3, color='r',
                     stat="probability", label='drug', legend=False)
        # sns.histplot(data=p_interaction[:, 1], bins=20, alpha=0.3, color='b',
        #              stat="probability", label='interaction', legend=False)
        ax.set_xlabel('p-value')
        # plt.legend(loc='upper right')
        ax.axvline(x=0.05, ymin=0, ymax=1, color='r', linestyle='--', lw=0.5)
        ax.set_xlim([0, 1])
        ax.set_ylim([0, 1])
        ax.set_xticks((0, 0.5, 1))
        # ax.set_yticks((0, 0.2, 0.4, 0.6, 0.8, 1.0))
        plt.tight_layout()
        plt.savefig(os.path.join(
            save_path, f"{experiment_name}_{variable}_sims10.svg"), transparent=True)

        prct_replicated = np.mean(p_drug < 0.05, axis=0)*100
        replication_df = pd.DataFrame({"n_samples": n_samples, "replication": prct_replicated})
        fig, ax = plt.subplots(figsize=(1.3, 1.3))
        sns.barplot(x="n_samples", y="replication", data=replication_df, hue="n_samples", palette=new_cmap,
                    legend=False)
        ax.set_ylabel("Replication %")
        ax.set_yticks([0, 20, 40, 60, 80, 100])
        ax.set_ylim([0, 100])
        ax.set_xlabel("Sample size")
        plt.tight_layout()
        plt.savefig(os.path.join(
            save_path, f"{experiment_name}_{variable}_sampleXrep.svg"), transparent=True)
        plt.show()
