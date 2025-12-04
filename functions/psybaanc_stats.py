# -*- coding: utf-8 -*-
"""
Created on Tue Jul  1 15:40:52 2025

@author: olu
Functions used for psy-baanc statistical analyses

"""
# %% Import packages
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from scipy import stats
import scikit_posthocs as sp

labs = ["Stanford", "Berkeley 1", "Berkeley 2", "UCSF 1", "UCSF 2"]
pval_limits = [0.05, 0.01, 0.001, 0.0001]
pval_stars = ["ns", "*", "**", "***"]


# %% Divide data into groups
def get_group_data(data_final, col):
    """
    Divides the data into four groups based on
    Treatment (Saline vs. Psilocybin) and Sex (Male vs. Female).

    Parameters:
    - data_final : pd.DataFrame containing the data to be analyzed.
    - col : str. The name of the column in the DataFrame to be tested.
    Returns:
    - m_sal : list of values for male saline group.
    - f_sal : list of values for female saline group.
    - m_psi : list of values for male psilocybin group.
    - f_psi : list of values for female psilocybin group.
    """
    m_sal = data_final[(data_final["Sex"] == "M") & (data_final["Treatment"] == "S")][col].tolist()
    f_sal = data_final[(data_final["Sex"] == "F") & (data_final["Treatment"] == "S")][col].tolist()
    m_psi = data_final[(data_final["Sex"] == "M") & (data_final["Treatment"] == "P")][col].tolist()
    f_psi = data_final[(data_final["Sex"] == "F") & (data_final["Treatment"] == "P")][col].tolist()

    return m_sal, f_sal, m_psi, f_psi


def sample_size(data_y):
    """
    Calculates the sample size for each group and the total sample size.

    Parameters:
    - data_y : dict containing lists of values for each group.
    Returns:
    - statement : str summarizing the sample sizes for each group.
    - n_total : int total sample size across all groups.
    """

    m_sal_n = len(data_y["m_sal"])
    f_sal_n = len(data_y["f_sal"])
    m_psi_n = len(data_y["m_psi"])
    f_psi_n = len(data_y["f_psi"])

    n_total = m_sal_n + f_sal_n + m_psi_n + f_psi_n

    statement = ("Sample size: m, sal = " + str(m_sal_n) +
                 "; f, sal = " + str(f_sal_n) +
                 "; m, psi = " + str(m_psi_n) +
                 "; f, psi = " + str(f_psi_n))

    return statement, n_total


# %% Test for homogeneity
def levenes_test(data_y, n_total, groups=4):
    """
    Performs Levene's test for homogeneity of variances across multiple groups.

    Parameters:
    - data_y : dict containing lists of values for each group.
    - n_total : int total sample size across all groups.
    - groups : int number of groups being compared (default is 4).

    Returns:
    - statistic : float. Levene's test statistic.
    - pvalue : float. The p-value associated with the test statistic.
    - statement : str summarizing the test results.
    """

    statistic, pvalue = stats.levene(data_y["m_sal"], data_y["f_sal"],
                                     data_y["m_psi"], data_y["f_psi"])

    df1 = str(groups-1)
    df2 = str(n_total - groups)
    statement = f"Levene's test: F({df1},{df2})={round(statistic, 3)}, p={round(pvalue, 3)}"
    print(statement)

    return statistic, pvalue, statement


# %% Stats tests difference between means
def sm_ANOVA(formula, data):
    """
    Performs a two-way ANOVA using statsmodels.

    Parameters:
    - formula : str. The formula for the ANOVA model.
    - data : pd.DataFrame containing the data to be tested.
    Returns:
    - table : pd.DataFrame containing the ANOVA results.
    """

    ols_lm = smf.ols(formula, data=data)
    fit = ols_lm.fit()
    table = sm.stats.anova_lm(fit, typ=2)
    return table


def kruskal_wallis(data_y):
    """
    Performs the Kruskal-Wallis H-test for independent samples.
    Parameters:
    - data_y : dict containing lists of values for each group.
    Returns:
    - kruskal_df : pd.DataFrame containing the Kruskal-Wallis test results
    """

    h_statistic, p_value = stats.kruskal(data_y["m_sal"], data_y["f_sal"],
                                         data_y["m_psi"], data_y["f_psi"])

    kruskal_results = np.array((h_statistic, p_value)).reshape(1, -1)
    kruskal_df = pd.DataFrame(data=kruskal_results,
                              columns=["H", "PR(>H)"],
                              index=["Kruskal"])

    return kruskal_df


# %% Post hocs
def dunns(data_y):
    """
    Performs Dunn's post-hoc test with Holm's correction after a Kruskal-Wallis test.
    Parameters:
    - data_y : dict containing lists of values for each group.
    Returns:
    - holm_df : pd.DataFrame containing the Holm-corrected p-values for pairwise comparisons
    """

    dunn_results = sp.posthoc_dunn([data_y["m_sal"], data_y["f_sal"],
                                    data_y["m_psi"], data_y["f_psi"]], p_adjust=None)

    comparisons = [(1, 3), (2, 4), (1, 2), (3, 4)]
    pvals = [dunn_results.loc[g1, g2] for g1, g2 in comparisons]
    pvals_holm = multipletests(pvals, method='holm')[1].reshape(-1, 1)

    holm_df = pd.DataFrame(data=pvals_holm, columns=["p_corrected"],
                           index=["m,sal v. m,psi", "f,sal v. f,psi",
                                  "m,sal v. f,sal", "m,psi v. f,psi"])

    return holm_df


def sidaks(data_y):
    """
    Performs Sidak's post-hoc test after a two-way ANOVA with significant interaction.
    Parameters:
    - data_y : dict containing lists of values for each group.
    Returns:
    - results : pd.DataFrame containing the Sidak-corrected p-values for pairwise comparisons
    """

    sub_data = [data_y["m_sal"], data_y["f_sal"], data_y["m_psi"], data_y["f_psi"]]

    comparisons = [(0, 2), (1, 3), (0, 1), (2, 3)]
    comparison_names = ["m,sal v. m,psi", "f,sal v. f,psi", "m,sal v. f,sal", "m,psi v. f,psi"]

    t_stat = np.empty(len(comparisons))
    p_val = np.empty(len(comparisons))
    for comparison_no, comparison in enumerate(comparisons):
        g1 = sub_data[comparison[0]]
        g2 = sub_data[comparison[1]]
        t_stat[comparison_no], p_val[comparison_no] = stats.ttest_ind(g1, g2)

    pvals_sidak = multipletests(p_val, method='sidak')[1]

    results_dict = {
        "t_stat": t_stat,
        "p_val": p_val,
        "p_corrected": pvals_sidak}
    results = pd.DataFrame(data=results_dict, index=comparison_names)

    return results


# %% Asterisks
def get_stars(results):
    """
    Assigns significance stars based on p-values in the results DataFrame.
    Parameters:
    - results : pd.DataFrame containing the statistical test results with p-values.
    Returns:
    - stars : dict mapping each test to its significance star.
    """
    keys = results.index.tolist()
    vals = ["n/a"]*len(keys)
    stars = dict(zip(keys, vals))

    for col in stars.keys():
        pval_to_check = results.loc[col].values[-1]
        if "Residual" in col:
            stars[col] = np.nan
        else:
            stars[col] = next((pval_stars[x]
                              for x, val in enumerate(pval_limits) if val < pval_to_check),
                              pval_stars[-1])

    return stars


# %% Full analysis pipelines
def stats_treatment_sex(data, variable):
    """
    Performs statistical analysis for the effect of Treatment and Sex on a given variable.

    Parameters:
    - data : pd.DataFrame containing the data to be analyzed.
    - variable : str. The name of the column in the DataFrame to be tested.
    Returns:
    - results : pd.DataFrame containing the results of the statistical tests performed.
    """
    print(f"\n--------- Treatment*Sex analysis for: {variable} (summary data) -------------")
    cols_to_keep = ["Treatment", "Sex", "Institution"] + [variable]
    data = data[cols_to_keep]
    m_sal, f_sal, m_psi, f_psi = get_group_data(data, variable)
    data_y = {"m_sal": m_sal, "f_sal": f_sal, "m_psi": m_psi, "f_psi": f_psi}

    # first, print sample size and do levene's test for homogeneity of variances
    samples_statement, n_total = sample_size(data_y)
    print(samples_statement)
    levenes_p = levenes_test(data_y, n_total)[1]

    # if levenes failed, run kruskal wallis.
    if levenes_p < 0.05:
        kruskal_df = kruskal_wallis(data_y)
        stats_print = "Kruskal Wallis "

        # if kruskal significant, run dunns posthocs with holms corrections
        if kruskal_df.loc["Kruskal", "PR(>H)"] < 0.05:
            holm_df = dunns(data_y)
            kruskal_df = pd.concat([kruskal_df, holm_df])
            stats_print = stats_print + "w/ Dunn's posthocs and Holm corrections"
        results = kruskal_df

    # if levenes passed, run two-way ANOVA with treatment, sex, interactions as factors.
    else:
        formula = variable + ' ~ Treatment*Sex'
        results = sm_ANOVA(formula=formula, data=data)
        stats_print = "2-Way ANOVA (Treatment*Sex) "

        # if interaction significant, run sidaks posthoc.
        if results.loc["Treatment:Sex", "PR(>F)"] < 0.05:
            results_sidaks = sidaks(data_y)
            results = pd.concat([results, results_sidaks])
            stats_print = stats_print + "w/ Sidak's posthocs"

    results["significance"] = get_stars(results)

    print(stats_print)
    print(results)

    return results


def stats_treatment_sex_lab(data, variable):
    """
    Performs a mixed model analysis for the effect of Treatment and Sex on a given variable,
    with Institution as a random effect.

    Parameters:
    - data : pd.DataFrame containing the data to be analyzed.
    - variable : str. The name of the column in the DataFrame to be tested.

    Returns:
    - results : pd.DataFrame containing the results of the mixed model analysis.
    """

    print(f"\n--- Treatment*Sex Mixed model for: {variable} (Random effect: Institution) ---")

    data = organize_categories(data)
    model = smf.mixedlm(f"{variable} ~ Treatment*Sex",
                        data=data,
                        groups=data["Institution"])
    result = model.fit()
    results = result.summary().tables[1]
    results = results.iloc[:-1, :].astype(float)
    results = rename_results_indices(results)

    # if interaction is significant, run separate mixed models within subsets of data
    result_pvals = rename_results_indices(result.pvalues)

    if result_pvals["Treatment:Sex"] < 0.05:
        abridged_data = []
        abridged_data.append(data[data["Sex"] == "M"])
        abridged_data.append(data[data["Sex"] == "F"])
        abridged_data.append(data[data["Treatment"] == "S"])
        abridged_data.append(data[data["Treatment"] == "P"])

        pvals_posthocs = []
        for i in range(len(abridged_data)):
            if i < 2:
                model = smf.mixedlm(f"{variable} ~ Treatment",
                                    data=abridged_data[i],
                                    groups=abridged_data[i]["Institution"])

            if i >= 2:
                model = smf.mixedlm(f"{variable} ~ Sex",
                                    data=abridged_data[i],
                                    groups=abridged_data[i]["Institution"])

            result = model.fit()
            pvals_posthocs.append(result_pvals.iloc[1])

        p_adjusted = multipletests(pvals_posthocs, alpha=0.05, method='holm')[1]
        holm_df = pd.DataFrame(data=p_adjusted, columns=["p_corrected"],
                               index=["m,sal v. m,psi", "f,sal v. f,psi",
                                      "m,sal v. f,sal", "m,psi v. f,psi"])
        results = pd.concat([results, holm_df])

    results["significance"] = get_stars(results)
    print(results)

    return results


def stats_treatment_sex_third(data, variable, third_factor, within_subject=True):
    """
    Performs a mixed model analysis for the effect of Treatment, Sex,
    and a third factor on a given variable,
    with Animal_ID as a random effect.

    Parameters:
    - data : pd.DataFrame containing the data to be analyzed.
    - variable : str. The name of the column in the DataFrame to be tested.
    - third_factor : str. The name of the third factor column in the DataFrame to be tested.

    Returns:
    - results : pd.DataFrame containing the results of the mixed model analysis.
    """

    print(f"\n--------- 3-factor mixed model for: {variable} -------------")

    data = organize_categories(data)

    if within_subject:
        re_formula = f"~{third_factor}"
    else:
        re_formula = "~1"
    model = smf.mixedlm(f"{variable} ~ Treatment * Sex * {third_factor}",
                        data,
                        groups=data["Animal_ID"],
                        re_formula=re_formula)

    result = model.fit()
    results = result.summary().tables[1]
    results = results.iloc[:-1, :].astype(float)
    results = rename_results_indices(results)

    # Check whether posthocs are necessary
    result_pvals = rename_results_indices(result.pvalues)
    do_sex_third = result_pvals[f"Treatment:Sex:{third_factor}"] < 0.05
    do_sex = result_pvals["Treatment:Sex"] < 0.05
    do_third = result_pvals[f"Treatment:{third_factor}"] < 0.05

    pvals_posthocs = []
    group = []
    third_factor_vals = np.unique(data[third_factor])
    if do_sex_third:
        for sex in ["M", "F"]:
            for third_val in third_factor_vals:
                abridged_data = data[(data["Sex"] == sex) &
                                     (data[third_factor] == third_val)]
                sal = abridged_data[abridged_data["Treatment"] == "S"]
                psi = abridged_data[abridged_data["Treatment"] == "P"]
                p = stats.ttest_ind(sal[variable], psi[variable], equal_var=False)[1]
                pvals_posthocs.append(p)
                group.append(f"sal v. psi [{sex}_{third_val}]")
        p_adjusted = multipletests(pvals_posthocs, alpha=0.05, method='holm')[1]
        holm_df = pd.DataFrame(data=p_adjusted, columns=["p_corrected"],
                               index=group)
        results = pd.concat([results, holm_df])

    elif do_sex or do_third:
        if do_sex:
            for sex in ["M", "F"]:
                abridged_data = data[data["Sex"] == sex]
                sal = abridged_data[abridged_data["Treatment"] == "S"]
                psi = abridged_data[abridged_data["Treatment"] == "P"]
                t, p = stats.ttest_ind(sal[variable], psi[variable], equal_var=False)
                pvals_posthocs.append(p)
                group.append(f"sal v. psi [{sex}]")
        if do_third:
            for third_val in third_factor_vals:
                abridged_data = data[data[third_factor] == third_val]
                sal = abridged_data[abridged_data["Treatment"] == "S"]
                psi = abridged_data[abridged_data["Treatment"] == "P"]
                t, p = stats.ttest_ind(sal[variable], psi[variable], equal_var=False)
                pvals_posthocs.append(p)
                group.append(f"sal v. psi [{third_val}]")
        p_adjusted = multipletests(pvals_posthocs, alpha=0.05, method='holm')[1]
        holm_df = pd.DataFrame(data=p_adjusted, columns=["p_corrected"],
                               index=group)
        results = pd.concat([results, holm_df])

    results["significance"] = get_stars(results)
    print(results)

    return results


# %% Misc
def organize_categories(data):
    """
    Organizes the 'Treatment' and 'Sex' columns in the DataFrame as categorical variables
    with specified order.
    Parameters:
    - data : pd.DataFrame containing the data to be organized.
    Returns:
    - data : pd.DataFrame with organized categorical columns.
    """
    data['Treatment'] = data['Treatment'].astype('category')
    data['Treatment'] = data['Treatment'].cat.reorder_categories(['S', 'P'], ordered=True)

    data['Sex'] = data['Sex'].astype('category')
    data['Sex'] = data['Sex'].cat.reorder_categories(['M', 'F'], ordered=True)

    return data


def rename_results_indices(results):
    """
    Renames the indices of the results DataFrame for better readability.
    Parameters:
    - results : pd.DataFrame containing the statistical test results.
    Returns:
    - results : pd.DataFrame with renamed indices.
    """
    for i, index in enumerate(results.index.tolist()):
        rename_treatment = "Treatment[T.P]" in index
        rename_sex = "Sex[T.F]" in index
        rename_stress = "Stress[T.Stress]" in index

        if rename_treatment:
            index = index.replace("Treatment[T.P]", "Treatment")

        if rename_sex:
            index = index.replace("Sex[T.F]", "Sex")

        if rename_stress:
            index = index.replace("Stress[T.Stress]", "Stress")

        results = results.rename(index={results.index.tolist()[i]: index})

    return results


def get_zscores(data_expt, variable):
    """
    Standardizes psi data to sal data within each lab and sex.
    Parameters:
    - data_expt: pd.DataFrame with data for all labs.
    - variable: str, column name of the variable to standardize.
    Returns:
    - psi_data: pd.DataFrame with standardized psi data.
    """
    psi_data = []
    for lab in labs:
        psi_data_lab = pd.DataFrame(columns=["Institution", "Treatment", "Sex", variable])
        data_expt_lab = data_expt[data_expt["Institution"] == lab]

        # get data
        sal_male_data, sal_fem_data, psi_male_data, psi_fem_data = get_group_data(data_expt_lab,
                                                                                  variable)

        # get standardization values
        sal_male_mean = np.mean(sal_male_data)
        sal_male_std = np.std(sal_male_data, ddof=1)
        sal_fem_mean = np.mean(sal_fem_data)
        sal_fem_std = np.std(sal_fem_data, ddof=1)

        # standardize values
        zscores_male = (psi_male_data - sal_male_mean)/sal_male_std
        zscores_female = (psi_fem_data - sal_fem_mean)/sal_fem_std
        zscores_all = np.concat((zscores_male, zscores_female))

        id_info = np.array([lab, "P", "M"]*len(zscores_male) +
                           [lab, "P", "F"]*len(zscores_female)).reshape(len(zscores_all), -1)
        psi_data_lab[variable] = zscores_all
        psi_data_lab[["Institution", "Treatment", "Sex"]] = id_info
        psi_data.append(psi_data_lab)
    psi_data = pd.concat(psi_data)

    return psi_data
