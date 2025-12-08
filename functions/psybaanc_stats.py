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

labs = ["Stanford", "Berkeley 1", "Berkeley 2", "UCSF 1", "UCSF 2"]
pval_limits = [0.05, 0.01, 0.001, 0.0001]
pval_stars = ["ns", "*", "**", "***"]

pd.set_option('display.float_format', lambda x: f'{x:.3g}')


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

    n_dict = {}
    n_total = 0
    statement = "Sample size: "
    for key in data_y.keys():
        n_dict[key] = len(data_y[key])
        n_total = n_total + n_dict[key]
        statement = statement + f"{key} = {n_dict[key]}; "

    statement = statement[:-2]

    return statement, n_total


# %% Test for homogeneity of variances
def levenes_test(data_y, n_total):
    """
    Performs Levene's test for homogeneity of variances across multiple groups.

    Parameters:
    - data_y : dict containing lists of values for each group.
    - n_total : int total sample size across all groups.

    Returns:
    - statistic : float. Levene's test statistic.
    - pvalue : float. The p-value associated with the test statistic.
    - statement : str summarizing the test results.
    """

    statistic, pvalue = stats.levene(*data_y.values())

    groups = len(data_y)
    df1 = str(groups-1)
    df2 = str(n_total - groups)
    statement = f"Levene's test: F({df1},{df2})={round(statistic, 3)}, p={round(pvalue, 3)}"
    print(statement)

    return statistic, pvalue, statement


# %% Stats tests and posthocs
def sm_ANOVA(formula, data, parametric=True):
    """
    Performs a type 3 ANOVA using statsmodels.

    Parameters:
    - formula : str. The formula for the ANOVA model.
    - data : pd.DataFrame containing the data to be tested.
    Returns:
    - table : pd.DataFrame containing the ANOVA results.
    """

    ols_lm = smf.ols(formula, data=data).fit()
    if parametric:
        table = sm.stats.anova_lm(ols_lm, typ=3)
        table = table.rename(columns={'PR(>F)': 'p>|F|'})
    else:
        # HC3 robust SEs
        robust = ols_lm.get_robustcov_results(cov_type="HC3")
        # Wald-type tests for each main effect and interaction
        table = robust.wald_test_terms(scalar=True).table
        table = table.rename(columns={"statistic": "Chi^2", "pvalue": "p>|Chi^2|",
                                      "df_constraint": "df_num"})
        table = table[["df_num", "df_denom", "Chi^2", "p>|Chi^2|"]]

    return table


def two_factor_posthocs(data_y, parametric, method='holm'):
    """
    Performs post-hoc tests after a two-way ANOVA with significant interaction.
    Parameters:
    - data_y : dict containing lists of values for each group.
    Returns:
    - results : pd.DataFrame containing the Holm-corrected p-values for pairwise comparisons
    """

    sub_data = [data_y["m_sal"], data_y["f_sal"], data_y["m_psi"], data_y["f_psi"]]

    comparisons = [(0, 2), (1, 3), (0, 1), (2, 3)]
    comparison_names = ["M,S v. M,P", "F,S v. F,P", "M,S v. F,S", "M,P v. F,P"]

    stat = np.empty(len(comparisons))
    p_val = np.empty(len(comparisons))
    for comparison_no, comparison in enumerate(comparisons):
        g1 = sub_data[comparison[0]]
        g2 = sub_data[comparison[1]]
        if parametric:
            stat[comparison_no], p_val[comparison_no] = stats.ttest_ind(g1, g2, equal_var=True)
        else:
            stat[comparison_no], p_val[comparison_no] = stats.ttest_ind(g1, g2, equal_var=False)

    pvals_corrected = multipletests(p_val, method=method)[1]

    results_dict = {
        "t-statistic": stat,
        "p>|t|": p_val,
        "p_holm": pvals_corrected}
    results = pd.DataFrame(data=results_dict, index=comparison_names)

    return results


def three_factor_posthocs(results, data, variable, third_factor, parametric=True,
                          third_factor_sex=False):
    """
    Performs post-hoc tests after a three-factor analysis.
    Parameters:
    - reults: pd.DataFrame containing p-values for various factors and interactions
    - data: pd.DataFrame containing data that the results are derived from.
    - variable: name of the variable, which is a column in data, you are analyzing.
    - third_factor: str of the additional variable to analyze (outside of Treatment and Sex)
    Returns:
    - holm_df : pd.DataFrame containing the Holm-corrected p-values for pairwise comparisons
    """
    stat_posthocs = []
    pvals_posthocs = []
    group = []
    third_factor_vals = np.unique(data[third_factor])
    if parametric:
        equal_var = True
    else:
        equal_var = False

    # Check which posthocs to do
    treatment_sex_third = results.loc[f"Treatment:Sex:{third_factor}"] < 0.05
    treatment_sex = results.loc["Treatment:Sex"] < 0.05
    treatment_third = results.loc[f"Treatment:{third_factor}"] < 0.05
    if third_factor_sex:
        third_sex = results.loc[f"Sex:{third_factor}"] < 0.05
    else:
        third_sex = False

    if treatment_sex_third:
        for sex in ["M", "F"]:
            for third_val in third_factor_vals:
                abridged_data = data[(data["Sex"] == sex) &
                                     (data[third_factor] == third_val)]
                sal = abridged_data[abridged_data["Treatment"] == "S"]
                psi = abridged_data[abridged_data["Treatment"] == "P"]
                stat, p = stats.ttest_ind(sal[variable], psi[variable], equal_var=equal_var)
                stat_posthocs.append(stat)
                pvals_posthocs.append(p)
                group.append(f"{third_val},S v. {third_val},P [{sex}]")
        if third_factor_sex:
            for sex in ["M", "F"]:
                for treatment in ["S", "P"]:
                    abridged_data = data[(data["Sex"] == sex) &
                                         (data["Treatment"] == treatment)]
                    third_1 = abridged_data[abridged_data[third_factor] == third_factor_vals[0]]
                    third_2 = abridged_data[abridged_data[third_factor] == third_factor_vals[1]]
                    stat, p = stats.ttest_ind(third_1[variable], third_2[variable],
                                              equal_var=equal_var)
                    stat_posthocs.append(stat)
                    pvals_posthocs.append(p)
                    group.append(f"{third_factor_vals[0]},{treatment} v. "
                                 f"{third_factor_vals[1]},{treatment} [{sex}]")
        p_adjusted = multipletests(pvals_posthocs, alpha=0.05, method='holm')[1]
        holm_df = pd.DataFrame(data=list(zip(stat_posthocs, pvals_posthocs, p_adjusted)),
                               columns=["t-statistic", "p>|t|", "p_holm"],
                               index=group)

    elif treatment_sex or treatment_third or third_sex:
        if treatment_sex:
            for sex in ["M", "F"]:
                abridged_data = data[data["Sex"] == sex]
                sal = abridged_data[abridged_data["Treatment"] == "S"]
                psi = abridged_data[abridged_data["Treatment"] == "P"]
                stat, p = stats.ttest_ind(sal[variable], psi[variable], equal_var=equal_var)
                stat_posthocs.append(stat)
                pvals_posthocs.append(p)
                group.append(f"{sex},S v. {sex},P")
        if treatment_third:
            for third_val in third_factor_vals:
                abridged_data = data[data[third_factor] == third_val]
                sal = abridged_data[abridged_data["Treatment"] == "S"]
                psi = abridged_data[abridged_data["Treatment"] == "P"]
                stat, p = stats.ttest_ind(sal[variable], psi[variable], equal_var=equal_var)
                stat_posthocs.append(stat)
                pvals_posthocs.append(p)
                group.append(f"{third_val},S v. {third_val},P")
        if third_sex:
            for sex in ["M", "F"]:
                abridged_data = data[data["Sex"] == sex]
                third_1 = abridged_data[abridged_data[third_factor] == third_factor_vals[0]]
                third_2 = abridged_data[abridged_data[third_factor] == third_factor_vals[1]]
                stat, p = stats.ttest_ind(third_1[variable], third_2[variable], equal_var=equal_var)
                stat_posthocs.append(stat)
                pvals_posthocs.append(p)
                group.append(f"{sex},{third_factor_vals[0]} v. {sex},{third_factor_vals[1]}")
        p_adjusted = multipletests(pvals_posthocs, alpha=0.05, method='holm')[1]
        holm_df = pd.DataFrame(data=list(zip(stat_posthocs, pvals_posthocs, p_adjusted)),
                               columns=["t-statistic", "p>|t|", "p_holm"],
                               index=group)

    if len(pvals_posthocs) == 0:
        return None
    else:
        return holm_df


# %% Asterisks
def get_stars(results, parametric=True):
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
        end_col = -1
        pval_to_check = np.nan
        while np.isnan(pval_to_check):
            pval_to_check = results.loc[col].values[end_col]
            end_col = end_col - 1

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
    institution = data["Institution"][0]
    print(f"\n------------- {institution}: {variable} (summary data) -------------")

    cols_to_keep = ["Treatment", "Sex", "Institution"] + [variable]
    data = data[cols_to_keep]
    m_sal, f_sal, m_psi, f_psi = get_group_data(data, variable)
    data_y = {"m_sal": m_sal, "f_sal": f_sal, "m_psi": m_psi, "f_psi": f_psi}

    # first, print sample size and do levene's test for homogeneity of variances
    samples_statement, n_total = sample_size(data_y)
    print(samples_statement)
    levenes_p = levenes_test(data_y, n_total)[1]

    # if levenes failed.
    if levenes_p < 0.05:
        parametric = False
        stats_print = "2-Way Robust ANOVA (Treatment*Sex) (Wald tests using HC3 covariance)"

    else:
        parametric = True
        stats_print = "2-Way ANOVA (Treatment*Sex)"

    # Run two-way ANOVA or non-parametric variant
    formula = variable + ' ~ C(Treatment, Sum) * C(Sex, Sum)'
    results = sm_ANOVA(formula=formula, data=data, parametric=parametric)
    results = rename_results_indices(results)

    # if interaction significant, run posthocs with holms corrections
    if results.loc["Treatment:Sex", results.columns[-1]] < 0.05:
        posthocs_df = two_factor_posthocs(data_y, parametric=parametric)
        results = pd.concat([results, posthocs_df])
        if parametric:
            stats_print = stats_print + "\nPosthocs: T-test + Holm's corrections"
        else:
            stats_print = stats_print + "\nPosthocs: Welch's test + Holm's corrections"

    results["significance"] = get_stars(results)
    results = results.fillna("")

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
    print(f"\n------------- All labs: {variable} (summary data) -------------")

    stats_print = f"Linear Mixed Effects Model: {variable}~Treatment*Sex+(1|Institution)"
    data = organize_categories(data)
    model = smf.mixedlm(f"{variable} ~ Treatment*Sex",
                        data=data,
                        groups=data["Institution"])

    result = model.fit()
    results = result.summary().tables[1]
    results.replace('', np.nan, inplace=True)
    results = results.dropna()
    results = results.iloc[:, :-2].astype(float)
    results = rename_results_indices(results)

    # if interaction is significant, run separate mixed models within subsets of data
    result_pvals = rename_results_indices(result.pvalues)
    if result_pvals["Treatment:Sex"] < 0.05:
        stats_print = stats_print + ("\nPosthocs: Linear Mixed Effects Model: "
                                     f"\n{variable}~Treatment+(1|Institution) in sex-separated data"
                                     f"\n{variable}~Sex+(1|Institution) in drug-separated data")
        abridged_data = []
        abridged_data.append(data[data["Sex"] == "M"])
        abridged_data.append(data[data["Sex"] == "F"])
        abridged_data.append(data[data["Treatment"] == "S"])
        abridged_data.append(data[data["Treatment"] == "P"])

        pvals_posthocs = []
        stat_posthocs = []
        for i in range(len(abridged_data)):
            if i < 2:
                model = smf.mixedlm(f"{variable} ~ Treatment",
                                    data=abridged_data[i],
                                    groups=abridged_data[i]["Institution"])

            if i >= 2:
                model = smf.mixedlm(f"{variable} ~ Sex",
                                    data=abridged_data[i],
                                    groups=abridged_data[i]["Institution"])

            results_posthoc = model.fit()
            pvals_posthocs.append(results_posthoc.iloc[:, -1])
            stat_posthocs.append(results_posthoc.iloc[:, -2])

        p_adjusted = multipletests(pvals_posthocs, alpha=0.05, method='holm')[1]
        holm_df = pd.DataFrame(data=p_adjusted, columns=["z", "p>|z|", "p_holm"],
                               index=["M,S v. M,P", "F,S v. F,P",
                                      "M,S v. F,S", "M,P v. F,P"])
        results = pd.concat([results, holm_df])

    results["significance"] = get_stars(results)
    results = results.fillna("")
    print(stats_print)
    print(results)

    return results


def stats_treatment_sex_third(data, variable, third_factor, third_factor_sex=False):
    """
    Performs a mixed model analysis for the effect of Treatment, Sex,
    and a third factor + all interactions on a given variable.
    Animal ID is included as a random effect with a random slope for [third factor].

    Parameters:
    - data : pd.DataFrame containing the data to be analyzed.
    - variable : str. The name of the column in the DataFrame to be tested.
    - third_factor : str. The name of the third factor column in the DataFrame to be tested.

    Returns:
    - results : pd.DataFrame containing the results of the mixed model analysis.
    """

    institution = data["Institution"][0]
    print(f"\n------------- {institution}: {variable} -------------")

    data = organize_categories(data)
    data_y = {}
    third_factor_vals = np.unique(data[third_factor])
    for val in third_factor_vals:
        group_names = [f"m_sal_{val}", f"f_sal_{val}", f"m_psi_{val}", f"f_psi_{val}"]
        (data_y[group_names[0]], data_y[group_names[1]],
         data_y[group_names[2]], data_y[group_names[3]]) = (
             get_group_data(data[data[third_factor] == val], variable))

    # first, print sample size and do levene's test for homogeneity of variances
    samples_statement, n_total = sample_size(data_y)
    print(samples_statement)
    levenes_p = levenes_test(data_y, n_total)[1]
    parametric = levenes_p > 0.05

    if len(third_factor_vals) <= 2:
        re_formula = "~1"
        stats_print = ("\nLinear Mixed Effects Model: "
                       f"{variable}~Treatment*Sex*{third_factor}+(1+|Animal ID)")
    else:
        re_formula = f"~{third_factor}"
        stats_print = ("\nLinear Mixed Effects Model: "
                       f"{variable}~Treatment*Sex*{third_factor}+(1+{third_factor}|Animal ID)")
    model = smf.mixedlm(f"{variable} ~ Treatment * Sex * {third_factor}",
                        data,
                        groups=data["Animal_ID"],
                        re_formula=re_formula)

    result = model.fit()
    results = result.summary().tables[1]
    results.replace('', np.nan, inplace=True)
    results = results.dropna()
    results = results.iloc[:, :-2].astype(float)
    results = rename_results_indices(results)

    # Do necessary posthocs
    result_pvals = rename_results_indices(result.pvalues)
    holm_df = three_factor_posthocs(result_pvals, data, variable,
                                    third_factor, parametric, third_factor_sex)
    if holm_df is not None:
        results = pd.concat([results, holm_df])
        if parametric:
            stats_print = stats_print + "\nPosthocs: T-test + Holm's corrections"
        else:
            stats_print = stats_print + "\nPosthocs: Welch's test + Holm's corrections"

    results["significance"] = get_stars(results)
    results = results.fillna("")
    print(stats_print)
    print(results)

    return results


def stats_treatment_sex_stress(data, variable):
    """
     Performs statistical analysis for the effect of Treatment, Sex, and Stress on a given variable.

     Parameters:
     - data : pd.DataFrame containing the data to be analyzed.
     - variable : str. The name of the column in the DataFrame to be tested.
     Returns:
     - results : pd.DataFrame containing the results of the statistical tests performed.
    """

    institution = data["Institution"][0]
    print(f"\n------------- {institution}: {variable} -------------")

    data = organize_categories(data)
    cols_to_keep = ["Treatment", "Sex", "Stress"] + [variable]
    data = data[cols_to_keep]

    data_ctrl = data[data["Stress"] == "Ctrl"]
    data_stress = data[data["Stress"] == "Stress"]

    m_sal_ctrl, f_sal_ctrl, m_psi_ctrl, f_psi_ctrl = get_group_data(data_ctrl, variable)
    m_sal_stress, f_sal_stress, m_psi_stress, f_psi_stress = get_group_data(data_stress, variable)

    data_y = {"m_sal_ctrl": m_sal_ctrl, "f_sal_ctrl": f_sal_ctrl,
              "m_psi_ctrl": m_psi_ctrl, "f_psi_ctrl": f_psi_ctrl,
              "m_sal_stress": m_sal_stress, "f_sal_stress": f_sal_stress,
              "m_psi_stress": m_psi_stress, "f_psi_stress": f_psi_stress}

    # first, print sample size and do levene's test for homogeneity of variances
    samples_statement, n_total = sample_size(data_y)
    print(samples_statement)
    levenes_p = levenes_test(data_y, n_total)[1]

    # if levenes failed.
    if levenes_p < 0.05:
        parametric = False
        stats_print = "3-Way Robust ANOVA (Treatment*Sex*Stress) (Wald tests using HC3 covariance) "

    else:
        parametric = True
        stats_print = "3-Way ANOVA (Treatment*Sex*Stress) "

    # Set up the model
    formula = variable + ' ~ C(Treatment, Sum)*C(Sex, Sum)*C(Stress, Sum)'
    results = sm_ANOVA(formula=formula, data=data, parametric=parametric)
    results = rename_results_indices(results)
    result_pvals = results[results.columns[-1]]
    holm_df = three_factor_posthocs(result_pvals, data, variable, "Stress", parametric,
                                    third_factor_sex=True)
    if holm_df is not None:
        results = pd.concat([results, holm_df])
        if parametric:
            stats_print = stats_print + "\nPosthocs: T-test + Holm's Corrections"
        else:
            stats_print = stats_print + "\nPosthocs: Welch's Test + Holm's Corrections"
    results["significance"] = get_stars(results)
    results = results.fillna("")
    print(stats_print)
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

        if "Treatment[T.P]" in index:
            index = index.replace("Treatment[T.P]", "Treatment")

        if "Sex[T.F]" in index:
            index = index.replace("Sex[T.F]", "Sex")

        if "Stress[T.Stress]" in index:
            index = index.replace("Stress[T.Stress]", "Stress")

        if "C(Treatment, Sum)" in index:
            index = index.replace("C(Treatment, Sum)", "Treatment")

        if "C(Sex, Sum)" in index:
            index = index.replace("C(Sex, Sum)", "Sex")

        if "C(Stress, Sum)" in index:
            index = index.replace("C(Stress, Sum)", "Stress")

        if "Zone_type[T.social]" in index:
            index = index.replace("Zone_type[T.social]", "Zone_type")

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
