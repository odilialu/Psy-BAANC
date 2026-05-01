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
import marginaleffects as me

labs = ["Stanford", "Berkeley 1", "Berkeley 2", "UCSF 1", "UCSF 2"]
pval_limits = [0.05, 0.01, 0.001, 0.0001]
pval_stars = ["ns", "*", "**", "***"]
pd.set_option('display.float_format', lambda x: f'{x:.3g}')


# %% Get sample sizes
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
    statement = statement + " mice"

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

    return statistic, pvalue, statement


# %% Stats tests and posthocs pipelines
def students_T(g1, g2, results, comparison_name, parametric=True, paired=False):
    """ Performs a T-test and returns critical stats, i.e., t-statistic w/ 95% confidence interval,
    df, effect size in Hedges g, and pvalue.

    Parameters:
    - g1 : group 1 values
    - g2 : group 2 values
    - results : results dataframe to add new comparisons to
    - comparison_name : name of the comparison, used as index added to results df.
    - parametric : True / False. For independent groups, Student's T-test if True. Otherwise Welch's.
    - paired : True / False. Option to do a paired T-test. 

    Returns:
    - results : dataframe containing results for the comparison. 
    """

    if paired:
        result = stats.ttest_rel(g2, g1)
        comparison_name = comparison_name + " (Paired)"
    else:
        if parametric:
            result = stats.ttest_ind(g2, g1, equal_var=True)
        else:
            result = stats.ttest_ind(g2, g1, equal_var=False)

    results.loc[comparison_name, "t_df"] = result.df
    results.loc[comparison_name, "t_statistic"] = result.statistic
    ci = result.confidence_interval(confidence_level=0.95)
    results.loc[comparison_name, "t_0.025"] = ci.low
    results.loc[comparison_name, "t_0.975"] = ci.high
    results.loc[comparison_name, "Hedges' g"] = calculate_hedges_g(g2, g1)
    results.loc[comparison_name, "P>|t|"] = result.pvalue
    return results


def sm_ANOVA(formula, data, parametric=True):
    """
    Performs a type 3 ANOVA using statsmodels.

    Parameters:
    - formula : str. The formula for the ANOVA model.
    - data : pd.DataFrame containing the data to be tested.
    - parametric : True / False. If true, regular ANOVAs are used.
        - If False, robust ANOVA is used (Wald tests using HC3 covariance')
    Returns:
    - results : pd.DataFrame containing the ANOVA results.
    """

    ols_lm = smf.ols(formula, data=data).fit()
    if parametric:
        results = sm.stats.anova_lm(ols_lm, typ=3)
        results = results.rename(columns={'PR(>F)': 'P>|F|'})
        coefficients = ols_lm.params
        results["Coef."] = np.append(coefficients.values, np.nan)
        ci_intervals = ols_lm.conf_int(alpha=0.05)
        results["[0.025"] = np.append(ci_intervals[0].values, np.nan)
        results["0.975]"] = np.append(ci_intervals[1].values, np.nan)
        results["partial_eta^2"] = (results["sum_sq"]/(results["sum_sq"] +
                                                       results.loc["Residual", "sum_sq"]))
        results = results[["Coef.", "[0.025", "0.975]",
                           "sum_sq", "partial_eta^2", "df", "F", "P>|F|"]]

    else:  # HC3 robust SEs
        # Wald-type tests for each main effect and interaction
        robust = ols_lm.get_robustcov_results(cov_type="HC3")
        results = robust.wald_test_terms(scalar=True).table
        results = results.rename(columns={"statistic": "Chi^2", "pvalue": "P>|Chi^2|",
                                          "df_constraint": "df_num"})
        results["Coef."] = robust.params
        conf_intervals = robust.conf_int(alpha=0.05)
        results["[0.025"] = conf_intervals[:, 0]
        results["0.975]"] = conf_intervals[:, 1]
        temp_results = sm.stats.anova_lm(ols_lm, typ=3)
        results["partial_eta^2"] = (temp_results["sum_sq"] /
                                    (temp_results["sum_sq"] + temp_results.loc["Residual", "sum_sq"]))
        results = results[["Coef.", "[0.025", "0.975]", "partial_eta^2",
                           "df_num", "df_denom", "Chi^2", "P>|Chi^2|"]]
    return results


def two_factor_posthocs(data_y, parametric=True, method='holm'):
    """
    Performs post-hoc tests after a two-way ANOVA with significant interaction.
    Parameters:
    - data_y : dict containing lists of values for each group.
    - parametric : True / False
    Returns:
    - results : pd.DataFrame containing the Holm-corrected p-values for pairwise comparisons
    """

    sub_data = [data_y["m_sal"], data_y["f_sal"],
                data_y["m_psi"], data_y["f_psi"]]
    comparisons = [(0, 2), (1, 3), (0, 1), (2, 3)]
    comparison_names = ["M,S v. M,P", "F,S v. F,P", "M,S v. F,S", "M,P v. F,P"]

    results = pd.DataFrame(index=comparison_names, columns=[
        "t_df", "t_statistic", "t_0.025", "t_0.975", "Hedges' g",
        "P>|t|", "p_holm"])
    for comparison_no, comparison in enumerate(comparisons):
        g1 = sub_data[comparison[0]]
        g2 = sub_data[comparison[1]]
        comparison_name = comparison_names[comparison_no]
        results = students_T(g1, g2, results, comparison_name, parametric=parametric)

    results["p_holm"] = multipletests(results["P>|t|"], method=method)[1]

    return results


def three_factor_posthocs(result, results, data, variable, third_factor, parametric=True,
                          third_factor_compare=True, third_paired=True, all_labs=False):
    """
    Performs post-hoc tests after a three-factor analysis.
    Parameters:
    - results: pd.DataFrame containing p-values for various factors and interactions
    - data: pd.DataFrame containing data that the results are derived from.
    - variable: name of the variable, which is a column in data, you are analyzing.
    - third_factor: str of the additional variable to analyze (outside of treatment and sex)
    Returns:
    - holm_df : pd.DataFrame containing the Holm-corrected p-values for pairwise comparisons
    """
    third_factor_vals = np.unique(data[third_factor])

    # Set up posthocs dataframe.
    if all_labs:
        holm_df = pd.DataFrame(columns=["Contrast Est.", "Contrast SE", "Contrast_0.025", "Contrast_0.975",
                                        "Contrast z", "Contrast P>|z|"])
    else:
        holm_df = pd.DataFrame(columns=[
            "t_df", "t_statistic", "t_0.025", "t_0.975", "Hedges' g",
            "P>|t|"])

    # Check which posthocs to do
    treatment_sex_third = results.loc[f"treatment:sex:{third_factor}"] < 0.05
    treatment_sex = results.loc["treatment:sex"] < 0.05
    treatment_third = results.loc[f"treatment:{third_factor}"] < 0.05
    if third_factor_compare:
        third_sex = results.loc[f"sex:{third_factor}"] < 0.05
    else:
        third_sex = False

    # Do posthoc comparisons
    if treatment_sex_third:
        if all_labs:
            posthoc_treatment = pd.DataFrame(me.comparisons(
                result, variables={"treatment": ["S", "P"]}, by=[f"{third_factor}", "sex"]))
        for sex in ["M", "F"]:
            for third_val in third_factor_vals:
                comparison_name = f"{third_val},S v. {third_val},P [{sex}]"
                if all_labs:
                    data_temp = posthoc_treatment[(posthoc_treatment[0] == third_val) & (
                        posthoc_treatment[1] == sex)].iloc[:, 4:11]
                    holm_df.loc[comparison_name] = data_temp.reindex(
                        columns=[4, 5, 9, 10, 6, 7]).values[0]
                else:
                    abridged_data = data[(data["sex"] == sex) &
                                         (data[third_factor] == third_val)]
                    sal = abridged_data[abridged_data["treatment"] == "S"]
                    psi = abridged_data[abridged_data["treatment"] == "P"]
                    holm_df = students_T(sal[variable], psi[variable],
                                         holm_df, comparison_name, parametric=parametric)
        if third_factor_compare:
            if all_labs:
                posthoc_third = pd.DataFrame(me.comparisons(
                    result, variables={f"{third_factor}": [
                        f"{third_factor_vals[0]}", f"{third_factor_vals[1]}"]},
                    by=["treatment", "sex"]))
            for sex in ["M", "F"]:
                for treatment in ["S", "P"]:
                    comparison_name = f"{third_factor_vals[0]},{treatment} v. {third_factor_vals[1]},{treatment} [{sex}]"
                    if all_labs:
                        data_temp = posthoc_third[(posthoc_third[0] == treatment) & (
                            posthoc_third[1] == sex)].iloc[:, 4:11]
                        holm_df.loc[comparison_name] = data_temp.reindex(
                            columns=[4, 5, 9, 10, 6, 7]).values[0]
                    else:
                        abridged_data = data[(data["sex"] == sex) &
                                             (data["treatment"] == treatment)]
                        third_1 = abridged_data[abridged_data[third_factor]
                                                == third_factor_vals[0]]
                        third_2 = abridged_data[abridged_data[third_factor]
                                                == third_factor_vals[1]]
                        holm_df = students_T(third_1[variable], third_2[variable],
                                             holm_df, comparison_name, parametric=parametric, paired=third_paired)

    elif treatment_sex or treatment_third or third_sex:
        if treatment_sex:
            if all_labs:
                posthoc_treatment = pd.DataFrame(me.comparisons(
                    result, variables={"treatment": ["S", "P"]}, by=["sex"]))
            for sex in ["M", "F"]:
                comparison_name = f"{sex},S v. {sex},P"
                if all_labs:
                    data_temp = posthoc_treatment[(posthoc_treatment[0] == sex)].iloc[:, 3:10]
                    holm_df.loc[comparison_name] = data_temp.reindex(
                        columns=[3, 4, 8, 9, 5, 6]).values[0]
                else:
                    abridged_data = data[data["sex"] == sex]
                    sal = abridged_data[abridged_data["treatment"] == "S"]
                    psi = abridged_data[abridged_data["treatment"] == "P"]
                    holm_df = students_T(sal[variable], psi[variable],
                                         holm_df, comparison_name, parametric=parametric)
        if treatment_third:
            if all_labs:
                posthoc_treatment = pd.DataFrame(me.comparisons(
                    result, variables={"treatment": ["S", "P"]}, by=[f"{third_factor}"]))
            for third_val in third_factor_vals:
                comparison_name = f"{third_val},S v. {third_val},P"
                if all_labs:
                    data_temp = posthoc_treatment[(posthoc_treatment[0] == third_val)].iloc[:, 3:10]
                    holm_df.loc[comparison_name] = data_temp.reindex(
                        columns=[3, 4, 8, 9, 5, 6]).values[0]
                else:
                    abridged_data = data[data[third_factor] == third_val]
                    sal = abridged_data[abridged_data["treatment"] == "S"]
                    psi = abridged_data[abridged_data["treatment"] == "P"]
                    holm_df = students_T(sal[variable], psi[variable],
                                         holm_df, comparison_name, parametric=parametric)
            if third_factor_compare:
                if all_labs:
                    posthoc_third = pd.DataFrame(me.comparisons(
                        result, variables={f"{third_factor}": [
                            f"{third_factor_vals[0]}", f"{third_factor_vals[1]}"]},
                        by=["treatment"]))
                for treatment in ["S", "P"]:  # effect of third depends on treatment.
                    comparison_name = f"{treatment},{third_factor_vals[0]} v. {treatment},{third_factor_vals[1]}"
                    if all_labs:
                        data_temp = posthoc_third[(posthoc_third[0] == treatment)].iloc[:, 3:10]
                        holm_df.loc[comparison_name] = data_temp.reindex(
                            columns=[3, 4, 8, 9, 5, 6]).values[0]
                    else:
                        abridged_data = data[data["treatment"] == treatment]
                        third_1 = abridged_data[abridged_data[third_factor] == third_factor_vals[0]]
                        third_2 = abridged_data[abridged_data[third_factor] == third_factor_vals[1]]
                        holm_df = students_T(third_1[variable], third_2[variable],
                                             holm_df, comparison_name, parametric=parametric, paired=third_paired)
        if third_sex:
            if all_labs:
                posthoc_third = pd.DataFrame(me.comparisons(
                    result, variables={f"{third_factor}": [
                        f"{third_factor_vals[0]}", f"{third_factor_vals[1]}"]},
                    by=["sex"]))
            for sex in ["M", "F"]:
                comparison_name = f"{sex},{third_factor_vals[0]} v. {sex},{third_factor_vals[1]}"
                if all_labs:
                    data_temp = posthoc_third[(posthoc_third[0] == sex)].iloc[:, 3:10]
                    holm_df.loc[comparison_name] = data_temp.reindex(
                        columns=[3, 4, 8, 9, 5, 6]).values[0]
                else:
                    abridged_data = data[data["sex"] == sex]
                    third_1 = abridged_data[abridged_data[third_factor]
                                            == third_factor_vals[0]]
                    third_2 = abridged_data[abridged_data[third_factor]
                                            == third_factor_vals[1]]
                    holm_df = students_T(third_1[variable], third_2[variable],
                                         holm_df, comparison_name, parametric=parametric, paired=third_paired)

    if len(holm_df) == 0:
        return None
    else:
        holm_df["p_holm"] = multipletests(holm_df.iloc[:, -1], alpha=0.05, method='holm')[1]
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
            if end_col < (len(results.loc[col].values)*-1):
                pval_to_check = 1
                break

        if col in ["Residual", "Group Var"] or "Var" in col or "Cov" in col:
            stars[col] = np.nan
        else:
            stars[col] = next((pval_stars[x]
                              for x, val in enumerate(pval_limits) if val < pval_to_check),
                              pval_stars[-1])

    return stars


# %% Full analysis pipelines
def stats_treatment_sex(data, variable):
    """
    Performs statistical analysis for the effect of treatment and sex on a given variable.

    Parameters:
    - data : pd.DataFrame containing the data to be analyzed.
    - variable : str. The name of the column in the DataFrame to be tested.
    Returns:
    - results : pd.DataFrame containing the results of the statistical tests performed.
    - results_print : list of information relevant to the statistical test
    """

    results_print = []
    institution = data["institution"].iloc[0]
    analysis_name = f"------------- {institution}: {variable} (summary data) -------------"
    results_print.append(f"\n{analysis_name}")
    data = organize_categories(data)

    cols_to_keep = ["treatment", "sex", "institution"] + [variable]
    data = data[cols_to_keep]
    m_sal, f_sal, m_psi, f_psi = get_group_data(data, variable)
    data_y = {"m_sal": m_sal, "f_sal": f_sal, "m_psi": m_psi, "f_psi": f_psi}

    # first, print sample size and do levene's test for homogeneity of variances
    samples_statement, n_total = sample_size(data_y)
    results_print.append(samples_statement)
    levenes_p, levenes_statement = levenes_test(data_y, n_total)[1:]
    results_print.append(levenes_statement)

    # if levenes failed.
    if levenes_p < 0.05:
        parametric = False
        results_print.append(
            "2-Way Robust ANOVA (treatment*sex) (Wald tests using HC3 covariance). ")

    else:
        parametric = True
        results_print.append("2-Way ANOVA (treatment*sex). ")

    # Run two-way ANOVA or non-parametric variant
    formula = variable + ' ~ C(treatment, Sum)*C(sex, Sum)'
    results = sm_ANOVA(formula=formula, data=data, parametric=parametric)
    results = rename_results_indices(results)

    # if interaction significant, run posthocs with holms corrections
    if results.loc["treatment:sex", results.columns[-1]] < 0.05:
        posthocs_df = two_factor_posthocs(data_y, parametric=parametric)
        results = pd.concat([results, posthocs_df])
        if parametric:
            results_print.append("Posthocs: T-test + Holm's corrections. ")
        else:
            results_print.append("Posthocs: Welch's test + Holm's corrections. ")

    results["significance"] = get_stars(results)
    results = results.fillna("")

    for line in results_print:
        print(line)
    print(results)

    return results, results_print


def stats_treatment_sex_lab(data, variable):
    """
    Performs a mixed model analysis for the effect of treatment and sex on a given variable,
    with institution as a random effect.

    Parameters:
    - data : pd.DataFrame containing the data to be analyzed.
    - variable : str. The name of the column in the DataFrame to be tested.

    Returns:
    - results : pd.DataFrame containing the results of the mixed model analysis.
    """
    results_print = []
    analysis_name = f"------------- All labs: {variable} (summary data) -------------"
    results_print.append(f"\n{analysis_name}")

    # Get sample sizes and perform Levene's
    cols_to_keep = ["mouse_ID", "treatment", "sex", "institution"] + [variable]
    data = data[cols_to_keep]
    m_sal, f_sal, m_psi, f_psi = get_group_data(data, variable)
    data_y = {"m_sal": m_sal, "f_sal": f_sal, "m_psi": m_psi, "f_psi": f_psi}
    samples_statement, n_total = sample_size(data_y)
    results_print.append(samples_statement)
    levenes_p, levenes_statement = levenes_test(data_y, n_total)[1:]
    results_print.append(levenes_statement)
    if levenes_p < 0.05:
        parametric = False
    else:
        parametric = True

    data = organize_categories(data, lme=True)
    model = smf.mixedlm(f"{variable} ~ C(treatment, Sum)*C(sex, Sum)",
                        data=data,
                        groups=data["institution"])

    result = model.fit(method="lbfgs", reml=True)
    results_print.append(f"LME model: {variable}~treatment*sex+(1|institution). "
                         f"llf (REML): {result.llf}. Converged: {result.converged}. ")
    results = result.summary().tables[1]
    results.replace('', np.nan, inplace=True)
    results = results[["Coef.", "Std.Err.", "[0.025", "0.975]", "z", "P>|z|"]]
    results = rename_results_indices(results).astype(float)

    # if interaction is significant, run separate mixed models within subsets of data
    if float(results.loc["treatment:sex", "P>|z|"]) < 0.05:
        if float(results.loc["Group Var", "Coef."]) < 0.1:
            holm_df = two_factor_posthocs(data_y, parametric=parametric, method='holm')
            if parametric:
                results_print.append("Posthocs: T-test + Holm's corrections. ")
            else:
                results_print.append("Posthocs: Welch's test + Holm's corrections. ")
        else:
            results_print.append(
                "Posthocs: Contrasts of estimated marginal means + Holm's corrections")
            comparisons = ["M,S v. M,P", "F,S v. F,P", "M,S v. F,S", "M,P v. F,P"]
            posthoc_treatment = pd.DataFrame(me.comparisons(
                result, variables={"treatment": ["S", "P"]}, by="sex"))
            posthoc_sex = pd.DataFrame(me.comparisons(
                result, variables={"sex": ["M", "F"]}, by="treatment"))
            males_SvP = posthoc_treatment[posthoc_treatment[0] == "M"].iloc[:, 3:10]
            females_SvP = posthoc_treatment[posthoc_treatment[0] == "F"].iloc[:, 3:10]
            sal_MvF = posthoc_sex[posthoc_sex[0] == "S"].iloc[:, 3:10]
            psi_MvF = posthoc_sex[posthoc_sex[0] == "P"].iloc[:, 3:10]
            holm_df = pd.DataFrame(columns=["Contrast Est.", "Contrast SE", "Contrast_0.025", "Contrast_0.975",
                                            "z", "P>|z|"], index=comparisons)
            holm_df.loc[comparisons[0]] = males_SvP.reindex(columns=[3, 4, 8, 9, 5, 6]).values[0]
            holm_df.loc[comparisons[1]] = females_SvP.reindex(columns=[3, 4, 8, 9, 5, 6]).values[0]
            holm_df.loc[comparisons[2]] = sal_MvF.reindex(columns=[3, 4, 8, 9, 5, 6]).values[0]
            holm_df.loc[comparisons[3]] = psi_MvF.reindex(columns=[3, 4, 8, 9, 5, 6]).values[0]

        holm_df["p_holm"] = multipletests(
            holm_df.iloc[:, -1].astype(float), alpha=0.05, method='holm')[1]
        results = pd.concat([results, holm_df]).astype(float)

    results["significance"] = get_stars(results)
    results = results.fillna("")
    for line in results_print:
        print(line)
    print(results)

    return results, results_print


def stats_treatment_sex_third(data, variable, third_factor, third_factor_compare, third_paired):
    """
    Performs a mixed model analysis for the effect of treatment, sex,
    and a third factor + all interactions on a given variable.
    Animal ID is included as a random effect with an optional random slope for [third factor].

    Parameters:
    - data : pd.DataFrame containing the data to be analyzed.
    - variable : str. The name of the column in the DataFrame to be tested.
    - third_factor : str. The name of the third factor column in the DataFrame to be tested.

    Returns:
    - results : pd.DataFrame containing the results of the mixed model analysis.
    """
    results_print = []
    institutions = np.unique(data["institution"])
    if len(institutions) > 1:
        institution = "All labs"
    else:
        institution = institutions[0]

    analysis_name = f"------------- {institution}: {variable}~treatment*sex*{third_factor} -------------"
    results_print.append(f"\n{analysis_name}")

    data_y = {}
    third_factor_vals = np.unique(data[third_factor])
    for val in third_factor_vals:
        group_names = [f"m_sal_{val}", f"f_sal_{val}",
                       f"m_psi_{val}", f"f_psi_{val}"]
        (data_y[group_names[0]], data_y[group_names[1]],
         data_y[group_names[2]], data_y[group_names[3]]) = (
             get_group_data(data[data[third_factor] == val], variable))

    # first, print sample size and do levene's test for homogeneity of variances
    samples_statement, n_total = sample_size(data_y)
    results_print.append(samples_statement)
    levenes_p, levenes_statement = levenes_test(data_y, n_total)[1:]
    results_print.append(levenes_statement)
    parametric = levenes_p > 0.05

    vc_formula = None
    if (len(third_factor_vals) == 2) or (third_factor == 'cue' and variable == "retrieval"):
        data = data.copy()
        data[f"{third_factor}_num"] = (data[third_factor].map(
            {f"{third_factor_vals[0]}": -0.5,
             f"{third_factor_vals[1]}": 0.5})).astype(float)
        re_formula = f"0 + {third_factor}_num"
        formula = f"{variable} ~ C(treatment, Sum) * C(sex, Sum) * C({third_factor}, Sum)"
        stats_print = f"LME Model: {variable}~treatment*sex*{third_factor}+(0+{third_factor}|mouse)"

    else:
        re_formula = f"1 + {third_factor}"
        formula = f"{variable} ~ C(treatment, Sum) * C(sex, Sum) * {third_factor}"
        stats_print = f"LME Model: {variable}~treatment*sex*{third_factor}+(1+{third_factor}|mouse)"

    if institution == "All labs":
        if vc_formula is not None:
            vc_formula["institution"] = "0 + C(institution)"
        else:
            vc_formula = {"institution": "0 + C(institution)"}
        stats_print = stats_print + "+(1|institution)"

    data = organize_categories(data, lme=True)
    model = smf.mixedlm(formula,
                        data,
                        groups=data["mouse_ID"],
                        re_formula=re_formula,
                        vc_formula=vc_formula)
    result = model.fit()
    results_print.append(stats_print)
    results_print.append(f"llf (REML): {result.llf}. Converged: {result.converged}. ")
    results = result.summary().tables[1]
    results.replace('', np.nan, inplace=True)
    results = results[["Coef.", "Std.Err.", "[0.025", "0.975]", "z", "P>|z|"]]
    results = rename_results_indices(results).astype(float)

    # Do necessary posthocs
    result_pvals = rename_results_indices(result.pvalues)
    all_labs = False
    if institution == "All labs":
        if float(results.loc["institution Var", "Coef."]) >= 0.1:
            all_labs = True
    holm_df = three_factor_posthocs(result, result_pvals, data, variable,
                                    third_factor, parametric, third_factor_compare, third_paired, all_labs=all_labs)

    if holm_df is not None:
        results = pd.concat([results, holm_df])
        if all_labs:
            results_print.append(
                "Posthocs: Contrasts of estimated marginal means + Holm's corrections")
        elif parametric:
            results_print.append(
                "Posthocs: T-test (Unpaired) or Paired T-tests as indicated + Holm's Corrections. ")
        else:
            results_print.append(
                "Posthocs: Welch's Test (Unpaired) or Paired T-tests as indicated + Holm's Corrections. ")

    results["significance"] = get_stars(results)
    results = results.fillna("")
    for line in results_print:
        print(line)
    print(results)

    return results, results_print


def stats_treatment_sex_stress(data, variable):
    """
     Performs statistical analysis for the effect of treatment, sex, and stress on a given variable.

     Parameters:
     - data : pd.DataFrame containing the data to be analyzed.
     - variable : str. The name of the column in the DataFrame to be tested.
     Returns:
     - results : pd.DataFrame containing the results of the statistical tests performed.
    """
    all_labs = False
    results_print = []

    institutions = np.unique(data["institution"])
    if len(institutions) > 1:
        institution = "All labs"
    else:
        institution = institutions[0]

    analysis_name = f"------------- {institution}: {variable}~treatment*sex*stress -------------"
    results_print.append(f"\n{analysis_name}")

    data = organize_categories(data)
    cols_to_keep = ["mouse_ID", "institution", "treatment", "sex", "stress"] + [variable]
    data = data[cols_to_keep]

    data_ctrl = data[data["stress"] == "Ctrl"]
    data_stress = data[data["stress"] == "Stress"]

    m_sal_ctrl, f_sal_ctrl, m_psi_ctrl, f_psi_ctrl = get_group_data(
        data_ctrl, variable)
    m_sal_stress, f_sal_stress, m_psi_stress, f_psi_stress = get_group_data(
        data_stress, variable)

    data_y = {"m_sal_ctrl": m_sal_ctrl, "f_sal_ctrl": f_sal_ctrl,
              "m_psi_ctrl": m_psi_ctrl, "f_psi_ctrl": f_psi_ctrl,
              "m_sal_stress": m_sal_stress, "f_sal_stress": f_sal_stress,
              "m_psi_stress": m_psi_stress, "f_psi_stress": f_psi_stress}

    # first, print sample size and do levene's test for homogeneity of variances
    samples_statement, n_total = sample_size(data_y)
    results_print.append(samples_statement)
    levenes_p, levenes_statement = levenes_test(data_y, n_total)[1:]
    results_print.append(levenes_statement)

    if levenes_p < 0.05:
        parametric = False
    else:
        parametric = True

    # Set up the model
    formula = variable + ' ~ C(treatment, Sum)*C(sex, Sum)*C(stress, Sum)'
    if institution != "All labs":
        # if levenes failed.
        if parametric:
            results_print.append("3-Way ANOVA (treatment*sex*stress). ")
        else:
            results_print.append(
                "3-Way Robust ANOVA (treatment*sex*stress) (Wald tests using HC3 covariance). ")

        results = sm_ANOVA(formula=formula, data=data, parametric=parametric)
        results = rename_results_indices(results)
        result_pvals = results[results.columns[-1]]
    else:
        data = organize_categories(data, lme=True)
        model = smf.mixedlm(formula,
                            data=data,
                            groups=data["institution"])

        result = model.fit(method="lbfgs", reml=True)
        results_print.append(f"LME model: {variable}~treatment*sex*stress+(1|institution). "
                             f"llf (REML): {result.llf}. Converged: {result.converged}. ")
        results = result.summary().tables[1]
        results.replace('', np.nan, inplace=True)
        results = results.iloc[:, :-2].astype(float)
        results = rename_results_indices(results)
        result_pvals = rename_results_indices(result.pvalues)
        if float(results.loc["Group Var", "Coef."]) >= 0.1:
            all_labs = True

    holm_df = three_factor_posthocs(None, result_pvals, data, variable,
                                    "stress", parametric,
                                    third_factor_compare=True, third_paired=False, all_labs=all_labs)

    if holm_df is not None:
        results = pd.concat([results, holm_df])
        if all_labs:
            results_print.append(
                "Posthocs: Contrasts of estimated marginal means + Holm's corrections")
        elif parametric:
            results_print.append(
                "Posthocs: T-test + Holm's Corrections. ")
        else:
            results_print.append(
                "Posthocs: Welch's Test + Holm's Corrections. ")

    results["significance"] = get_stars(results)
    results = results.fillna("")
    for line in results_print:
        print(line)
    print(results)

    return results, results_print


def paired_t(data_expt_long, variable, pairing_var, group_names=None, method='holm'):
    groups = ["S", "P"]
    results_print = []
    results = pd.DataFrame(columns=[
        "t_df", "t_statistic", "t_0.025", "t_0.975", "Hedges' g",
        "P>|t|", "p_holm"])

    institutions = np.unique(data_expt_long["institution"])
    if len(institutions) > 1:
        institution = "All labs"
    else:
        institution = institutions[0]

    analysis_name = f"------------- {institution}: {variable}, paired t-test -------------"
    results_print.append(f"{analysis_name}")

    if institution == "All labs":
        data_to_analyze = data_expt_long
    else:
        data_to_analyze = data_expt_long[data_expt_long["institution"] == institution]

    for group_no, group in enumerate(groups):
        if group_names is not None:
            group_name = group_names[group_no]
        else:
            group_name = group
        comparison_name = f"Day 1 vs. Day 5 [{group_name}]"
        data = data_to_analyze[data_to_analyze["treatment"] == group]
        pair_vals = np.sort(np.unique(data[pairing_var]))
        g1 = data[variable][data[pairing_var] == pair_vals[0]]  # Day 1
        g2 = data[variable][data[pairing_var] == pair_vals[1]]  # Day 5

        # get sample size
        samples_statement = sample_size({pair_vals[0]: g1, pair_vals[1]: g2})[0]
        samples_statement = f"{group_name} {samples_statement}"
        results_print.append(samples_statement)

        # Results
        results = students_T(g1, g2, results, comparison_name, parametric=True, paired=True)
    results["p_holm"] = multipletests(results["P>|t|"], method=method)[1]
    results["significance"] = get_stars(results)

    for line in results_print:
        print(line)
    print(results)
    return results, results_print


# %% Data organization

def get_group_data(data_final, col):
    """
    Divides the data into four groups based on
    treatment (Saline vs. Psilocybin) and sex (Male vs. Female).

    Parameters:
    - data_final : pd.DataFrame containing the data to be analyzed.
    - col : str. The name of the column in data_final that is the outcome variable.
    Returns:
    - m_sal : list of values for male saline group.
    - f_sal : list of values for female saline group.
    - m_psi : list of values for male psilocybin group.
    - f_psi : list of values for female psilocybin group.
    """
    m_sal = data_final[(data_final["sex"] == "M") & (
        data_final["treatment"] == "S")][col].tolist()
    f_sal = data_final[(data_final["sex"] == "F") & (
        data_final["treatment"] == "S")][col].tolist()
    m_psi = data_final[(data_final["sex"] == "M") & (
        data_final["treatment"] == "P")][col].tolist()
    f_psi = data_final[(data_final["sex"] == "F") & (
        data_final["treatment"] == "P")][col].tolist()

    return m_sal, f_sal, m_psi, f_psi


def organize_categories(data, lme=True, plot=False):
    """
    Organizes the 'treatment' and 'sex' columns in the DataFrame as categorical variables
    with specified order.
    Parameters:
    - data : pd.DataFrame containing the data to be organized.
    Returns:
    - data : pd.DataFrame with organized categorical columns.
    """
    treatment_order = ['P', 'S']
    sex_order = ['F', 'M']
    stress_order = ['Stress', 'Ctrl']
    zone_order = ['social', 'empty']

    if plot:
        treatment_order = ['S', 'P']
        sex_order = ['M', 'F']
        stress_order = ['Ctrl', 'Stress']
        zone_order = ['empty', 'social']

    data = data.copy()
    data["mouse_ID"] = data["mouse_ID"].astype(str).astype("category")
    data['treatment'] = data['treatment'].astype('category')
    data['treatment'] = data['treatment'].cat.reorder_categories(
        treatment_order, ordered=True
    )

    data['sex'] = data['sex'].astype('category')
    data['sex'] = data['sex'].cat.reorder_categories(
        sex_order, ordered=True
    )

    if "stress" in data.columns:
        data['stress'] = data['stress'].astype('category')
        data['stress'] = data['stress'].cat.reorder_categories(
            stress_order, ordered=True
        )

    if "zone" in data.columns:
        data['zone'] = data['zone'].astype('category')
        data['zone'] = data['zone'].cat.reorder_categories(
            zone_order, ordered=True
        )

    if "object" in data.columns:
        data['object'] = data['object'].astype('category')
        data['object'] = data['object'].cat.reorder_categories(
            ['2', '1'], ordered=True
        )

    if "day" in data.columns:
        data['day'] = data['day'].astype('category')
        data['day'] = data['day'].cat.reorder_categories(
            ['3', '1'], ordered=True
        )

    if "day_scpp" in data.columns:
        data['day_scpp'] = data['day_scpp'].astype('category')
        data['day_scpp'] = data['day_scpp'].cat.reorder_categories(
            ['5', '1'], ordered=True
        )

    if "cue" in data.columns and np.max(data["cue"].astype(int)) == 2:
        data['cue'] = data['cue'].astype('category')
        data['cue'] = data['cue'].cat.reorder_categories(
            ['2', '1'], ordered=True
        )
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

        if "treatment[T.P]" in index:
            index = index.replace("treatment[T.P]", "treatment")

        if "sex[T.F]" in index:
            index = index.replace("sex[T.F]", "sex")

        if "stress[T.stress]" in index:
            index = index.replace("stress[T.stress]", "stress")

        if "C(treatment, Sum)" in index:
            index = index.replace("C(treatment, Sum)", "treatment")
            if "[S.P]" in index:
                index = index.replace("[S.P]", "")

        if "C(sex, Sum)" in index:
            index = index.replace("C(sex, Sum)", "sex")
            if "[S.F]" in index:
                index = index.replace("[S.F]", "")

        if "C(stress, Sum)" in index:
            index = index.replace("C(stress, Sum)", "stress")
            if "[S.Stress]" in index:
                index = index.replace("[S.Stress]", "")

        if "C(zone, Sum)[S.social]" in index:
            index = index.replace("C(zone, Sum)[S.social]", "zone")

        if "C(institution, Sum)" in index:
            index = index.replace("C(institution, Sum)", "institution")

        if "C(object, Sum)[S.2]" in index:
            index = index.replace("C(object, Sum)[S.2]", "object")

        if "C(day, Sum)[S.3]" in index:
            index = index.replace("C(day, Sum)[S.3]", "day")

        if "C(day_scpp, Sum)[S.5]" in index:
            index = index.replace("C(day_scpp, Sum)[S.5]", "day_scpp")

        if "C(cue, Sum)[S.2]" in index:
            index = index.replace("C(cue, Sum)[S.2]", "cue")
        results = results.rename(index={results.index.tolist()[i]: index})

    return results


# %% Effect sizes
def get_zscores(data_expt, variable, stress=False):
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
        data_expt_lab = data_expt[data_expt["institution"] == lab]

        # initialize dataframe
        psi_data_lab = pd.DataFrame(
            columns=["institution", "treatment", "sex", "stress", variable])

        # get data
        if stress:
            data_expt_lab = data_expt_lab[data_expt_lab["stress"] == "Ctrl"]
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

        id_info = np.array([lab, "P", "M", "Ctrl"]*len(zscores_male) +
                           [lab, "P", "F", "Ctrl"]*len(zscores_female)).reshape(len(zscores_all), -1)
        psi_data_lab[variable] = zscores_all
        psi_data_lab[["institution", "treatment", "sex", "stress"]] = id_info
        psi_data.append(psi_data_lab)

        if stress:
            psi_data_lab = pd.DataFrame(
                columns=["institution", "treatment", "sex", "stress", variable])
            data_expt_lab = data_expt[data_expt["institution"] == lab]
            data_expt_lab = data_expt_lab[data_expt_lab["stress"] == "Stress"]
            sal_male_data, sal_fem_data, psi_male_data, psi_fem_data = get_group_data(data_expt_lab,
                                                                                      variable)
            zscores_male = (psi_male_data - sal_male_mean)/sal_male_std
            zscores_female = (psi_fem_data - sal_fem_mean)/sal_fem_std
            zscores_all = np.concat((zscores_male, zscores_female))
            id_info = np.array([lab, "P", "M", "Stress"]*len(zscores_male) +
                               [lab, "P", "F", "Stress"]*len(zscores_female)).reshape(len(zscores_all), -1)
            psi_data_lab[variable] = zscores_all
            psi_data_lab[["institution", "treatment", "sex", "stress"]] = id_info
            psi_data.append(psi_data_lab)

            psi_data_lab = pd.DataFrame(
                columns=["institution", "treatment", "sex", "stress", variable])
            zscores_male = (sal_male_data - sal_male_mean)/sal_male_std
            zscores_female = (sal_fem_data - sal_fem_mean)/sal_fem_std
            zscores_all = np.concat((zscores_male, zscores_female))
            id_info = np.array([lab, "S", "M", "Stress"]*len(zscores_male) +
                               [lab, "S", "F", "Stress"]*len(zscores_female)).reshape(len(zscores_all), -1)
            psi_data_lab[variable] = zscores_all
            psi_data_lab[["institution", "treatment", "sex", "stress"]] = id_info
            psi_data.append(psi_data_lab)

    psi_data = pd.concat(psi_data)

    return psi_data


def calculate_hedges_g(group2, group1):
    """Calculates Hedges' g effect size."""
    n1 = len(group1)
    n2 = len(group2)

    df = n1 + n2 - 2

    mean_diff = np.mean(group2) - np.mean(group1)
    s1 = np.std(group1, ddof=1)
    s2 = np.std(group2, ddof=1)
    pooled_std = np.sqrt(((n1 - 1) * s1**2 + (n2 - 1) * s2**2) / df)

    # Cohen's d
    cohens_d = mean_diff / pooled_std

    # Hedges' g correction factor (bias correction)
    # Using the correction formula: J(df) = 1 - (3 / (4*df - 1))
    correction = 1 - (3 / (4 * df - 1))

    hedges_g = cohens_d * correction

    return hedges_g
