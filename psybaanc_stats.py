# -*- coding: utf-8 -*-
"""
Created on Tue Jul  1 15:40:52 2025

@author: olu
Functions used for psy-baanc statistical analyses 
"""
#%%
# Import packages
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from scipy import stats
import scikit_posthocs as sp 

#%%
def levenes_test_dataframe(data_final):
    """
    Performs Levene's test for homogeneity of variance across groups in a DataFrame.
    
    Parameters:
    - data_final : pd.DataFrame containing the data to be tested. The first three columns are identity columns and not data columns.
        
    Returns:
    - levenes_results : pd.DataFrame containing the results of Levene's test for each measure.'
    """
    statistic = [None]*len(data_final.columns[3:])
    pvalue = [None]*len(data_final.columns[3:])
    statement = [None]*len(data_final.columns[3:])
    column_names = data_final.columns[3:].tolist()
    col_index = 0
    for col in data_final.columns[3:]:
        m_sal_col = data_final[(data_final["Sex"] == "M") & (data_final["Treatment"] == "S")][col].tolist()
        f_sal_col = data_final[(data_final["Sex"] == "F") & (data_final["Treatment"] == "S")][col].tolist()
        m_psi_col = data_final[(data_final["Sex"] == "M") & (data_final["Treatment"] == "P")][col].tolist()
        f_psi_col = data_final[(data_final["Sex"] == "F") & (data_final["Treatment"] == "P")][col].tolist()
        
        statistic[col_index], pvalue[col_index] = stats.levene(m_sal_col, f_sal_col, m_psi_col, f_psi_col)
        
        groups = 4
        df1 = str(groups-1)
        df2 = str(len(data_final) - groups)
        statement[col_index] = col + "; Levene's test: F(" + df1 + ", " + df2 + ") = " + str(round(statistic[col_index], 3)) + ", p = " + str(round(pvalue[col_index], 3))
        print(statement[col_index])
        col_index = col_index+1
    
    levenes_results = pd.DataFrame({
        "Measure": column_names,
        "Statistic": statistic,
        "P_value": pvalue,
        "Statement": statement
        })    
    return levenes_results
    
#%%
def sm_ANOVA(formula, data):
    """
    Performs a two-way ANOVA using statsmodels.
    
    Parameters
    - formula : str formula for the ANOVA model
    - data : pd.DataFrame containing the data to be analyzed.
    
    Returns:
    - table : pd.DataFrame containing the ANOVA results, including sum of squares, degrees of freedom, F-statistic, and p-value.
    """

    # create model object. NOTE: Model hasn't been fit yet, so there are no results yet!
    ols_lm = smf.ols(formula, data=data)

    # call this to fit the linear model
    fit = ols_lm.fit()

    # Grab the ANOVA results from the linear model
    table = sm.stats.anova_lm(fit, typ=2)
    return table

#%%
def kruskal_wallis(data_final, col):
    """    
    Performs the Kruskal-Wallis H-test for independent samples.
    
    Parameters:
    - data_final : pd.DataFrame containing the data to be tested.
    col : str
        The name of the column in the DataFrame to be tested.

    Returns:
    - h_statistic : float. Kruskal-Wallis H statistic.
    - p_value : float. The p-value associated with the H statistic.
    """
    m_sal_col = data_final[(data_final["Sex"] == "M") & (data_final["Treatment"] == "S")][col].tolist()
    f_sal_col = data_final[(data_final["Sex"] == "F") & (data_final["Treatment"] == "S")][col].tolist()
    m_psi_col = data_final[(data_final["Sex"] == "M") & (data_final["Treatment"] == "P")][col].tolist()
    f_psi_col = data_final[(data_final["Sex"] == "F") & (data_final["Treatment"] == "P")][col].tolist()
    
    h_statistic, p_value = stats.kruskal(m_sal_col, f_sal_col, m_psi_col, f_psi_col)

    return h_statistic, p_value

#%%
def dunns(data_final, col):
    """
    Performs Dunn's post-hoc test after a Kruskal-Wallis test.
    
    Parameters:
    - data_final : pd.DataFrame containing the data to be tested. The first three columns are identity columns and not data columns.
    - col : str. The name of the column in the DataFrame to be tested.
    
    Returns
    - pvals_holm : list of p-values adjusted using the Holm method for multiple comparisons.   
    """

    m_sal_col = data_final[(data_final["Sex"] == "M") & (data_final["Treatment"] == "S")][col].tolist()
    f_sal_col = data_final[(data_final["Sex"] == "F") & (data_final["Treatment"] == "S")][col].tolist()
    m_psi_col = data_final[(data_final["Sex"] == "M") & (data_final["Treatment"] == "P")][col].tolist()
    f_psi_col = data_final[(data_final["Sex"] == "F") & (data_final["Treatment"] == "P")][col].tolist()
    
    dunn_results = sp.posthoc_dunn([m_sal_col, f_sal_col, m_psi_col, f_psi_col], p_adjust=None)
    comparisons = [(1, 3), (2, 4), (1, 2), (3, 4)]
    pvals = [dunn_results.loc[g1, g2] for g1, g2 in comparisons]
    rejected, pvals_holm, _, _ = multipletests(pvals, method='holm')
    
    return pvals_holm

#%%
def sidaks(data_final, col):
    """
    Performs Sidak's post-hoc test after a t-test.

    Parameters:
    - data_final : pd.DataFrame containing the data to be tested. The first three columns are identity columns and not data columns.
    - col : str. The name of the column in the DataFrame to be tested.

    Returns:
    - results : pd.DataFrame containing the t-statistics, p-values, and Sidak-adjusted p-values for each comparison.
    """

    m_sal_col = data_final[(data_final["Sex"] == "M") & (data_final["Treatment"] == "S")][col].tolist()
    f_sal_col = data_final[(data_final["Sex"] == "F") & (data_final["Treatment"] == "S")][col].tolist()
    m_psi_col = data_final[(data_final["Sex"] == "M") & (data_final["Treatment"] == "P")][col].tolist()
    f_psi_col = data_final[(data_final["Sex"] == "F") & (data_final["Treatment"] == "P")][col].tolist()
    
    sub_data = [m_sal_col, f_sal_col, m_psi_col, f_psi_col]
    comparisons = [(0, 2), (1, 3), (0, 1), (2, 3)]
    comparison_names = ["m,sal v. m,psi", "f,sal v. f,psi", "m,sal v. f,sal", "m,psi v. f,psi"]
    t_stat = np.empty(len(comparisons))
    p_val = np.empty(len(comparisons))
    for comparison_no, comparison in enumerate(comparisons):
        g1 = sub_data[comparison[0]]
        g2 = sub_data[comparison[1]]
        t_stat[comparison_no], p_val[comparison_no] = stats.ttest_ind(g1, g2)
        
    _, pvals_sidak, _, _ = multipletests(p_val, method='sidak')
    
    results_dict = {
        "t_stat": t_stat,
        "p_val": p_val,
        "p_val_sidak": pvals_sidak}
    results = pd.DataFrame(data=results_dict, index=comparison_names)
    
    return results

#%% 
def sample_size(data_final):
    """
    Calculates the sample size for each group in the DataFrame.
    
    Parameters:
    - data_final : pd.DataFrame containing the data to be analyzed. The first three columns are identity columns and not data columns.
    
    Returns:
    - statement : A string summarizing the sample sizes for each group.
    """
    
    m_sal_n = str(len(data_final[(data_final["Sex"] == "M") & (data_final["Treatment"] == "S")]))
    f_sal_n = str(len(data_final[(data_final["Sex"] == "F") & (data_final["Treatment"] == "S")]))
    m_psi_n = str(len(data_final[(data_final["Sex"] == "M") & (data_final["Treatment"] == "P")]))
    f_psi_n = str(len(data_final[(data_final["Sex"] == "F") & (data_final["Treatment"] == "P")]))
    
    statement = "Sample size: m, sal = " + m_sal_n + "; f, sal = " + f_sal_n + "; m, psi = " + m_psi_n + "; f, psi = " + f_psi_n

    return statement