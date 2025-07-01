# -*- coding: utf-8 -*-
"""
Created on Tue Jul  1 15:40:52 2025

@author: olu
"""
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from scipy import stats
import scikit_posthocs as sp 

#%%
def levenes_test_dataframe(data_final):
    
    statistic = [None]*len(data_final.columns[3:])
    pvalue = [None]*len(data_final.columns[3:])
    column_names = data_final.columns[3:].tolist()
    col_index = 0
    for col in data_final.columns[3:]:
        m_sal_col = data_final[(data_final["Sex"] == "M") & (data_final["Treatment"] == "S")][col].tolist()
        f_sal_col = data_final[(data_final["Sex"] == "F") & (data_final["Treatment"] == "S")][col].tolist()
        m_psi_col = data_final[(data_final["Sex"] == "M") & (data_final["Treatment"] == "P")][col].tolist()
        f_psi_col = data_final[(data_final["Sex"] == "F") & (data_final["Treatment"] == "P")][col].tolist()
        
        statistic[col_index], pvalue[col_index] = stats.levene(m_sal_col, f_sal_col, m_psi_col, f_psi_col)
        col_index = col_index+1
        
    levenes_results = pd.DataFrame({
        "Measure": column_names,
        "Statistic": statistic,
        "P_value": pvalue
        })    
    return levenes_results
    
#%%
def sm_ANOVA(formula, data):
    """ Uses the statsmodels api to perform an ANOVA. From the previous notebook """

    # create model object. NOTE: Model hasn't been fit yet, so there are no results yet!
    ols_lm = smf.ols(formula, data=data)

    # call this to fit the linear model
    fit = ols_lm.fit()

    # Grab the ANOVA results from the linear model
    table = sm.stats.anova_lm(fit, typ=2)
    return table

#%%
def kruskal_wallis(data_final, col):
    m_sal_col = data_final[(data_final["Sex"] == "M") & (data_final["Treatment"] == "S")][col].tolist()
    f_sal_col = data_final[(data_final["Sex"] == "F") & (data_final["Treatment"] == "S")][col].tolist()
    m_psi_col = data_final[(data_final["Sex"] == "M") & (data_final["Treatment"] == "P")][col].tolist()
    f_psi_col = data_final[(data_final["Sex"] == "F") & (data_final["Treatment"] == "P")][col].tolist()
    
    h_statistic, p_value = stats.kruskal(m_sal_col, f_sal_col, m_psi_col, f_psi_col)

    return h_statistic, p_value

#%%
def dunns(data_final, col):
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