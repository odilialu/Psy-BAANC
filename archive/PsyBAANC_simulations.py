# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 18:08:04 2024
Psy-BAANC simulations 
@author: olu
"""

import seaborn as sns
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats

filePath = r"C:/Users/olu/Documents/Psy-BAANC/Data/Simulations/noeAcute-ZScore.xlsx"
nLabs = 5
labID = ["Stanford", "Berkeley 1", "Berkeley 2", "UCSF 1", "UCSF 2"]
nSims = 10000
sampleSize = 10

df = pd.read_excel(filePath)

############# remove the strings
for runs in range(10):
    for row in range(len(df)):
        for col in range(len(df.columns)):
            if type(df.iloc[row, col]) == str:
                df.iloc[row, col:] = df.iloc[row, col:].shift(-1)
                
####################### separate data into sal and psi 
nanCols = df.isnull().values.all(axis=0)

for idx in range(len(nanCols)):
    if sum(nanCols[idx:idx+5]) == 5:
        salPsiColCutoff = idx
        break

salData = df.iloc[:, 0:salPsiColCutoff]
psiData = df.iloc[:, salPsiColCutoff:]


################ 1. Format data into groups 

def formatData(data, row, sex, drug, drugSex):
    pooledData = []
    for col in range(len(data.columns)):
        if pd.isna(data.iloc[row, col]) == 0:
            pooledData.append(data.iloc[row, col])
            
    drugString = [drug]*len(pooledData)
    sexString = [sex]*len(pooledData)
    drugSexString = [drugSex]*len(pooledData)
    pooledData = np.array((pooledData, sexString, drugString, drugSexString)).T

    return pooledData

pooledDataSalMales = formatData(salData, 0, "males", "sal", "sal males")
pooledDataSalFemales = formatData(salData, 1, "females", "sal", "sal females")
pooledDataPsiMales = formatData(psiData, 0, "males", "psi", "psi males")
pooledDataPsiFemales = formatData(psiData, 1, "females", "psi", "psi females")

kdeSalMales = stats.gaussian_kde(pooledDataSalMales[:, 0].astype(np.float64))
kdeSalFemales = stats.gaussian_kde(pooledDataSalFemales[:, 0].astype(np.float64))
kdePsiMales = stats.gaussian_kde(pooledDataPsiMales[:, 0].astype(np.float64))
kdePsiFemales = stats.gaussian_kde(pooledDataPsiFemales[:, 0].astype(np.float64))

pooledData = np.concatenate((pooledDataSalMales, pooledDataPsiMales, pooledDataSalFemales, pooledDataPsiFemales))
pooledDataDF = pd.DataFrame(data=pooledData, columns=["zScore", "sex", "drug", "drugSex"])
pooledDataDF.zScore = pooledDataDF.zScore.astype(float)


plt.rcParams['figure.dpi'] = 600
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 8

plt.figure(figsize=(1.2, 1.2))
colors = ["#808080", "#2E8B57", "#C0C0C0","#50C878"]
sns.set_palette(colors)
ax =sns.kdeplot(pooledDataDF, x="zScore", hue="drugSex", common_norm=False, legend=False)
plt.gca().set_position([0, 0, 1, 1])
plt.xticks([-5, 0, 5])
plt.yticks([0, 0.2, 0.4, 0.6])
plt.ylim([0, 0.6])
plt.xlim([-5, 5])
plt.xlabel('z-score')
figurePath = filePath.replace("xlsx", "svg")
plt.savefig(figurePath)

plt.show()


# 2. Resample 10,000 times 

def allocateFactors(kdeDrugSex, sampleSize, drug, sex):
    sampleDrugSex = kdeDrugSex.resample(sampleSize).T
    # sampleDrugSex = random.sample(pooledDataDrugSex, n)
    drugString = np.array([drug]*len(sampleDrugSex))
    sexString = np.array([sex]*len(sampleDrugSex))
    stringArray = np.column_stack((drugString, sexString))
    sampleArray = np.concatenate((sampleDrugSex, stringArray), axis=1)
    
    return sampleArray

def sm_ANOVA(formula, data):
    """ Uses the statsmodels api to perform an ANOVA. From the previous notebook """

    # create model object. NOTE: Model hasn't been fit yet, so there are no results yet!
    ols_lm = smf.ols(formula, data=data)

    # call this to fit the linear model
    fit = ols_lm.fit()

    # Grab the ANOVA results from the linear model
    table = sm.stats.anova_lm(fit, typ=2)
    return table
    
resultsP = np.empty([nSims, 3])
resultsP[:] = np.nan
for run in range(nSims):
    sampleArraySalMales = allocateFactors(kdeSalMales, sampleSize, "sal", "males")
    sampleArraySalFemales = allocateFactors(kdeSalFemales, sampleSize, "sal", "females")
    sampleArrayPsiMales= allocateFactors(kdePsiMales, sampleSize, "psi", "males")
    sampleArrayPsiFemales = allocateFactors(kdePsiFemales, sampleSize, "psi", "females")
    sampleArray = np.concatenate((sampleArraySalMales, sampleArraySalFemales, sampleArrayPsiMales, sampleArrayPsiFemales))
    sampleDF = pd.DataFrame(data=sampleArray, columns=["y", "drug", "sex"])
    sampleDF.y = sampleDF.y.astype(float)
    
    results = sm_ANOVA(formula='y ~ drug*sex', data=sampleDF)
    resultsP[run, 0] = results["PR(>F)"][0]
    resultsP[run, 1] = results["PR(>F)"][1]
    resultsP[run, 2] = results["PR(>F)"][2]
    

plt.figure(figsize=(1.5, 1.5))
ax = sns.histplot(data=resultsP[:,0], bins=20, alpha=0.3, color='r', stat="probability", label='main drug effect', legend=False)
sns.histplot(data=resultsP[:, 2], bins=20, alpha=0.3, color='b', stat="probability", label='interaction effect', legend=False)
plt.xlabel('p-value')
# plt.legend(loc='upper right')
plt.axvline(x=0.05, ymin=0, ymax=0.95, color='r', linestyle=':')
plt.xlim([0, 1])
plt.ylim([0, 0.6])
ax.set_xticks((0, 0.5, 1))
ax.set_yticks((0, 0.2, 0.4, 0.6))
mainEffectP = round(np.mean(resultsP[:, 0]<0.05)*100)
mainEffectPString = str(mainEffectP) + "% p<α"

interactEffectP = round(np.mean(resultsP[:, 2]<0.05)*100)
interactEffectPString = str(interactEffectP) + "% p<α"

# plt.text(0.45, 0.97, mainEffectPString, va='top', transform=ax.transAxes)
# plt.text(0.45, 0.89, interactEffectPString, va='top', transform=ax.transAxes)
figurePath = filePath.replace(".xlsx", "-simResult.svg")
plt.savefig(figurePath)
plt.show()
# sns.histplot(data=resultsP[:, 1], bins=50, alpha=0.3, color='b', stat="probability")
print(mainEffectPString)
