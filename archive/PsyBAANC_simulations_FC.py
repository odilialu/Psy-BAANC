# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 12:41:45 2024

@author: olu
"""


import seaborn as sns
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats
from itertools import groupby
from operator import itemgetter

filePath = r"C:\Users\olu\Documents\Psy-BAANC\Data\Simulations\FCRetrieval-ZScore.xlsx"
nLabs = 5
labID = ["Stanford", "Berkeley 1", "Berkeley 2", "UCSF 1", "UCSF 2"]
nSims = 10000
sampleSize = 10
nGroups = 4
nRows = 2

plt.rcParams['figure.dpi'] = 600
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 8

df = pd.read_excel(filePath)

################ clean up the data #########################

# remove the strings (happens with excluded data)
for runs in range(10):
    for row in range(len(df)):
        for col in range(len(df.columns)):
            if type(df.iloc[row, col]) == str:
                df.iloc[row, col:] = df.iloc[row, col:].shift(-1) # if string, delete and shift left 
                
################## SEPARATE DATA INTO GROUPS (AKA COLS IN PRISM) ##############
# separate data into groups
nanCols = df.isnull().values.all(axis=0)
nanColsIdx = [idx for idx, x in enumerate(nanCols) if x]

groups=[]
for k, g in groupby(enumerate(nanColsIdx), lambda i_x: i_x[0] - i_x[1]):
    groups.append(list(map(itemgetter(1), g)))
    
groupCutoff = []
for i in range(len(groups)):
    if len(groups[i]) > 1:
        groupCutoff.append(groups[i][0])
groupCutoff.insert(0, 0)
groupCutoff.insert(nGroups, len(df.columns))
        
groupsData=[]
for i in range(nGroups):
    groupsData.append(df.iloc[:, groupCutoff[i]:groupCutoff[i+1]])


########################## 1. Format data  ############################

def format_nonans_factors(data, drug, sex, drugsex):
    data_formatted = []
    for row in range(nRows):
        tempList = []
        for col in range(len(data.columns)):
            if pd.isna(data.iloc[row, col]) == 0: # get rid of nans
                tempList.append(data.iloc[row, col])
        drugString = [drug]*len(tempList) # add factors to each data point
        sexString = [sex]*len(tempList)
        drugsexString = [drugsex]*len(tempList)
        
        data_formatted_array = np.array((tempList, sexString, drugString, drugsexString)).T
                
        data_formatted.append(data_formatted_array)

    return data_formatted

salMalesFormatted = format_nonans_factors(groupsData[0], 'sal', 'males', 'sal males')
salFemalesFormatted = format_nonans_factors(groupsData[2], 'sal', 'females', 'sal females')
psiMalesFormatted = format_nonans_factors(groupsData[1], 'psi', 'males', 'psi males')
psiFemalesFormatted = format_nonans_factors(groupsData[3], 'psi', 'females', 'psi females')


## create KDE plot 
pooledDataFormatted = []
for i in range(nRows):
    tempDF = np.concatenate((salMalesFormatted[i], psiMalesFormatted[i], salFemalesFormatted[i], psiFemalesFormatted[i]))
    tempDF = pd.DataFrame(data=tempDF, columns=['zScore', 'sex', 'drug', 'drugsex'])
    tempDF.zScore = tempDF.zScore.astype(float)
    
    pooledDataFormatted.append(tempDF)

    plt.figure(figsize=(1.4, 1.4))
    colors = ["#808080", "#2E8B57", "#C0C0C0","#50C878"]
    sns.set_palette(colors)
    ax = sns.kdeplot(tempDF, x="zScore", hue="drugsex", common_norm=False, legend=False)
    plt.xlabel("z-score")
    plt.gca().set_position([0, 0, 1, 1])
    plt.savefig("C:/Users/olu/Documents/Psy-BAANC/Data/figures/retrievalKDE.svg")
    plt.show()


########################## 2. Resample from KDE distribution 10,000 times 
# get KDE distributions
kdeSalMales = []
kdeSalFemales = []
kdePsiMales = []
kdePsiFemales = []
for i in range(nRows):
    kdeSalMales.append(stats.gaussian_kde(salMalesFormatted[i][:, 0].astype(np.float64)))
    kdeSalFemales.append(stats.gaussian_kde(salFemalesFormatted[i][:, 0].astype(np.float64)))
    kdePsiMales.append(stats.gaussian_kde(psiMalesFormatted[i][:, 0].astype(np.float64)))
    kdePsiFemales.append(stats.gaussian_kde(psiFemalesFormatted[i][:, 0].astype(np.float64)))

# allocate factors
def randomSample(kdeDrugSex, sampleSize, drug, sex):
    sampleArray = []
    for i in range(nRows): 
        sampleDrugSex = kdeDrugSex[i].resample(sampleSize).T
        drugString = np.array([drug]*len(sampleDrugSex))
        sexString = np.array([sex]*len(sampleDrugSex))
        timeString = np.array([i+1]*len(sampleDrugSex))
        stringArray = np.column_stack((drugString, sexString, timeString))
        sampleArray.append(np.concatenate((sampleDrugSex, stringArray), axis=1))
        
    sample = np.concatenate((sampleArray[0], sampleArray[1]))
            
    return sample

def sm_ANOVA(formula, data):

    # create model object. NOTE: Model hasn't been fit yet, so there are no results yet!
    ols_lm = smf.ols(formula, data=data)

    # call this to fit the linear model
    fit = ols_lm.fit()

    # Grab the ANOVA results from the linear model
    table = sm.stats.anova_lm(fit, typ=2)
    return table
    
resultsP = np.empty([nSims, 2])
resultsP[:] = np.nan
for run in range(nSims):
    sampleSalMales = randomSample(kdeSalMales, sampleSize, "sal", "males")
    sampleSalFemales = randomSample(kdeSalFemales, sampleSize, "sal", "females")
    samplePsiMales= randomSample(kdePsiMales, sampleSize, "psi", "males")
    samplePsiFemales = randomSample(kdePsiFemales, sampleSize, "psi", "females")
    sampleAll = np.concatenate((sampleSalMales, sampleSalFemales, samplePsiMales, samplePsiFemales))
    sampleDF = pd.DataFrame(data=sampleAll, columns=["y", "drug", "sex", "cue"])
    sampleDF.y = sampleDF.y.astype(float)
    
    results = sm_ANOVA(formula='y ~ drug*sex*cue', data=sampleDF)
    resultsP[run, 0] = results["PR(>F)"][0]
    resultsP[run, 1] = results["PR(>F)"][3]
    

plt.figure(figsize=(1.5,1.5))
ax = sns.histplot(data=resultsP[:,0], bins=10, alpha=0.3, color='r', stat="probability", label='main drug effect')
# sns.histplot(data=resultsP[:, 1], bins=40, alpha=0.3, color='b', stat="probability", label='interaction effect')
# xticks = np.linspace(0, 1, 11)
plt.xlabel('p-value')
# plt.legend(loc='upper right')

plt.xlim([0, 0.1])
plt.ylim([0, 1])
ax.set_xticks((0, 0.05, 0.1))
ax.set_yticks((0, 0.2, 0.4, 0.6, 0.8, 1))
plt.axvline(x=0.05, ymin=0, ymax=0.95, color='r', linestyle=':')
# ax.set_xticks(xticks)

mainEffectP = round(np.mean(resultsP[:, 0]<0.05)*100)

mainEffectPString = str(mainEffectP) + "% p<α"

# interactEffectP = round(np.mean(resultsP[:, 1]<0.05)*100)
# interactEffectPString = str(interactEffectP) + "% p<α"

# plt.text(0.15, 0.94, mainEffectPString, va='top', transform=ax.transAxes)
# plt.text(0.2, 0.84, interactEffectPString, va='top', transform=ax.transAxes)

plt.savefig("C:/Users/olu/Documents/Psy-BAANC/Data/figures/retrievalPValues.svg", dpi=300)
plt.show()
# sns.histplot(data=resultsP[:, 1], bins=50, alpha=0.3, color='b', stat="probability")
