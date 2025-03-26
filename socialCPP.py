# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 09:08:13 2024

social CPP analysis, based on body coordinate and line dividing two chambers.  

@author: olu
"""

import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
import os
import cv2 as cv
from IPython import get_ipython


directory=r'Y:/PsyBAANC/paperExperiments/social CPP/SkinnerBoxExperiments/Adult salpsi experiments/fourthCohort/preTest/reanalyze/leftSocial'
vidType = '.mp4'
csvType = 'csv'
bodyX = 1
framerate = 30
totalTime = 30 # in min




def getFilePaths(directory, stringEnd):
    filePaths = []
    fileNames = os.listdir(directory)
    for fileNo in range(len(fileNames)):
        if fileNames[fileNo].endswith(stringEnd):
            joined = os.path.join(directory, fileNames[fileNo])
            filePaths.append(joined)
            
    return filePaths

csvPaths = getFilePaths(directory, csvType)
vidPaths = getFilePaths(directory, vidType)
         


##### ANALYSIS ##########################################
leftTime = np.empty(len(vidPaths))
leftTime[:] = np.nan
rightTime = np.empty(len(vidPaths))
rightTime[:] = np.nan
ambiguousTime = np.empty(len(vidPaths))
ambiguousTime[:] = np.nan

leftTimePrct = np.empty(len(vidPaths))
leftTimePrct[:] = np.nan
rightTimePrct = np.empty(len(vidPaths))
rightTimePrct[:] = np.nan



for i in range(len(vidPaths)): 
    
    ## define social and isolation zones
    cap = cv.VideoCapture(vidPaths[i])
    totalTime = cap.get(cv.CAP_PROP_FRAME_COUNT)/framerate
    cap.set(cv.CAP_PROP_POS_MSEC, 100)
    success, image = cap.read()
    get_ipython().run_line_magic('matplotlib', 'qt')
    plt.figure()
    plt.imshow(image)
    [lineTop] =plt.ginput(1)
    [lineBottom] = plt.ginput(1)
    plt.close()
        
    ## read in body data 
    df = pd.read_csv(csvPaths[i])
    dfTrunc = df.iloc[2:54002, bodyX:bodyX+2].astype(float)
    
    ## find cross product 
    vX = lineBottom[0] - lineTop[0]
    vY = lineBottom[1] - lineTop[1]
    
    yCoordMinusLineTop = dfTrunc.iloc[:, 1] - lineTop[1]
    xCoordMinusLineTop = dfTrunc.iloc[:, 0] - lineTop[0]
    
    crossProduct = (vX*yCoordMinusLineTop) - (vY*xCoordMinusLineTop)
    crossProductLogical = crossProduct>0
    
    #### if cross product greater than 0, then point is to the left of a line. 
    leftTime[i] = sum(crossProduct>0)/framerate
    rightTime[i] = sum(crossProduct<0)/framerate
    ambiguousTime[i] = sum(crossProduct==0)/framerate
    
    leftTime[i] = leftTime[i] + (ambiguousTime[i]/2)
    rightTime[i] = rightTime[i] + (ambiguousTime[i]/2)
    
    # leftTimePrct[i] = leftTime[i]/totalTime*100
    # rightTimePrct[i] = rightTime[i]/totalTime*100
    
    leftTimePrct[i] = leftTime[i]/30/60*100
    rightTimePrct[i] = rightTime[i]/30/60*100
    

    ## plot body coordinate track
    get_ipython().run_line_magic('matplotlib', 'inline')

    col=[]
    for i in range(len(crossProductLogical)):
        if crossProductLogical.iloc[i] == True:
            col.append('r')
        else:
            col.append('b')
            
            
    plt.rcParams['figure.dpi'] = 600
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 8
    
    plt.figure()
    plt.imshow(image)
    plt.scatter(dfTrunc.iloc[:, 0], dfTrunc.iloc[:, 1], c=col, s=0.2, alpha=0.5)
    plt.show()
    
    plt.figure()
    plt.imshow(image)
    plt.plot(dfTrunc.iloc[:, 0], dfTrunc.iloc[:, 1], 'k', alpha=0.5, linewidth=0.2)
    plt.show()
    


prefScore = leftTime/rightTime
prefScore2 = rightTime/leftTime
