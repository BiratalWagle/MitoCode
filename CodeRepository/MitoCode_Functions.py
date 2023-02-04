# MitoCode_Functions

import pandas as pd
import numpy as np
import os

def Dist(pointOne, pointTwo): # get the distance between two points
    a = 0 
    for i in range(len(pointOne)):  # this allows us to consider it for an n dimensional array
        a += (pointOne[i] - pointTwo[i])**2
    return np.sqrt(a)

def GetData(k): 
    dir = k
    for file in os.listdir(dir):
        if file.endswith(".txt"):
            skelePos = pd.read_table(os.path.join(dir, file))
        if file.endswith(".gnet"):
            nodeDist = pd.read_table(os.path.join(dir, file)) # a list of the nodes and the distance b/w them
        if file.endswith(".coo"): 
            nodeList = pd.read_table(os.path.join(dir, file)) # a list of the nodes and their coordinates
        if file.endswith(".cc"): 
            compList = pd.read_table(os.path.join(dir, file)) # a list of the nodes and their coordinates
    return nodeDist, nodeList, skelePos, compList

def findDistofNodes(node1, node2, temp):  # temp is a clean node distance thing
    t1 = temp.loc[(temp.X == node1) & (temp.Y == node2)]
    if len(t1) < 1:
        t1 = temp.loc[(temp.X == node2) & (temp.Y  == node1)]
    if len(t1) < 1:
        print("Node's are not connected directly, so no distance between them") 
        return None
    return t1.dist.values[0]
    
def findLength(test, temp): 
    # test = bigDf1.loc[bigDf1.cc == k]
    linids = []
    lens = []
    endparts = []

    for i in set(test.line_id): 
        secTest = test.loc[test.line_id == i]
        b = secTest.loc[secTest.nodeState == True]
        if len(b) < 2: 
            endparts.append(secTest)
        else: 
            linids.append(b)
        
    for i in linids: 
        nodesInid = i.node.values
        lens.append(findDistofNodes(nodesInid[0], nodesInid[1], temp))

    for j in endparts: 
        lens.append(findLength2(j))
    
    return np.sum(lens)

def findAvgWidth(df): 
    widths = df['width_(um)'].values
    return np.mean(widths)

def findNumberofNodes(df): 
    return len(df.loc[df.nodeState == True])

def findPixelIntensity(df): 
    pixint = df.pix_ratio.mean()
    return pixint

def findVolume2(df): 
    vals = df[['x', 'y', 'z']].values
    widths = df['width_(um)'].values
    volume = []
    for i in range(len(vals) - 1):
        diff_of_vals = vals[i + 1] - vals[i]
        avg_width = (widths[i + 1] + widths[i])/2
        volume.append(np.sqrt((diff_of_vals * avg_width)**2))
    return np.sum(volume)

def findLength2(df): # find's length of the entire thing
    X = set(df.line_id.values)
    bigDisList = []
    for i in X: 
        bigDisList = []
        tf = df.loc[df.line_id == i]
        vals = df[['x', 'y', 'z']].values
        disList = []
        for i in range(len(vals) - 1): 
            disList.append(Dist(vals[i], vals[i + 1]))
        bigDisList.append(np.sum(disList))
    return np.sum(bigDisList)