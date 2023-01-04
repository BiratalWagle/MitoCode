# MitoCode_Functions

import pandas as pd
import numpy as np
import os

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

    # nodeDist = pd.read_table(str(k) + ".gnet") # a list of the nodes and the distance b/w them
    # nodeList = pd.read_table(str(k) + ".coo") # a list of the nodes and their coordinates
    # # mitoAttList = pd.read_table(str(k) + ".mitograph") # shows the mitochondrial attributes
    # skelePos = pd.read_table(str(k) + ".txt") # a list of the skeleton coordinates
                                        # a list of the values of the pixel intensities
    return nodeDist, nodeList, skelePos, compList

def GetVolume(table): 
    X, Y, Z = table['x'], table['y'], table['z']
    width = table['width_(um)']
    VolDiff = []

    X = np.array(X)
    Y = np.array(Y)
    Z = np.array(Z)
    width = np.array(width)

    for i in range(len(X)):
        try:
            VolDiff.append(np.sqrt((X[i+1]-X[i])**2 + (Y[i+1]-Y[i])**2 + (Z[i+1]-Z[i])**2) * width[i] )
        except IndexError:
            VolDiff.append(0)
    
    return VolDiff
    
