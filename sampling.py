# -*- coding: utf-8 -*-

import json
import numpy as np
import matplotlib.pyplot as plt

# read files
with open('bigBenResults_1.json', 'r') as myBigBenfile:
    bigBenData = myBigBenfile.read()
    
with open('uttResults_1.json', 'r') as myUttfile:
    uttData = myUttfile.read()

# parse files
bigBenObj = json.loads(bigBenData)
uttObj = json.loads(uttData)

# create data list objects
bigBenDataList = []

for line in bigBenObj:
    data = []
    data.append(float("{0:.2f}".format(line["TIMESTAMP"])))
    data.append(line["topLeftKg"])
    data.append(line["topRightKg"])
    data.append(line["bottomLeftKg"])
    data.append(line["bottomRightKg"])
    data.append(line["weightKg"]-bigBenObj[0]["weightKg"])
    data.append(line["gravity"]["X"])
    data.append(line["gravity"]["Y"])
    bigBenDataList.append(data)

    
uttAcquisitions = uttObj["acquisitions"]
uttDataList = []
g_const = 9.8

for line in uttAcquisitions:
    data = []
    data.append(line["time_stamp"])
    data.append((line["pressure"][0]-uttAcquisitions[0]["pressure"][0])/g_const)
    data.append((line["pressure"][1]-uttAcquisitions[0]["pressure"][1])/g_const)
    data.append((line["pressure"][3]-uttAcquisitions[0]["pressure"][3])/g_const)
    data.append((line["pressure"][2]-uttAcquisitions[0]["pressure"][2])/g_const)
    uttDataList.append(data)
    
bigBenMat = np.array(bigBenDataList)
uttMat = np.array(uttDataList)

somme = np.sum(uttMat[::, 1:5],axis=1)
uttMat=np.append(uttMat,somme[:,None],axis=1)

bigBenStartTime = np.where(bigBenMat[:,5]>0.25)
uttStartTime = np.where(uttMat[:,5]>0.25)

uttMat = uttMat[uttStartTime,:]
bigBenMat = bigBenMat[bigBenStartTime,:]

timestamp_array = bigBenMat[0][:, 0]
unique_keys, indices = np.unique(timestamp_array, return_index=True)
bigBenMat = bigBenMat[:, np.sort(indices)]

delta_time = bigBenMat[0][0][0] - uttMat[0][0][0]
bigBenMat[0][:, 0] = bigBenMat[0][:, 0]-delta_time

def getInterpolation(firstLine, lastLine):
    steps = int(round((lastLine[0]-firstLine[0])/0.01));
    interpolations = []
    for i in range(1,steps):
        newLine = firstLine + (((lastLine-firstLine)/steps)*i)
        interpolations.append(newLine)
    #print(interpolations)
    return np.array(interpolations)

newBigBenMat = []
bigBenMat = bigBenMat[0]
for i in range(np.shape(bigBenMat)[0]-2):
    newBigBenMat.append(bigBenMat[i])
    if(bigBenMat[i+1][0]-bigBenMat[i][0] > 0.015):
        inter = getInterpolation(bigBenMat[i], bigBenMat[i+1])
        newBigBenMat.append(inter)
    #print(bigBenMat[i][0])

newBigBenMat = np.vstack(newBigBenMat )
#newBigBenMat = np.asarray(newBigBenMat)
#time = newBigBenMat[:,0]
plt.plot(newBigBenMat[:,0], newBigBenMat[:,5], 'k')
plt.plot(uttMat[0][:,0], uttMat[0][:,5], 'r--')
plt.show()

"""
def getCoP_x(uttDataFrame):
    topLeft = uttDataFrame[1]
    topRight = uttDataFrame[2]
    bottomLeft = uttDataFrame[3]
    bottomRight = uttDataFrame[4]
    
    if (topRight+bottomRight+topLeft+bottomLeft) > 5 :
        cop_x = (topRight+bottomRight)-(topLeft+bottomLeft)
        cop_x = cop_x / (topRight+bottomRight+topLeft+bottomLeft)
        cop_x = (270*cop_x)/2        
        return cop_x
    else :
        return 0
    
def getCoP_y(uttDataFrame):
    topLeft = uttDataFrame[1]
    topRight = uttDataFrame[2]
    bottomLeft = uttDataFrame[3]
    bottomRight = uttDataFrame[4]
    
    if (topRight+bottomRight+topLeft+bottomLeft) > 5 :
        cop_x = (topRight+topLeft)-(bottomRight+bottomLeft)
        cop_x = cop_x / (topRight+bottomRight+topLeft+bottomLeft)
        cop_x = (270*cop_x)/2        
        return cop_x
    else :
        return 0
    
cop_x = np.apply_along_axis(getCoP_x, 1, uttMat)
cop_y = np.apply_along_axis(getCoP_y, 1, uttMat)
"""