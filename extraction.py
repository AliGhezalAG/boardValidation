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

uttStartTime = float(uttObj["general_infos"]["start_time"][17:22])
bigBenStartTime = float(bigBenObj[0]["Horodate"][18:23])

delta_time = uttStartTime-bigBenStartTime

# create data list objects
bigBenDataList = []

for line in bigBenObj:
    data = []
    data.append(float("{0:.2f}".format(line["TIMESTAMP"]+bigBenStartTime-delta_time)))
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
    data.append(line["time_stamp"]+uttStartTime)
    data.append(line["pressure"][0]/g_const)
    data.append(line["pressure"][1]/g_const)
    data.append(line["pressure"][3]/g_const)
    data.append(line["pressure"][2]/g_const)
    uttDataList.append(data)
    
bigBenMat = np.array(bigBenDataList)
uttMat = np.array(uttDataList)

somme = np.sum(uttMat[::, 1:5],axis=1)
uttMat=np.append(uttMat,somme[:,None],axis=1)

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

non_zero_cop_x = np.nonzero(cop_x)[0]
non_zero_cop_y = np.nonzero(cop_y)[0]
non_zero_indexes = np.concatenate((non_zero_cop_x, non_zero_cop_y))
_, i = np.unique(non_zero_indexes, return_index=True)
non_zero_indexes = non_zero_indexes[np.sort(i)]

cop_x = cop_x[non_zero_indexes]
cop_y = cop_y[non_zero_indexes]

plt.plot(cop_x, cop_y)
plt.show()

utt_acceleration_x = []
utt_acceleration_y = []

def getAcceleration(pos_vector):
    a_x = pos_vector[i-1]+pos_vector[i+1]-(2*pos_vector[i])
    a_x = a_x/(0.01*0.01)
    return a_x

for i in range (1, cop_x.size-1):
    utt_acceleration_x.append(getAcceleration(cop_x))
    utt_acceleration_y.append(getAcceleration(cop_y))