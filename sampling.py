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

def getCoP_x(uttDataFrame):
    topLeft = uttDataFrame[1]
    topRight = uttDataFrame[2]
    bottomLeft = uttDataFrame[3]
    bottomRight = uttDataFrame[4]
    
    cop_x = (topRight+bottomRight)-(topLeft+bottomLeft)
    cop_x = cop_x / (topRight+bottomRight+topLeft+bottomLeft)
    cop_x = (0.270*cop_x)/2        
    return cop_x
    
def getCoP_y(uttDataFrame):
    topLeft = uttDataFrame[1]
    topRight = uttDataFrame[2]
    bottomLeft = uttDataFrame[3]
    bottomRight = uttDataFrame[4]
    
    cop_x = (topRight+topLeft)-(bottomRight+bottomLeft)
    cop_x = cop_x / (topRight+bottomRight+topLeft+bottomLeft)
    cop_x = (0.270*cop_x)/2        
    return cop_x
    
def getCoP_x2(uttDataFrame):
    topLeft = uttDataFrame[1]
    topRight = uttDataFrame[2]
    bottomLeft = uttDataFrame[3]
    bottomRight = uttDataFrame[4]
    
    cop_x = (topRight+bottomRight)-(topLeft+bottomLeft)
    cop_x = cop_x / (topRight+bottomRight+topLeft+bottomLeft)
    cop_x = (0.440*cop_x)/2        
    return cop_x
    
def getCoP_y2(uttDataFrame):
    topLeft = uttDataFrame[1]
    topRight = uttDataFrame[2]
    bottomLeft = uttDataFrame[3]
    bottomRight = uttDataFrame[4]
    
    cop_x = (topRight+topLeft)-(bottomRight+bottomLeft)
    cop_x = cop_x / (topRight+bottomRight+topLeft+bottomLeft)
    cop_x = (0.260*cop_x)/2        
    return cop_x

utt_cop_x = np.apply_along_axis(getCoP_x, 1, uttMat[0])
utt_cop_y = np.apply_along_axis(getCoP_y, 1, uttMat[0])

bigBen_cop_x = np.apply_along_axis(getCoP_x2, 1, newBigBenMat)
bigBen_cop_y = np.apply_along_axis(getCoP_y2, 1, newBigBenMat)

plt.plot(utt_cop_x, utt_cop_y)
plt.plot(bigBen_cop_x, bigBen_cop_y)
plt.show()

def getAcceleration(pos_vector, i):
    a_x = pos_vector[i-1]+pos_vector[i+1]-(2*pos_vector[i])
    a_x = a_x/(0.01*0.01)
    return a_x

utt_acceleration_x = []
utt_acceleration_y = []

for i in range (1, utt_cop_x.size-1):
    utt_acceleration_x.append(getAcceleration(utt_cop_x, i))
    utt_acceleration_y.append(getAcceleration(utt_cop_y, i))

plt.plot(utt_acceleration_x[100:350])
plt.show()

plt.plot(utt_acceleration_y[100:350])
plt.show()

weight = np.mean(uttMat[0][:,5])
temp = uttMat[0][1:-1,5]

cop_x_utt = utt_cop_x[1:-1]*uttMat[0][1:-1,5]*9.8
cop_y_utt = utt_cop_y[1:-1]*uttMat[0][1:-1,5]*9.8

F_x = uttMat[0][1:-1,5]*np.array(utt_acceleration_x)*0.003
F_y = uttMat[0][1:-1,5]*np.array(utt_acceleration_y)*0.003

F_z = (newBigBenMat[0:761,5] + bigBenObj[0]["weightKg"])*9.8

F_x_z = F_x / F_z

utt_fz = utt_cop_x[1:-1]*F_z
#bigBen_cop_x_th = utt_cop_x[1:-1]*uttMat[0][1:-1,5]*9.8
#bigBen_cop_x_th = bigBen_cop_x_th + (np.array(utt_acceleration_x) * weight *1)
bigBen_cop_x_th = ((utt_cop_x[1:-1]*newBigBenMat[0:761,5]*9.8) - F_x) / F_z
bigBen_cop_y_th = ((utt_cop_y[1:-1]*newBigBenMat[0:761,5]*9.8) - F_y) / F_z

plt.plot(bigBen_cop_x_th, bigBen_cop_y_th)
plt.plot(bigBen_cop_x, bigBen_cop_y)
plt.show()

plt.plot(bigBen_cop_x_th, bigBen_cop_y_th)
plt.plot(utt_cop_x, utt_cop_y)
plt.show()

error_x = bigBen_cop_x_th - bigBen_cop_x[0:761]
error_y = bigBen_cop_y_th - bigBen_cop_y[0:761]

plt.plot(error_x)
plt.show()

plt.plot(error_y)
plt.show()