# -*- coding: utf-8 -*-
import json
import pandas as pd
import matplotlib.pyplot as plt

# read files
with open('bigBenResults_2.json', 'r') as myBigBenfile:
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
    
    
# Create the pandas DataFrames
bigBenDataFrame = pd.DataFrame(bigBenDataList, columns = ['timestamp', 
                                                          'topLeft', 
                                                          'topRight', 
                                                          'bottomLeft', 
                                                          'bottomRight', 
                                                          'weight', 
                                                          'X', 
                                                          'Y']) 
    
uttDataFrame = pd.DataFrame(uttDataList, columns = ['timestamp', 
                                                      'topLeft', 
                                                      'topRight', 
                                                      'bottomLeft', 
                                                      'bottomRight']) 
    

#bigBenDataFrame = bigBenDataFrame[bigBenDataFrame["timestamp"].isin(uttDataFrame["timestamp"])]

uttWeightDatalist= list(uttDataFrame)
uttWeightDatalist.remove('timestamp')
uttDataFrame['weight'] = uttDataFrame[uttWeightDatalist].sum(axis=1)
val2substract = uttDataFrame.head(1)['weight'][0]
uttDataFrame['weight'] = uttDataFrame['weight'] - val2substract

uttDataFrame.plot(x ='timestamp', y='weight', kind = 'line')	
bigBenDataFrame.plot(x ='timestamp', y='weight', kind = 'line')

#fig = plt.figure()

#for frame in [uttDataFrame, bigBenDataFrame]:
#    plt.plot(frame['timestamp'], frame['weight'])

#plt.xlim(25,45)
#plt.ylim(0,100)
#plt.show()

def getCoP_x(uttDataFrame):
    topLeft = uttDataFrame['topLeft']
    topRight = uttDataFrame['topRight']
    bottomLeft = uttDataFrame['bottomLeft']
    bottomRight = uttDataFrame['bottomRight']
    
    if (topRight+bottomRight+topLeft+bottomLeft) > 5 :
        cop_x = (topRight+bottomRight)-(topLeft+bottomLeft)
        cop_x = cop_x / (topRight+bottomRight+topLeft+bottomLeft)
        cop_x = (270*cop_x)/2        
        return cop_x
    else :
        return 0
    
def getCoP_y(uttDataFrame):
    topLeft = uttDataFrame['topLeft']
    topRight = uttDataFrame['topRight']
    bottomLeft = uttDataFrame['bottomLeft']
    bottomRight = uttDataFrame['bottomRight']
    
    if (topRight+bottomRight+topLeft+bottomLeft) > 5 :
        cop_x = (topRight+topLeft)-(bottomRight+bottomLeft)
        cop_x = cop_x / (topRight+bottomRight+topLeft+bottomLeft)
        cop_x = (270*cop_x)/2        
        return cop_x
    else :
        return 0
    
def isMeasure(dataFrame, threshold):
    return dataFrame['topLeft'] + dataFrame['topRight'] + dataFrame['bottomLeft'] + dataFrame['bottomRight'] > threshold

uttDataFrame['cop_x'] = uttDataFrame.apply(getCoP_x, axis=1)
uttDataFrame['cop_y'] = uttDataFrame.apply(getCoP_y, axis=1)

bigBenDataFrame['cop_x'] = bigBenDataFrame.apply(getCoP_x, axis=1)
bigBenDataFrame['cop_y'] = bigBenDataFrame.apply(getCoP_y, axis=1)

uttDataFrame = uttDataFrame[isMeasure(uttDataFrame, 5)]
bigBenDataFrame = bigBenDataFrame[isMeasure(bigBenDataFrame, 25)]

uttDataFrame.plot(x ='cop_x', y='cop_y', kind = 'line')
bigBenDataFrame.plot(x ='cop_x', y='cop_y', kind = 'line')
bigBenDataFrame.plot(x ='X', y='Y', kind = 'line')
