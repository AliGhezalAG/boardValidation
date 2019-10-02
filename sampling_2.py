# -*- coding: utf-8 -*-

import json
import numpy as np
import matplotlib.pyplot as plt

class Local_SWARII:
    """
    Implementation of the Sliding Windows Weighted Averaged Interpolation method

    How To use :
        First instantiate the class with the desired parameters
        Then call resample on the desired signal

    """

    def __init__(self, window_size=1, desired_frequency=25):
        """
        Instantiate SWARII

        Parameters :
            desired_frequency : The frequency desired for the output signal,
                                after the resampling.
            window_size : The size of the sliding window, in seconds.
        """
        self.desired_frequency = desired_frequency
        self.window_size = window_size

    def resample(self, time, signal):
        """
        Apply the SWARII to resample a given signal.

        Input :
            time:   The time stamps of the the data point in the signal. A 1-d
                    array of shape n, where n is the number of points in the
                    signal. The unit is seconds.
            signal: The data points representing the signal. A k-d array of
                    shape (n,k), where n is the number of points in the signal,
                    and k is the dimension of the signal (e.g. 2 for a
                    statokinesigram).

        Output:
            resampled_time : The time stamps of the signal after the resampling
            resampled_signal : The resampled signal.
        """

        a_signal = np.array(signal)
        current_time = max(0., time[0])
        # print current_time
        output_time = []
        output_signal = []

        while current_time < time[-1]:

            relevant_times = [t for t in range(len(time)) if abs(
                time[t] - current_time) < self.window_size * 0.5]
            if len(relevant_times) == 0:
                print
                "Trying to interpolate an empty window ! at time ", current_time
                pass

            else:
                if len(relevant_times) == 1:
                    value = a_signal[relevant_times[0]]

                else:
                    value = 0
                    weight = 0

                    for i, t in enumerate(relevant_times):
                        if i == 0 or t == 0:
                            left_border = max(
                                time[0], (current_time - self.window_size * 0.5))

                        else:
                            left_border = 0.5 * (time[t] + time[t - 1])

                        if i == len(relevant_times) - 1:
                            right_border = min(
                                time[-1], current_time + self.window_size * 0.5)
                        else:
                            right_border = 0.5 * (time[t + 1] + time[t])

                        w = right_border - left_border

                        value += a_signal[t] * w
                        weight += w

                    value /= weight
                output_time.append(current_time)
                output_signal.append(value)
            current_time += 1. / self.desired_frequency

        return np.array(output_time), np.array(output_signal)

    @staticmethod
    def purge_artefact(time, signal, threshold_up=2, threshold_down=0.5):
        asignal = np.array(signal)
        nsignal = []
        ntime = []
        n_artefact = 0

        for t in range(1, len(time) - 1):
            if time[t] < 0.1:
                pass
            elif (len(ntime) > 0 and t > 0 and t < len(time) - 1 and np.sum(
                    np.abs(asignal[t + 1] - nsignal[-1])) < threshold_down and np.sum(
                    np.abs(asignal[t] - nsignal[-1])) > threshold_up):
                n_artefact += 1
                pass
            elif (len(ntime) > 0 and t > 1 and t < len(time) - 2 and np.sum(
                    np.abs(asignal[t + 2] - nsignal[-1])) < threshold_down and np.sum(
                    np.abs(asignal[t] - nsignal[-1])) > threshold_up):
                n_artefact += 1
                pass
            else:
                ntime.append(time[t])
                nsignal.append(signal[t])
        print
        "skipped", n_artefact, "artefacts"
        return ntime, nsignal

swarii = Local_SWARII(window_size=3, desired_frequency=100)

# read files
with open('bigBenResults_4.json', 'r') as myBigBenfile:
    bigBenData = myBigBenfile.read()
    
with open('uttResults_4.json', 'r') as myUttfile:
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
    data.append((line["pressure"][0])/g_const)
    data.append((line["pressure"][1])/g_const)
    data.append((line["pressure"][3])/g_const)
    data.append((line["pressure"][2])/g_const)
    uttDataList.append(data)
    
bigBenMat = np.array(bigBenDataList)
uttMat = np.array(uttDataList)

somme = np.sum(uttMat[::, 1:5],axis=1)
uttMat=np.append(uttMat,somme[:,None],axis=1)

bigBenStartTime = np.where(bigBenMat[:,5]>0.25)
uttStartTime = np.where(uttMat[:,5]>9)

uttMat = uttMat[uttStartTime,:]
bigBenMat = bigBenMat[bigBenStartTime,:]

timestamp_array = bigBenMat[0][:, 0]
unique_keys, indices = np.unique(timestamp_array, return_index=True)
bigBenMat = bigBenMat[:, np.sort(indices)]

delta_time = bigBenMat[0][0][0] - uttMat[0][0][0]
bigBenMat[0][:, 0] = bigBenMat[0][:, 0]-delta_time

bigBenEndTime = np.where(bigBenMat[0][:, 0] <= uttMat[0][-1, 0])
bigBenMat = bigBenMat[0][bigBenEndTime,:]

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

"""
plt.plot(nnt, nnsignal, 'k')
plt.plot(utt_nnt, utt_nnsignal, 'r--')
plt.show()
"""

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


nnt, utt_cop_x = swarii.resample(newBigBenMat[:,0], utt_cop_x)
nnt, utt_cop_y = swarii.resample(newBigBenMat[:,0], utt_cop_y)

utt_nnt, bigBen_cop_x = swarii.resample(uttMat[0][:,0], bigBen_cop_x)
utt_nnt, bigBen_cop_y = swarii.resample(uttMat[0][:,0], bigBen_cop_y)

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
    
F_z = (uttMat[0][1:-1,5] - 6.1)*9.8
F_x = -1 * uttMat[0][1:-1,5] * np.array(utt_acceleration_x) * 0.07
F_y = -1 * uttMat[0][1:-1,5] * np.array(utt_acceleration_y) * 0.07

test_x = utt_cop_x[1:-1]*uttMat[0][1:-1,5]*9.8

bigBen_cop_x_th = ((utt_cop_x[1:-1]*uttMat[0][1:-1,5]*9.8) - F_x) / F_z
bigBen_cop_y_th = ((utt_cop_y[1:-1]*uttMat[0][1:-1,5]*9.8) - F_y) / F_z

plt.plot(bigBen_cop_x_th, bigBen_cop_y_th)
plt.plot(bigBen_cop_x, bigBen_cop_y)
plt.show()

error_x = bigBen_cop_x_th - bigBen_cop_x[2:-2]
error_y = bigBen_cop_y_th - bigBen_cop_y[2:-2]

plt.plot(error_x)
plt.show()

plt.plot(error_y)
plt.show()