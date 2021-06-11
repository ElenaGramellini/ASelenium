import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
import glob
from Read_Scope import *
import numpy as np
from scipy.integrate import simps
from scipy import integrate
from numpy import trapz
import copy


import argparse

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument("name", nargs='?', default = "Xenon_1.6um", type = str, help="insert first  FileName")
parser.add_argument("temperature", nargs='?', default = "293k", type = str, help="insert first  FileName")
args = parser.parse_args()
SUBFOLDER = args.name
TEMP      = args.temperature # "293k" # <-- temperature subfolder
HEAD      = "/Users/elenag/Desktop/ORNL/ASelenium/ORNL/"+SUBFOLDER+"/%s/" % TEMP 

dictionaryFolders = {"1_500v":"500 V", "2_250v":"250 V", "3_0v":"0 V", "4_-250v":"-250 V", "5_-500v":"-500 V"}

#FOLDERS = glob.glob("/Users/elenag/Desktop/ORNL/ASelenium/ORNL/Xenon_1.6um/293k/1_500v/")

FOLDERS = glob.glob("/Users/elenag/Desktop/ORNL/ASelenium/ORNL/Xenon_Bare/*/*/")
FOLDERS = glob.glob("/Users/elenag/Desktop/ORNL/ASelenium/ORNL/Xenon_1.6um/293k/*/")
FOLDERS = glob.glob("/Users/elenag/Desktop/ORNL/ASelenium/ORNL/Xenon_3.2um_wow/297k/5_-500v/")

FOLDERS.sort()

'''
print(FOLDERS)
FOLDERS = ["/Users/elenag/Desktop/ORNL/ASelenium/ORNL/Xenon_1.6um/78k/5_-500v/",
           "/Users/elenag/Desktop/ORNL/ASelenium/ORNL/Xenon_1.6um/87k/5_-500v/",
           "/Users/elenag/Desktop/ORNL/ASelenium/ORNL/Xenon_1.6um/160k/5_-500v/",
           "/Users/elenag/Desktop/ORNL/ASelenium/ORNL/Xenon_1.6um/230k/5_-500v/",
           "/Users/elenag/Desktop/ORNL/ASelenium/ORNL/Xenon_1.6um/293k/5_-500v/"]
'''

print("board temp volt area areaStd")

class DATA:
    X = np.array([])
    Y = np.array([])
    Ys = np.array([])
    Ys2 = np.array([])
    Y_Std = np.array([])
    name = str()
    
Data = dict()
counter = 0
length = 25002 # Number of time ticks
#print("board temp volt	area lowerTick upperTick lowerTime upperTime")
fig, axs = plt.subplots(2)
for ff, FOLDER in enumerate(FOLDERS):
    positiveVoltage = True
    BOARD = FOLDER.split('/')[-4]
    TEMP = FOLDER.split('/')[-3]
    VOLT = ((FOLDER.split('/')[-2]).split("_")[1])[:-1]
    #print(VOLT)
    if "-" in  FOLDER.split('/')[-2]:
        positiveVoltage = False
    #print("positiveVoltage",positiveVoltage)
    TMP_FILES = glob.glob(FOLDER+"*.trc")
    TMP_FILES.sort()    
    DAT = DATA()
    DAT.X     = np.zeros(length);
    DAT.Y     = np.zeros(length);
    DAT.Y_Std = np.zeros(length)
    yArray2D  = np.zeros((len(TMP_FILES),length))
    # Read in all waveforms
    for i, FILES in enumerate(TMP_FILES):  # Sums all the waveforms in a 
        datX, datY, _ = readTrc(FILES)     # folder and returns the mean... but I also want std...
        yArray2D[i,:] = datY               # make 2D array of voltages, so you can calculate STD for each point after
        DAT.X += datX

    # Raw waveforms
    DAT.X /= len(TMP_FILES)

    # Let's do some filtering
    filtered_Y2D = copy.deepcopy(yArray2D)
    endIntegral = np.zeros(250)
    cumulative  = np.zeros_like(filtered_Y2D)
    xLessThanZeroPortion_loc    = np.max(np.where( DAT.X  < 0.0))
    avg_minimum_loc = 0    
    for i in range(250):
        filtered_Y2D[i] = savgol_filter(filtered_Y2D[i], 41, 3) #not really sure what this does, but whatevs
        if positiveVoltage:
            minimum          = np.min(filtered_Y2D[i][7049:]) ##############
            minimum_loc      = np.argmin(filtered_Y2D[i][7049:]) ##############
            minimum_loc = -1
            avg_minimum_loc  += minimum_loc
            #filtered_Y2D[i] -= minimum
            filtered_Y2D[i] -= np.mean(filtered_Y2D[i][:5030]) 
            #endIntegral[i]   = trapz(filtered_Y2D[i][xLessThanZeroPortion_loc: minimum_loc] )
        else:
            minimum          = np.max(filtered_Y2D[i][7049:]) ##############
            minimum_loc      = np.argmax(filtered_Y2D[i][7049:]) ##############
            minimum_loc = -1
            avg_minimum_loc  += minimum_loc
            #filtered_Y2D[i] -= minimum
            filtered_Y2D[i] -= np.mean(filtered_Y2D[i][:5030]) 
            filtered_Y2D[i] *= -1
        endIntegral[i]   = trapz(filtered_Y2D[i][xLessThanZeroPortion_loc: minimum_loc] )
        cumulative[i] = np.cumsum(filtered_Y2D[i] )

    avg_integrated_signal = np.mean(endIntegral)
    std_integrated_signal = np.std(endIntegral)
    DAT.Y     = np.mean(filtered_Y2D,axis=0)
    DAT.Y_Std = np.std(filtered_Y2D,axis=0)
    cumY    = np.mean(cumulative,axis=0)
    cumYStd = np.std(cumulative,axis=0)
    
    avg_minimum_loc = int(avg_minimum_loc/250)
    area = trapz( DAT.Y[xLessThanZeroPortion_loc:avg_minimum_loc] )
    #print(avg_integrated_signal, std_integrated_signal, area)
    axs[0].plot(DAT.X, DAT.Y,label=VOLT)  # plot after 2st offset
    axs[1].plot(DAT.X, cumY, label=VOLT)  # plot after 2st offset
    #axs.fill_between(DAT.X, DAT.Y-DAT.Y_Std, DAT.Y+DAT.Y_Std,  alpha=0.5)
    #axs.plot(DAT.X[5049:15050] , cumY[5049:15050] , label = "T = "+TEMP)
    #axs.fill_between(DAT.X[5049:15050], cumY[5049:15050]-cumYStd[5049:15050], cumY[5049:15050]+cumYStd[5049:15050],  alpha=0.50)

    #axs.plot(DAT.X[5049:15050] , np.cumsum( DAT.Y[5049:15050] ) , label = "T = "+TEMP)  
    #axs.set_xlabel(r'Time [s]'); axs.set_ylabel('Integrated Voltage [V*s]'); axs.set_title('Xenon Flash Lamp, ASe Thickness: 1.6 um ')

    plt.legend()

    print(BOARD, TEMP, VOLT, avg_integrated_signal, std_integrated_signal)
plt.grid()    
plt.show()
exit()

    
'''
    print(signal_IntegralY)
    area     = np.mean(signal_IntegralY)
    area_Std = np.std(signal_IntegralY)
    print(area,area_Std)
    print("signal_IntegralY", signal_IntegralY.shape)
    print("yArray2D", yArray2D.shape)
    #maxY = np.max(yArray2D[xLessGreaterZeroPortion], axis = 1)

    
    #Offset     = np.mean(DAT.Y[0:Offset_Loc])   # (data, coefficients, order)




    #plt.show()
    #fig, axs = plt.subplots()
    #axs.plot(xArray2D.T, yArray2D.T)
    
    DAT.Ys2    = savgol_filter(DAT.Y, 41, 3)    # Applies Savitzky-Golay
    DAT.Ys     = savgol_filter(DAT.Y, 41, 3)    # Applies Savitzky-Golay
    Offset_Loc = np.where(DAT.X<0.0)[0][-10]    # filter to y (voltage)
    Offset     = np.mean(DAT.Y[0:Offset_Loc])   # (data, coefficients, order)
    DAT.Ys  -= Offset
    Offset2 = DAT.Y[DAT.X > 0.0001 ].min()
    if sum(DAT.Ys) < 0:
        Offset2 = DAT.Y[DAT.X > 0.0001 ].max()
        DAT.Ys2 -= Offset2
        DAT.Ys2 = -1*DAT.Ys2
    else:
        DAT.Ys2 -= Offset2
    
    Data[counter] = DAT

    ############## this is a terrible hack
    zero_crossings  = np.where(np.diff(np.sign(Data[counter].Ys)))[0]
    spaces = zero_crossings[1:] - zero_crossings[:-1]
    table = np.where( spaces > 700 )
    
    index = table[0][0]

    if zero_crossings[index] < 5049:
        index = 12
        #print(zero_crossings[index], index, table, "ooooooooooooooooooooo")


    lowerBound = zero_crossings[index]
    upperBound = zero_crossings[index+1]
    smallYPortion = Data[counter].Ys[lowerBound:upperBound]
    smallXPortion = Data[counter].X[lowerBound:upperBound]
    area = trapz(smallYPortion) #Data[ct].Ys, [0,0.002], dx=5)
    print("end area",area)
    labels = FOLDER.split("/")[6:-1]
    w = (labels[3])[:-1]
    w = (w.split("_"))[1]
    
    #print(labels[0],(labels[1])[:-1],w, area, lowerBound, upperBound, Data[counter].X[lowerBound], Data[counter].X[upperBound])

    
        
    #plt.errorbar(Data[counter].X, Data[counter].Ys, Data[counter].Y_Std)
    #plt.plot(Data[counter].X, Data[counter].Ys,label="full waveform "+labels[2],linewidth=2.0)


    #axs.plot(Data[counter].X, Data[counter].Y,label="full waveform OG"+labels[2],linewidth=1.0)
    # axs.plot(Data[counter].X, Data[counter].Ys,label="full waveform OG offset" ,linewidth=1.0)
    #axs.plot(Data[counter].X, Data[counter].Ys2,label="full waveform NEW" ,linewidth=1.0)
    #plt.plot(smallXPortion, smallYPortion,'k',label="integral portion "+labels[2],linewidth=2.0)
    axs.set_xlabel(r'Time [s]'); axs.set_ylabel('Voltage [V]'); axs.set_title('Xenon Flash Lamp, T = 293k, ASe Thickness: 1.6 um ')
    plt.xticks(); plt.yticks(); plt.grid()

    counter += 1

plt.legend()
'''
