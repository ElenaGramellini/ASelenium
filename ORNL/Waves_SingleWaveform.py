import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
import glob
from Read_Scope import *
import numpy as np
from scipy.integrate import simps
from scipy import integrate
from numpy import trapz



import argparse

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument("name", nargs='?', default = "Xenon_1.6um", type = str, help="insert first  FileName")
parser.add_argument("temperature", nargs='?', default = "293k", type = str, help="insert first  FileName")
args = parser.parse_args()
SUBFOLDER = args.name
TEMP      = args.temperature # "293k" # <-- temperature subfolder
HEAD      = "/Users/elenag/Desktop/ORNL/ASelenium/ORNL/"+SUBFOLDER+"/%s/" % TEMP 

dictionaryFolders = {"1_500v":"500 V", "2_250v":"250 V", "3_0v":"0 V", "4_-250v":"-250 V", "5_-500v":"-500 V"}

FOLDERS = glob.glob("/Users/elenag/Desktop/ORNL/ASelenium/ORNL/Xenon_1.6um/293k/1_500v/")
#FOLDERS = glob.glob("/Users/elenag/Desktop/ORNL/ASelenium/ORNL/Xenon_1.6um/293k/5_-500v/")
FOLDERS.sort()
print(FOLDERS)



class DATA:
    X = np.array([])
    Y = np.array([])
    Ys = np.array([])
    Ys2 = np.array([])
    Y_Std = np.array([])
    name = str()
    
Data = dict()
counter = 0
length = 25002 
#print("board temp volt	area lowerTick upperTick lowerTime upperTime")
fig, axs = plt.subplots()
for FOLDER in FOLDERS:
    positiveVoltage = True
    if "-" in  FOLDER.split('/')[-2]:
        positiveVoltage = False
    print("positiveVoltage",positiveVoltage)
    TMP_FILES = glob.glob(FOLDER+"*.trc")
    TMP_FILES.sort()    
    DAT = DATA()
    DAT.X = np.zeros(length); DAT.Y = np.zeros(length); DAT.Y_Std = np.zeros(length)
    yArray2D = np.zeros((len(TMP_FILES),length))
    # Read in all waveforms
    for i, FILES in enumerate(TMP_FILES):  # Sums all the waveforms in a 
        datX, datY, _ = readTrc(FILES)     # folder and returns the mean... but I also want std...
        yArray2D[i,:] = datY               # make 2D array of voltages, so you can calculate STD for each point after
        DAT.X += datX
        DAT.Y += datY        

    DAT.Y /= len(TMP_FILES)                # calculate mean
    DAT.X /= len(TMP_FILES)



    # Single file data Analysis
    xLessThanZeroPortion_loc    = np.max(np.where( DAT.X  < 0.0))
    xGreaterZeroPortion_loc     = np.min(np.where( DAT.X  > 0.0004))
    print(xGreaterZeroPortion_loc)
    print("xLessThanZeroPortion_loc",xLessThanZeroPortion_loc.shape, xGreaterZeroPortion_loc.shape)
    # Normalize to mean first
    minY  = np.mean(yArray2D[ :, : xLessThanZeroPortion_loc-10  ], axis = 1)
    yArray2D -= (minY.reshape(250,1))


    #yArray2D1 = yArray2D - (minY.reshape(250,1))      # plot normalized to baseline
    #minY  = np.min(yArray2D1[ :, xGreaterZeroPortion_loc :  ], axis = 1)
    #yArray2D1 -= (minY.reshape(250,1))
    
    minY  = np.min(yArray2D[ :, xGreaterZeroPortion_loc :  ], axis = 1)
    maxY  = np.max(yArray2D[ :, xGreaterZeroPortion_loc :  ], axis = 1)
    minY_loc  = np.argmin(yArray2D[ :, xGreaterZeroPortion_loc :  ], axis = 1)
    maxY_loc  = np.argmax(yArray2D[ :, xGreaterZeroPortion_loc :  ], axis = 1)
    print("maxY_loc",minY_loc)
    print("xGreaterZeroPortion_loc",xGreaterZeroPortion_loc)
    for i in minY_loc:
        print(DAT.X[i])
    # Calculate Integral after zero
    signal_IntegralY     = np.zeros_like(minY_loc)
    print("signal_IntegralY"    , signal_IntegralY.shape)
    print("Positions, min, max ",np.mean(minY_loc), np.mean(maxY_loc) )
    if positiveVoltage :
        yArray2D -= (minY.reshape(250,1))
        for i, iLoc in enumerate(minY_loc):
            signal_IntegralY[i] = np.sum(yArray2D[i, 5049 : 7650]) #xGreaterZeroPortion_loc
    else:
        yArray2D -= (maxY.reshape(250,1))
        for i, iLoc in enumerate(maxY_loc):
            signal_IntegralY[i] = np.sum(yArray2D[i, xGreaterZeroPortion_loc: iLoc])

    #axs.plot(DAT.X, yArray2D.T)  # plot after 2st offset

    print(signal_IntegralY)
    area     = np.mean(signal_IntegralY)
    area_Std = np.std(signal_IntegralY)
    print(area,area_Std)
    print("signal_IntegralY", signal_IntegralY.shape)
    print("yArray2D", yArray2D.shape)
    #maxY = np.max(yArray2D[xLessGreaterZeroPortion], axis = 1)

    
    #Offset     = np.mean(DAT.Y[0:Offset_Loc])   # (data, coefficients, order)


    DAT.Y_Std = np.std (yArray2D,axis=0)   # calculate std in column

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

    #axs.fill_between(Data[counter].X, Data[counter].Ys2-Data[counter].Y_Std,Data[counter].Ys2+Data[counter].Y_Std, alpha=0.5)
    #axs.plot(Data[counter].X, Data[counter].Y,label="full waveform OG"+labels[2],linewidth=1.0)
    # axs.plot(Data[counter].X, Data[counter].Ys,label="full waveform OG offset" ,linewidth=1.0)
    #axs.plot(Data[counter].X, Data[counter].Ys2,label="full waveform NEW" ,linewidth=1.0)
    #plt.plot(smallXPortion, smallYPortion,'k',label="integral portion "+labels[2],linewidth=2.0)
    axs.set_xlabel(r'Time [s]'); axs.set_ylabel('Voltage [V]'); axs.set_title('Xenon Flash Lamp, T = 293k, ASe Thickness: 1.6 um ')
    plt.xticks(); plt.yticks(); plt.grid()

    counter += 1

plt.legend()
plt.show()





