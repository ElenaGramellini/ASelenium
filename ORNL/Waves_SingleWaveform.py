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

FOLDERS = glob.glob("/Users/elenag/Desktop/ORNL/ASelenium/ORNL/Xenon_1.6um/293K/*/")
FOLDERS.sort()
print(FOLDERS)



class DATA:
    X = np.array([])
    Y = np.array([])
    Ys = np.array([])
    Y_Std = np.array([])
    name = str()
    
Data = dict()
counter = 0
length = 25002 
#print("board temp volt	area lowerTick upperTick lowerTime upperTime")
fig, axs = plt.subplots()
for FOLDER in FOLDERS:
    TMP_FILES = glob.glob(FOLDER+"*.trc")
    TMP_FILES.sort()    
    DAT = DATA()
    DAT.X = np.zeros(length); DAT.Y = np.zeros(length); DAT.Y_Std = np.zeros(length)
    yArray2D = np.zeros((len(TMP_FILES),length))

    print(yArray2D.shape)
    print(len(TMP_FILES))
    for i, FILES in enumerate(TMP_FILES):  # Sums all the waveforms in a 
        datX, datY, _ = readTrc(FILES)     # folder and returns the mean... but I also want std...
        DAT.X += datX
        DAT.Y += datY
        yArray2D[i,:] = datY               # make 2D array of voltages, so you can calculate STD for each point after

    DAT.Y_Std = np.std (yArray2D,axis=0)   # calculate std
    DAT.Y /= len(TMP_FILES)                # calculate mean
    DAT.X /= len(TMP_FILES)


    DAT.Ys     = savgol_filter(DAT.Y, 41, 3)    # Applies Savitzky-Golay
    Offset_Loc = np.where(DAT.X<0.0)[0][-10]    # filter to y (voltage)
    Offset     = np.mean(DAT.Y[0:Offset_Loc])   # (data, coefficients, order)
    DAT.Ys -= Offset
    
  
    Data[counter] = DAT

    zero_crossings  = np.where(np.diff(np.sign(Data[counter].Ys)))[0]
    spaces = zero_crossings[1:] - zero_crossings[:-1]
    table = np.where( spaces > 700 )
    
    index = table[0][0]
    ############## this is a terrible hack
    if zero_crossings[index] < 5049:
        index = 12
        #print(zero_crossings[index], index, table, "ooooooooooooooooooooo")


    lowerBound = zero_crossings[index]
    upperBound = zero_crossings[index+1]
    smallYPortion = Data[counter].Ys[lowerBound:upperBound]
    smallXPortion = Data[counter].X[lowerBound:upperBound]
    area = trapz(smallYPortion) #Data[ct].Ys, [0,0.002], dx=5)
    labels = FOLDER.split("/")[6:-1]
    #lineToFile = labels[0],labels[1],labels[2], area, upperBound, upperBound-lowerBound
    w = (labels[3])[:-1]
    w = (w.split("_"))[1]
    #try:
    print(labels[0],(labels[1])[:-1],w, area, lowerBound, upperBound, Data[counter].X[lowerBound], Data[counter].X[upperBound])
    #except ValueError:
    #    x=1
        
    #plt.errorbar(Data[counter].X, Data[counter].Ys, Data[counter].Y_Std)
    #plt.plot(Data[counter].X, Data[counter].Ys,label="full waveform "+labels[2],linewidth=2.0)

    axs.fill_between(Data[counter].X, Data[counter].Ys-Data[counter].Y_Std,Data[counter].Ys+Data[counter].Y_Std, alpha=0.5)
    axs.plot(Data[counter].X, Data[counter].Ys,label="Voltage: "+dictionaryFolders[labels[3]] ,linewidth=1.0)
    #plt.plot(Data[counter].X, Data[counter].Y,label="full waveform "+labels[2],linewidth=1.0)
    #plt.plot(smallXPortion, smallYPortion,'k',label="integral portion "+labels[2],linewidth=2.0)
    axs.set_xlabel(r'Time [s]'); axs.set_ylabel('Voltage [V]'); axs.set_title('Xenon Flash Lamp, T = 293k, ASe Thickness: 1.6 um ')
    plt.xticks(); plt.yticks(); plt.grid()

    counter += 1

plt.legend()
plt.show()
fig.savefig("WaveForms.png")




