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
#FOLDERS = glob.glob("/Users/elenag/Desktop/ORNL/ASelenium/ORNL/Xenon_1.6um/*k/*/")
#FOLDERS = glob.glob("/Users/elenag/Desktop/ORNL/ASelenium/ORNL/Xenon_Bare/*/*/")


#FOLDERS.sort()
#FOLDERS = [  "/Users/elenag/Desktop/ORNL/ASelenium/ORNL/Xenon_1.6um/293k/3_0v/",
#             "/Users/elenag/Desktop/ORNL/ASelenium/ORNL/Xenon_1.6um/293k/5_-500v/"]

#FOLDERS = [  "/Users/elenag/Desktop/ORNL/ASelenium/ORNL/Xenon_Bare/297k/*/",
#             "/Users/elenag/Desktop/ORNL/ASelenium/ORNL/Xenon_1.6um/293k/*"]
#160k  230k  293_2 293k  78k   87k
#FOLDERS  = glob.glob("/Users/elenag/Desktop/ORNL/ASelenium/ORNL/Xenon_1.6um/87k/*/")
#FOLDERS2 = glob.glob("/Users/elenag/Desktop/ORNL/ASelenium/ORNL/Xenon_Bare/297k/*/")
#FOLDERS +=  FOLDERS2
extraShift =  0.0
#print(FOLDERS)

#FOLDERS.append(FOLDERS2)
#print(FOLDERS)


'''
FOLDERS = ["/Users/elenag/Desktop/ORNL/ASelenium/ORNL/Xenon_1.6um/78k/5_-500v/",
           "/Users/elenag/Desktop/ORNL/ASelenium/ORNL/Xenon_1.6um/87k/5_-500v/",
           "/Users/elenag/Desktop/ORNL/ASelenium/ORNL/Xenon_1.6um/160k/5_-500v/",
           "/Users/elenag/Desktop/ORNL/ASelenium/ORNL/Xenon_1.6um/230k/5_-500v/",
           "/Users/elenag/Desktop/ORNL/ASelenium/ORNL/Xenon_1.6um/293k/5_-500v/"]
'''

FOLDERS  = glob.glob("/Users/elenag/Desktop/ORNL/ASelenium/ORNL/Xenon_Bare/*/*/")

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
stdZero  = np.zeros((len(FOLDERS),length))
print("board temp volt area areastd")
#print("board temp volt	area lowerTick upperTick lowerTime upperTime")
fig, axs = plt.subplots()
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
        DAT.Y += datY

    # Raw waveforms
    DAT.X /= len(TMP_FILES)
    DAT.Y /= len(TMP_FILES)

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
            filtered_Y2D[i] -= (np.mean(filtered_Y2D[i][:5030]))
            filtered_Y2D[i] += extraShift 
            #endIntegral[i]   = trapz(filtered_Y2D[i][xLessThanZeroPortion_loc: minimum_loc] )
        else:
            minimum          = np.max(filtered_Y2D[i][7049:]) ##############
            minimum_loc      = np.argmax(filtered_Y2D[i][7049:]) ##############
            minimum_loc = -1
            avg_minimum_loc  += minimum_loc
            #filtered_Y2D[i] -= minimum
            filtered_Y2D[i] -= (np.mean(filtered_Y2D[i][:5030])+extraShift) 
            filtered_Y2D[i] *= -1
        zero_crossings = np.where(np.diff(  np.sign(filtered_Y2D[i][6100:]) ) )[0]
        zero_crossingsH = zero_crossings[0]+5100
        zero_crossings = np.where(np.diff(  np.sign(filtered_Y2D[i][:zero_crossingsH]) ) )[0]
        zero_crossingsL = zero_crossings[-1]
        #print(zero_crossingsL,zero_crossingsH)
        endIntegral[i] = trapz(filtered_Y2D[i][zero_crossingsL: zero_crossingsH] )
        cumulative[i]  = np.cumsum(filtered_Y2D[i] )



    DAT.Y     = np.mean(filtered_Y2D,axis=0)
    DAT.Y_Std = np.std(filtered_Y2D,axis=0)
    avg_integrated_signal = np.mean(endIntegral)
    std_integrated_signal = np.std(endIntegral)

    cumY    = np.mean(cumulative,axis=0)
    cumYStd = np.std(cumulative,axis=0)

    stdZero[ff] = DAT.Y
    avg_minimum_loc = int(avg_minimum_loc/250)
    #area = trapz( DAT.Y[5050:zero_crossings] )
    
    lw = 1.0
    if "1.6" in BOARD:
        lw = 1.0
    
    axs.fill_between(DAT.X, DAT.Y-DAT.Y_Std, DAT.Y+DAT.Y_Std,  alpha=0.5)
    axs.plot(DAT.X, DAT.Y,label= "Temp: "+TEMP+", Volt: "+VOLT+" V", linewidth=lw)  # plot after 2st offset
    axs.set_xlabel(r'Time [s]'); axs.set_ylabel('Voltage [V]'); axs.set_title('Xenon Flash Lamp,-500 V data, ASe thickness: 1.6 um  ')
    #axs[1].fill_between(DAT.X, cumY-cumYStd, cumY+cumYStd,  alpha=0.5)
    #axs[1].plot(DAT.X, cumY, label= "Temp: "+TEMP+", Volt: "+VOLT+" V", linewidth=lw)  # plot after 2st offset
    #axs[1].set_xlabel(r'Time [s]'); axs[1].set_ylabel('Integrated Signal [V*s]'); axs[1].set_title('Xenon Flash Lamp,-500 V data, ASe thickness: 1.6 um  ')
    ##axs.plot(DAT.X[zero_crossingsL:zero_crossingsH], DAT.Y[ zero_crossingsL:zero_crossingsH],label=BOARD+" "+TEMP+" "+VOLT+" V", linewidth=lw)  # plot after 2st offset
    #axs[0].grid(b=True, which='major', color='b', linestyle='-')
    #axs[1].grid(b=True, which='major', color='b', linestyle='-')

    print(BOARD, TEMP, VOLT,  np.max(cumY),  cumYStd[np.argmax(cumY)])


    #break

'''

#axs[0].plot(DAT.X[5049:15050] , cumY[5049:15050] , label = "T = "+TEMP)
#axs[0].fill_between(DAT.X[5049:15050], cumY[5049:15050]-cumYStd[5049:15050], cumY[5049:15050]+cumYStd[5049:15050],  alpha=0.50)
#axs[0].plot(DAT.X[5049:15050] , np.cumsum( DAT.Y[5049:15050] ) , label = "T = "+TEMP)  

stdDev = np.std(stdZero,axis=0)
print(stdDev.shape)
axs[1].plot(DAT.X, stdDev,label=VOLT)  # plot after 2st offset

print(zero_crossings)


fig.savefig(TEMP+"WaveForms.png")
'''
plt.legend()
plt.grid()    
plt.show()
exit()
