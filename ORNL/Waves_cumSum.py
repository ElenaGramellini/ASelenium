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
thick = "3.2"
#HEAD      = "/Users/elenag/Desktop/ORNL/ASelenium/ORNL/"+SUBFOLDER+"/%s/" % TEMP
HEAD    = "/Users/elenag/Desktop/ORNL/ASelenium/ORNL/Xenon_"+thick+"um/*k/*_-500*/" # % TEMP 
FOLDERS = glob.glob(HEAD)
FOLDERS.sort()
folderDic = {}
awful = []
for i in FOLDERS:
    key = int(i.split("/")[8][:-1])
    folderDic[key] = i
    awful.append(key)
awful.sort()
FOLDERS = []
for k in awful:
    FOLDERS.append(folderDic[k])



class DATA:
    X = np.array([])
    Y = np.array([])
    Ys = np.array([])
    name = str()
    
Data = dict()
counter = 0
length = 25002 
#print("board temp volt	area lowerTick upperTick lowerTime upperTime")

fig, axs = plt.subplots(1)  

for FOLDER in FOLDERS:
    TMP_FILES = glob.glob(FOLDER+"*.trc")
    TMP_FILES.sort()    
    DAT = DATA()
    DAT.X = np.zeros(length); DAT.Y = np.zeros(length)    
    for FILES in TMP_FILES:                 # Sums all the waveforms in a 
        datX, datY, _ = readTrc(FILES)      # folder and returns the mean
        DAT.X += datX
        DAT.Y += datY
    DAT.Y /= len(TMP_FILES)
    DAT.X /= len(TMP_FILES)    
    DAT.Ys = savgol_filter(DAT.Y, 41, 3)        # Applies Savitzky-Golay
    Offset_Loc = np.where(DAT.X<0.0)[0][-10]    # filter to y (voltage)
    Offset = np.mean(DAT.Y[0:Offset_Loc])       # (data, coefficients, order)
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
    smallYCum     = np.cumsum(np.abs(smallYPortion))

    area = trapz(smallYPortion) #Data[ct].Ys, [0,0.002], dx=5)
    labels = FOLDER.split("/")[6:-1]
    labels = labels[1][6:]+ " "+ labels[2] + " "+ (labels[3][:-1]).split("_")[1]+ " V"
    print(labels)
    '''
    w = (labels[3])[:-1]
    w = (w.split("_"))[2]
    try:
        print(labels[1],int((labels[2])[:-1]),int(w), area, lowerBound, upperBound, Data[counter].X[lowerBound], Data[counter].X[upperBound])
    except ValueError:
        x=1
    '''    

    axs.plot(smallXPortion, smallYCum,label=labels)
    axs.set_xlabel(r'Time [s]'); axs.set_ylabel('Integrated Voltage [V]');
    axs.set_title('Cumulatives  ')
    plt.xticks(); plt.yticks(); plt.grid()
    counter += 1
    
plt.legend()
plt.show()
fig.savefig("Cumulative"+thick+".png", dpi=fig.dpi)


# area1 = trapz(Data[0].Ys, [0,0.002], dx=5) 
# area2 = trapz(Data[1].Ys, [0,0.002], dx=5)
# area3 = trapz(Data[2].Ys, [0,0.002], dx=5)
# area4 = trapz(Data[3].Ys, [0,0.002], dx=5)
# area5 = trapz(Data[4].Ys, [0,0.002], dx=5)

# print(TEMP, area1, area2, area3,area4, area5)

#-----------------------------SELECTIVE PLOTS__-------------------------------#


'''
#plt.figure()
for ct in range(len(Data) -1 ):
    zero_crossings = np.where(np.diff(np.sign(Data[ct].Ys)))[0]
    index = np.argmax(zero_crossings[1:] - zero_crossings[:-1])
    lowerBound = zero_crossings[index]
    upperBound = zero_crossings[index+1]
    smallYPortion = Data[ct].Ys[lowerBound:upperBound]
    smallXPortion = Data[ct].X[lowerBound:upperBound]
    area1 = trapz(smallYPortion) #Data[ct].Ys, [0,0.002], dx=5)
    print(area1)
    #plt.plot(zero_crossings)
    #plt.show()
    #input()
    #ct = 0 # folder number in order 0 = first

    plt.plot(Data[ct].X, Data[ct].Ys,label="full waveform")
    plt.plot(smallXPortion, smallYPortion,'k',label="integral portion")
    plt.xlabel(r'Time [s]'); plt.ylabel('Voltaeg [V]'); plt.title('1.6 um ' + str(TEMP))
    plt.xticks(); plt.yticks(); plt.grid()
    #plt.legend(('+500', '+250', '0', '-250', '-500'))
    plt.legend()
    plt.show()
    #input()
'''
# #plt.fill_between(Data[ct].X, Data[ct].Y-error, Data[ct].Y+error)
#ct = 1
#plt.plot(Data[ct].X, Data[ct].Ys)

#ct = 2
#plt.plot(Data[ct].X, Data[ct].Ys)

#ct = 3
#plt.plot(Data[ct].X, Data[ct].Ys)

#ct = 4
#plt.plot(Data[ct].X, Data[ct].Ys)

# ct = 8
## plt.plot(Data[ct].X, Data[ct].Ys)

# ct = 10
## plt.plot(Data[ct].X, Data[ct].Ys)





#------------------------------MULTI-PLOT-------------------------------------#
# plt.figure()
# for ct in range(1,11):
#     plt.plot(Data[ct].X*1e6, -1*Data[ct].Ys*1e3, label = "-"+str(ct*100)+"V")

# plt.legend(loc='upper right',ncol=3)
# plt.xlabel(r'Time [$\mu$s]'); plt.ylabel('Voltaeg [mV]')
# plt.xticks(); plt.yticks(); plt.grid()
# plt.xlim(-100, 900)


#-------------------------------AREA PLOT------------------------------------#
# area = []
# bias = []
# WindowL = 0
# WindowR = 0.0001
# for ct in range(0,5):
#     LP = np.where(Data[ct].X<WindowL)[0][-1]
#     RP = np.where(Data[ct].X<WindowR)[0][-1]
#     bias.append(ct*100)
#     area.append(np.sum(-1*Data[ct].Ys[LP:RP]))
    
# plt.figure()
# plt.scatter(bias, area, s=50, color='k')
# plt.xlabel(r'Bias [-V]'); plt.ylabel('Area [Vs]')
# plt.xticks(); plt.yticks(); plt.grid()
# plt.xlim(-10, 410)
