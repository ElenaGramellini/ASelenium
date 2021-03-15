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
HEAD      = "/Users/elenag/Desktop/ORNL/ASelenium/ORNL/Photodiode_od4/"

class DATA:
    X = np.array([])
    Y = np.array([])
    Ys = np.array([])
    name = str()
    
Data = dict()
counter = 0
length = 25002 
#print("board temp volt	area lowerTick upperTick lowerTime upperTime")

TMP_FILES = glob.glob(HEAD+"*.trc")
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
area = trapz(smallYPortion) #Data[ct].Ys, [0,0.002], dx=5)
labels = HEAD.split("/")[6:-1]
#lineToFile = labels[0],labels[1],labels[2], area, upperBound, upperBound-lowerBound
#w = (labels[2])[:-1]
#w = (w.split("_"))[1]
try:
    print(labels[0],int((labels[1])[:-1]), area, lowerBound, upperBound, Data[counter].X[lowerBound], Data[counter].X[upperBound])
except ValueError:
    x=1
        
print(area)
plt.plot(Data[counter].X, Data[counter].Ys,label="Full Waveform PhotoDiode")
plt.text(0.0045, 0.4, f'Integral = {area:.2f} [V s]')
#plt.plot(smallXPortion, smallYPortion,'k',label="integral portion PhotoDiode")
plt.xlabel(r'Time [s]'); plt.ylabel('Voltage [V]'); plt.title('Photodiode')
plt.xticks(); plt.yticks(); plt.grid()
#plt.legend(('+500', '+250', '0', '-250', '-500'))
plt.legend()
plt.show()

input()
