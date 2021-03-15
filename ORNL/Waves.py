import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
import glob
from Read_Scope import *
import numpy as np
from scipy.integrate import simps
from scipy import integrate
from numpy import trapz

TEMP = '78k' # <-- temperature subfolder
HEAD = "/Users/elenag/Desktop/ORNL/ORNL/Xenon_1.6um/%s/" % TEMP 
FOLDERS = glob.glob(HEAD+"*/")
FOLDERS.sort()
FOLDERS

class DATA:
    X = np.array([])
    Y = np.array([])
    Ys = np.array([])
    name = str()
    
Data = dict()
counter = 0
length = 25002 

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
    counter += 1




# area1 = trapz(Data[0].Ys, [0,0.002], dx=5) 
# area2 = trapz(Data[1].Ys, [0,0.002], dx=5)
# area3 = trapz(Data[2].Ys, [0,0.002], dx=5)
# area4 = trapz(Data[3].Ys, [0,0.002], dx=5)
# area5 = trapz(Data[4].Ys, [0,0.002], dx=5)

# print(TEMP, area1, area2, area3,area4, area5)

#-----------------------------SELECTIVE PLOTS__-------------------------------#
#print(Data[0].Ys)
zero_crossings = np.where(np.diff(np.sign(Data[0].Ys)))[0]
print (zero_crossings)
plt.figure()
ct = 0 # folder number in order 0 = first 
plt.plot(Data[ct].X, Data[ct].Ys,'k')
# #plt.fill_between(Data[ct].X, Data[ct].Y-error, Data[ct].Y+error)
ct = 1
#plt.plot(Data[ct].X, Data[ct].Ys)

ct = 2
#plt.plot(Data[ct].X, Data[ct].Ys)

ct = 3
#plt.plot(Data[ct].X, Data[ct].Ys)

ct = 4
#plt.plot(Data[ct].X, Data[ct].Ys)

# ct = 8
## plt.plot(Data[ct].X, Data[ct].Ys)

# ct = 10
## plt.plot(Data[ct].X, Data[ct].Ys)

plt.xlabel(r'Time [s]'); plt.ylabel('Voltaeg [V]'); plt.title('3.2 um ' + str(TEMP))
plt.xticks(); plt.yticks(); plt.grid()
plt.legend(('+500', '+250', '0', '-250', '-500'))
plt.show()
input()



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
