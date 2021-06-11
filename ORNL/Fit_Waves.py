import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
import glob
from Read_Scope import *
import numpy as np
from scipy.integrate import simps
from scipy import integrate
from scipy import optimize
from numpy import trapz
import copy
import argparse
from scipy.optimize import curve_fit


# Define the Gaussian function
#def FitFunction(x, A, B, C, D, E,  O):
def FitFunction(x, *p):
    A, B, C, D, E,  O = p
    y = (A*x + O )*np.exp(-1*B*x**2) *np.exp(-1*x/C) *np.exp(x*D)*np.exp(x*E)
    return y

# Define model function to be used to fit to the data above:
def gauss(x, *p):
    #A, mu, sigma, O, C, D, E = p
    A, mu, sigma, O, C, D = p
    return (A*x+O)*np.exp(-(x-mu)**2/(2.*sigma**2))*np.exp(-1*x/C)*np.exp(x*D) #*np.exp(x*E) 

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument("name", nargs='?', default = "Xenon_1.6um", type = str, help="insert first  FileName")
parser.add_argument("temperature", nargs='?', default = "293k", type = str, help="insert first  FileName")
args = parser.parse_args()
SUBFOLDER = args.name
TEMP      = args.temperature # "293k" # <-- temperature subfolder
HEAD      = "/Users/elenag/Desktop/ORNL/ASelenium/ORNL/"+SUBFOLDER+"/%s/" % TEMP 

dictionaryFolders = {"1_500v":"500 V", "2_250v":"250 V", "3_0v":"0 V", "4_-250v":"-250 V", "5_-500v":"-500 V"}

FOLDERS  = glob.glob("/Users/elenag/Desktop/ORNL/ASelenium/ORNL/Xenon_1.6um/293k/5_-500v/")

extraShift =  0.0
print(FOLDERS)


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
        endIntegral[i]   = trapz(filtered_Y2D[i][xLessThanZeroPortion_loc: minimum_loc] )
        cumulative[i] = np.cumsum(filtered_Y2D[i] )

    avg_integrated_signal = np.mean(endIntegral)
    std_integrated_signal = np.std(endIntegral)
    print(BOARD, TEMP, VOLT, avg_integrated_signal, std_integrated_signal)
    
    DAT.Y     = np.mean(filtered_Y2D,axis=0)
    DAT.Y_Std = np.std(filtered_Y2D,axis=0)
    cumY    = np.mean(cumulative,axis=0)
    cumYStd = np.std(cumulative,axis=0)

    stdZero[ff] = DAT.Y
    avg_minimum_loc = int(avg_minimum_loc/250)
    area = trapz( DAT.Y[xLessThanZeroPortion_loc:avg_minimum_loc] )
    
    lw = 1.0
    if "1.6" in BOARD:
        lw = 2.0

    print("Max min", np.argmin(DAT.Y), np.argmax(DAT.Y))
    lowerFitRange = 5145 #5090
    upperFitRange = 9000
    print("lowerFitRange", DAT.X[lowerFitRange] -  DAT.X[lowerFitRange+1],  DAT.X[lowerFitRange])
    #Y1 = filtered_Y2D[:,lowerFitRange:upperFitRange].flatten()
    #X1 = np.repeat(DAT.X[lowerFitRange:upperFitRange],250)
    #print(X1.shape, Y1.shape)
    #exit()
    #p0 = [7886.710672498523, -20327020.898405362, -298563.4144920338, 16843.34936472954, -39709.04267968404, 0.08427022965045368]
    #parameters, covariance = curve_fit(FitFunction, DAT.X[lowerFitRange:upperFitRange], DAT.Y[lowerFitRange:upperFitRange],p0=p0)
    p0 = [1.,DAT.X[np.argmax(DAT.Y)], 3., 4., 5., 6.]#, 7.]
    parameters, covariance = curve_fit(gauss, DAT.X[lowerFitRange:upperFitRange], DAT.Y[lowerFitRange:upperFitRange],p0=p0)
    #parameters, covariance = curve_fit(FitFunction, X1, Y1)
    print("Fit Parameters ", parameters)

    #fit_y = FitFunction(DAT.X[lowerFitRange:upperFitRange], fit_A, fit_B, fit_C, fit_D, fit_E,fit_O)
    fit_y = gauss(DAT.X[lowerFitRange:upperFitRange], *parameters)
    #fit_y = FitFunction(DAT.X, fit_A, fit_B, fit_C, fit_D)
        
    axs.fill_between(DAT.X, DAT.Y-DAT.Y_Std, DAT.Y+DAT.Y_Std,  alpha=0.5)
    axs.plot        (DAT.X, DAT.Y,label=BOARD+" "+TEMP+" "+VOLT+" V", linewidth=1.0)  # plot after 2st offset
    #plt.plot(DAT.X[lowerFitRange:upperFitRange], fit_y, '-', label='fit')
    plt.plot(DAT.X[lowerFitRange:upperFitRange], fit_y, '-', label='fit')
    axs.set_xlabel(r'Time [s]'); axs.set_ylabel('Voltage [V]'); axs.set_title('Xenon Flash Lamp, Bare Boards & 0 V data ')


    


plt.legend()
plt.grid()    
plt.show()
#fig.savefig(TEMP+"WaveForms.png")
exit()

