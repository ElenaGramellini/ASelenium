import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
import glob
from Read_Scope import *
import numpy as np
from scipy.integrate import simps
from scipy import integrate
from numpy import trapz





tau       = 0.00013334972527252184 # RC constant in seconds
pulse     = 0.5                    # light pulse in micro seconds
amplitute = 143.64004507032658     # V * s area under -500V 1.6 um curve

fig, axs = plt.subplots()
for i in range(100):
    mu, sigma = 0.0002, 2e-5 # mean and standard deviation
    s = amplitute*np.random.normal(mu, sigma, 1000)

    x = np.linspace(0,0.002, 1000)
    exponentialTime = s * np.exp(-1*x/tau) 


    axs.plot(x, exponentialTime)

plt.legend()
plt.show()
