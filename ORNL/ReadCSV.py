import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np
import itertools


# Importing the full MC data frames for Run 1 and Run 3
df = pd.read_csv("ExperimentSummaryTable.csv",delim_whitespace=True)

df500 = df.query("volt == 500 and board==\"Xenon_1.6um\"")
df50m = df.query("volt == -500 and board==\"Xenon_1.6um\"")
df250 = df.query("volt == 250 and board==\"Xenon_1.6um\"")
df25m = df.query("volt == -250 and board==\"Xenon_1.6um\"")


df500_3 = df.query("volt == 500 and board==\"Xenon_3.2um\"")
df50m_3 = df.query("volt == -500 and board==\"Xenon_3.2um\"")
df250_3 = df.query("volt == 250 and board==\"Xenon_3.2um\"")
df25m_3 = df.query("volt == -250 and board==\"Xenon_3.2um\"")


fig, axs = plt.subplots(2)
axs[0].scatter(df500['temp'],df500['area'],label="500 V, Room Temp")
axs[0].scatter(df50m['temp'],-1*df50m['area'],label="-500 V, Room Temp")
axs[1].scatter(df250['temp'],df250['area'],label="250 V, Room Temp")
axs[1].scatter(df25m['temp'],-1*df25m['area'],label="-250 V, Room Temp")
axs[0].set_title('Xenon_1.6um')
axs[0].set_xlabel(r'temp [K]');
axs[0].set_ylabel('Pulse Area')
plt.legend()
plt.show()

fig, axs = plt.subplots(2)
axs[0].scatter(df500_3['temp'],df500_3['area'],label="500 V, Room Temp")
axs[0].scatter(df50m_3['temp'],-1*df50m_3['area'],label="-500 V, Room Temp")
axs[1].scatter(df250_3['temp'],df250_3['area'],label="250 V, Room Temp")
axs[1].scatter(df25m_3['temp'],-1*df25m_3['area'],label="-250 V, Room Temp")
axs[0].set_title('Xenon_3.2um')
axs[0].set_xlabel(r'temp [K]');
axs[0].set_ylabel('Pulse Area')
plt.legend()
plt.show()



df500 = df.query("temp > 280 and volt == 500")
df50m = df.query("temp > 280 and volt == -500")
df25m = df.query("temp > 280 and volt == -250")
df250 = df.query("temp > 280 and volt == 250")
df0 = df.query("temp > 280 and volt == 0")
print(df500)
plt.scatter(df500['board'],df500['area'],label="500 V, Room Temp")
plt.scatter(df250['board'],df250['area'],label="250 V, Room Temp")
plt.scatter(df0['board'],df0['area'],label="0 V, Room Temp")
plt.scatter(df25m['board'],df25m['area'],label="-250 V, Room Temp")
plt.scatter(df50m['board'],df50m['area'],label="-500 V, Room Temp")
plt.xlabel(r'Board'); plt.ylabel('Pulse Area')
plt.legend()
plt.show()

