import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np
import itertools


# Importing the full MC data frames for Run 1 and Run 3
df = pd.read_csv("ExperimentSummaryTable_4.csv",delim_whitespace=True)

# Comparison with Bare Board
df_bare = df.query("volt == -500 and board==\"Xenon_Bare\"")
df_16um = df.query("volt == -500 and board==\"Xenon_1.6um\"")
print(df_bare)
print(df_16um)

print(df_bare.columns)


df_bare = df_bare.astype({'temp': int})
df_16um = df_16um.astype({'temp': int})
df_bare = df_bare.sort_values(by='temp',axis=0)
df_16um = df_16um.sort_values(by='temp',axis=0)


fig, axs = plt.subplots(1)
axs.errorbar(df_bare['temp'],df_bare['area'],yerr=df_bare['areaStd'],label=" No Coating",marker="o",markersize=5.,lw=0.5)
axs.errorbar(df_16um['temp'],df_16um['area'],yerr=df_16um['areaStd'],label="ASe Coating",marker="o",markersize=5.,lw=0.5)
axs.set_title('Source: Xenon Flash Lamp -- ASe thickness 1.6 um -- V = -500 V')
axs.set_xlabel(r'temp [K]');
axs.set_ylabel('Pulse Area [V*s]')
plt.grid()
plt.legend()
plt.show()
fig.savefig('CompareBare.png', dpi=fig.dpi)


# Compare different voltages

df500 = df.query("volt == 500 and board==\"Xenon_1.6um\"")
df50m = df.query("volt == -500 and board==\"Xenon_1.6um\"")
df0   = df.query("volt == 0 and board==\"Xenon_1.6um\"")
df250 = df.query("volt == 250 and board==\"Xenon_1.6um\"")
df25m = df.query("volt == -250 and board==\"Xenon_1.6um\"")


df500 = df500.astype({'temp': int})
df500 = df500.sort_values(by='temp',axis=0)
df250 = df250.astype({'temp': int})
df250 = df250.sort_values(by='temp',axis=0)
df0   = df0  .astype({'temp': int})
df0   = df0  .sort_values(by='temp',axis=0)
df50m = df50m.astype({'temp': int})
df50m = df50m.sort_values(by='temp',axis=0)
df25m = df25m.astype({'temp': int})
df25m = df25m.sort_values(by='temp',axis=0)

fig, axs = plt.subplots()
axs.errorbar(df500['temp'], df500['area'],yerr=df500['areaStd'],label="500 V",marker="o",markersize=5.,lw=0.5)
axs.errorbar(df250['temp'], df250['area'],yerr=df250['areaStd'],label="250 V",marker="o",markersize=5.,lw=0.5)
axs.errorbar(df0  ['temp'],   df0['area'],yerr=  df0['areaStd'],label="0 V",color='black',marker="o",markersize=5.,lw=0.5)
axs.errorbar(df25m['temp'], df25m['area'],yerr=df25m['areaStd'],label="-250 V",marker="o",markersize=5.,lw=0.5)
axs.errorbar(df50m['temp'], df50m['area'],yerr=df50m['areaStd'],label="-500 V",marker="o",markersize=5.,lw=0.5)
axs.set_title('Source: Xenon Flash Lamp -- ASe thickness 1.6 um -- T = 293 k')
axs.set_xlabel(r'temp [K]');
axs.set_ylabel('Pulse Area [V*s]')
plt.legend()
plt.grid()
plt.show()
fig.savefig('VoltageCompare.png', dpi=fig.dpi)


'''
df500_3 = df.query("volt == 500 and board==\"Xenon_3.2um\"")
df50m_3 = df.query("volt == -500 and board==\"Xenon_3.2um\"")
df250_3 = df.query("volt == 250 and board==\"Xenon_3.2um\"")
df25m_3 = df.query("volt == -250 and board==\"Xenon_3.2um\"")
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

'''
