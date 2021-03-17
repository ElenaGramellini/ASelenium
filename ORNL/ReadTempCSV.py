import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np
import itertools


# Importing the full MC data frames for Run 1 and Run 3
df = pd.read_csv("TempCycle.csv",delim_whitespace=True)



'''
df50m_1 = df.query("volt == -500 and board==\"293k\"")
df250_1 = df.query("volt == 250 and board==\"293k\"")
df25m_1 = df.query("volt == -250 and board==\"293k\"")
df50m_2 = df.query("volt == -500 and board==\"293_2\"")
df250_2 = df.query("volt == 250 and board==\"293_2\"")
df25m_2 = df.query("volt == -250 and board==\"293_2\"")
'''

df_1 = df.query("temp==\"293k\"")
df_2 = df.query("temp==\"293_2\"")

fig, axs = plt.subplots()
axs.scatter(df_1['volt'],df_1['area'],label="293k, before cryo cyclying")
axs.scatter(df_2['volt'],df_2['area'],label="293k,  after cryo cyclying")
axs.set_title('Source: Xenon Flash Lamp -- ASe Thickness: 1.6um')
axs.set_xlabel(r'Volt [V]');
axs.set_ylabel('Pulse Area')
plt.grid()
plt.legend()
plt.show()



