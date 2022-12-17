# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 00:47:35 2022

Steven Shockley
homework #1
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.animation as ani
import matplotlib.style as mplstyle

#path = 'C:/Users/steve/Downloads/hw6_temp/'
path = 'C:/Users/steve/OneDrive/Documents/c-projects/steubN-body/data/high_density/'

mplstyle.use(['dark_background', 'fast'])

s = [(10*i)+9 for i in range(80)]
df = pd.read_csv(f'{path}data_s_{s[0]}.csv', names=
                 ['m', 'x', 'y', 'z', 'vx', 'vy', 'vz'], delimiter=' ')
data_frames = [pd.read_csv(f'{path}data_s_{i}.csv', 
                           names=['m', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'ax', 
                                  'ay', 'az'], 
                           delimiter=' ') for i in s]

maxs = np.zeros(len(data_frames))
mins = np.zeros(len(data_frames))

r_inner = 2.1
width = 1.1
ecc = 0.05

min_peri = 2.1*(1-ecc)
max_apo = 3.2*(1+ecc)

j = 0
for df in data_frames:
    r = np.zeros(len(df['x'])-6)
    for i in range(len(df['x'])):
        if i<6: continue
        r[i-6] = np.sqrt(df.at[i, 'x']**2 + df.at[i, 'y']**2 + df.at[i, 'z']**2)
    print(f'maximum: {np.max(r)}\nminimum: {np.min(r)}')
    maxs[j] = np.max(r)
    mins[j] = np.min(r)
    j += 1
    
fig = plt.figure(1, figsize=(12., 8.))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212, sharex=ax1)
ax1.plot(maxs, label='Maximum orbital radius')
ax2.plot(mins, label='Minimum orbital radius')
ax1.hlines(max_apo, 5, 75, colors='lightgreen', label='Expected apihelion')
ax2.hlines(min_peri, 5, 75, colors='peru', label='Expected perihelion')
ax2.set_xlabel('Unit time')
for ax in [ax1, ax2]:
    ax.legend()
    ax.set_ylabel('Radius (AU)')