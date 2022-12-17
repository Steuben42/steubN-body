import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.animation as ani
import matplotlib.style as mplstyle
from matplotlib.gridspec import GridSpec

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

def find_com(df):
    return np.array([np.sum(df['m']*df['x'])/np.sum(df['m']), 
                     np.sum(df['m']*df['y'])/np.sum(df['m']),
                     np.sum(df['m']*df['z'])/np.sum(df['m'])])

def find_E(df):
    com = find_com(df)
    E = 0
    for i in range(len(df)):
        E += (0.5)/(np.linalg.norm(df.loc[i, 'x':'z'] - com))
    return E
        
E = [find_E(df) for df in data_frames]

fig = plt.figure(1, figsize=(12., 8.))
plt.scatter([i for i in range(len(E))], E)
plt.ylim((0, 2.*np.max(E)))
plt.xlabel('Unit time')
plt.ylabel('Relative energy')