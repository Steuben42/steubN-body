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

def find_L(df, j, com):
    r = np.array([df.at[j, 'x'] - com[0], df.at[j, 'y'] - com[1],
                  df.at[j, 'z'] - com[2]]) 
    return np.array([-df.at[j, 'm']*np.dot(r, np.array([df.at[j, 'vx'], 
                                                        df.at[j, 'vy'],
                                                        df.at[j, 'vz']])), 
                     np.dot(r, np.array([df.at[j, 'vx'], df.at[j, 'vy'],
                                         df.at[j, 'vz']]))])
def find_com(df):
    return np.array([np.sum(df['m']*df['x'])/np.sum(df['m']), 
                     np.sum(df['m']*df['y'])/np.sum(df['m']),
                     np.sum(df['m']*df['z'])/np.sum(df['m'])])

sun_L = np.zeros((len(data_frames), 2))
mer_L = np.zeros((len(data_frames), 2))
ven_L = np.zeros((len(data_frames), 2))
ear_L = np.zeros((len(data_frames), 2))
mar_L = np.zeros((len(data_frames), 2))
jup_L = np.zeros((len(data_frames), 2))
ast_L = np.zeros((len(data_frames), 2))
com = np.zeros((len(data_frames), 3))
for i in range(len(data_frames)):
    com[i] = find_com(data_frames[i])
    sun_L[i] = find_L(data_frames[i], 0, com[i])
    mer_L[i] = find_L(data_frames[i], 1, com[i])
    ven_L[i] = find_L(data_frames[i], 2, com[i])
    ear_L[i] = find_L(data_frames[i], 3, com[i])
    mar_L[i] = find_L(data_frames[i], 4, com[i])
    jup_L[i] = find_L(data_frames[i], 5, com[i])
    for j in range(len(data_frames[i]['x'])):
        if j<6: continue
        ast_L[i] += find_L(data_frames[i], j, com[i])
        
offset = -np.min((np.min(sun_L[:, 0]), np.min(mer_L[:, 0]), np.min(ven_L[:, 0]), 
                  np.min(ear_L[:, 0]), np.min(mar_L[:, 0]), np.min(jup_L[:, 0]),
                  np.min(ast_L[:, 0])))*(1.01)      

if offset < 0.0: offset = 0.0  

fig = plt.figure(1, figsize=(12., 12.))
gs = GridSpec(2, 1)
ax1 = fig.add_subplot(gs[0, 0])
ax1.set_yscale('log')
ax2 = fig.add_subplot(gs[1, 0], sharex=ax1)
ax1.plot(sun_L[:, 0]+offset, label='Sun')
ax1.plot(mer_L[:, 0]+offset, label='Mercury')
ax1.plot(ven_L[:, 0]+offset, label='Venus')
ax1.plot(ear_L[:, 0]+offset, label='Earth')
ax1.plot(mar_L[:, 0]+offset, label='Mars')
ax1.plot(jup_L[:, 0]+offset, label='Jupiter')
ax1.plot(ast_L[:, 0]+offset, label='Asteroid Belt')
ax1.hlines(offset, 0, 80, colors='snow', lw=4., label='Zero point')
ax1.legend()
ax2.set_xlabel('Unit time')
ax1.set_ylabel('log(Unit L+min(L)*1.01)')
ax2.set_ylabel('log(Unit L/m+min(L/m)*1.01)')
ax2.plot(sun_L[:, 1], label='Sun')
ax2.plot(mer_L[:, 1], label='Mercury')
ax2.plot(ven_L[:, 1], label='Venus')
ax2.plot(ear_L[:, 1], label='Earth')
ax2.plot(mar_L[:, 1], label='Mars')
ax2.plot(jup_L[:, 1], label='Jupiter')
ax2.plot(ast_L[:, 1], label='Asteroid Belt')
ax2.hlines(offset, 0, 80, colors='snow', lw=4., label='Zero point')

fig2 = plt.figure(2, figsize=(10.,10.))
ax3 = fig2.add_subplot()
ax3.scatter(com[:, 0], com[:, 1], c=s)
ax3.set_xlabel('x (AU)')
ax3.set_ylabel('y (AU)')