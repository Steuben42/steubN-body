import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.animation as ani
import matplotlib.style as mplstyle
from matplotlib.gridspec import GridSpec

#path = 'C:/Users/steve/Downloads/hw6_temp/'
path = 'C:/Users/steve/OneDrive/Documents/c-projects/steubN-body/data/long_solar/'

mplstyle.use(['dark_background', 'fast'])

s = [(10*i)+9 for i in range(1072)]
df = pd.read_csv(f'{path}data_s_{s[0]}.csv', names=
                 ['m', 'x', 'y', 'z', 'vx', 'vy', 'vz'], delimiter=' ')
data_frames = [pd.read_csv(f'{path}data_s_{i}.csv', 
                           names=['m', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'ax1', 
                                  'ay', 'az'], 
                           delimiter=' ') for i in s]



def test(frame, pts):
    i = 0
    for point in pts:
        print((frame.at[i, 'x'],frame.at[i, 'y']))
        print(frame.at[i, 'z'])
        ++i


def plot_pts(frame, pts):
    for i in range(len(pts)):
        pts[i].set_offsets((frame.at[i, 'x'],frame.at[i, 'y']))
        pts[i].set_3d_properties(frame.at[i, 'z'], zdir='z')
        if i==0: 
            pts[i].set_color('gold')
            pts[i].set_sizes([42])
        elif i==1: pts[i].set_color('lightgrey')
        elif i==2: pts[i].set_color('navajowhite')
        elif i==3: pts[i].set_color('skyblue')
        elif i==4: pts[i].set_color('coral')
        elif i==5: pts[i].set_color('burlywood')
        elif i==6: pts[i].set_color('bisque')
        elif i==7: pts[i].set_color('cornflowerblue')
        elif i==8: pts[i].set_color('mediumturquoise')
        

fig = plt.figure(1, figsize=[8., 6.])
ax1 = fig.add_subplot(projection='3d')
ax1.set(xlim3d=(-5.,5.), xlabel='x')
ax1.set(ylim3d=(-5.,5.), ylabel='y')
ax1.set(zlim3d=(-5.,5.))
plt.axis('off')
ax1.view_init(elev=30., azim=45.)

pts = [ax1.scatter([],[],[]) for _ in range(7)]
   
fig_ani = ani.FuncAnimation(fig, plot_pts, frames=data_frames, fargs=(pts,), interval=200)
fig_ani.save("C:/Users/steve/Downloads/long_solar.gif")