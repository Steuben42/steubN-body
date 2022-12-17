import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.animation as ani
import matplotlib.style as mplstyle

#path = 'C:/Users/steve/Downloads/hw6_temp/'
# change this path string for a subdirectory containing the csvs
path = 'C:/Users/steve/OneDrive/Documents/c-projects/steubN-body/data/trojan_after_after/'

mplstyle.use(['dark_background', 'fast'])

# change the range value for the number of frames if using bh-mode and 
# not 500 steps
#
# remove +9 if not using bh-mode
#
# change 10 to the output frequency used
s = [(10*i)+9 for i in range(500)]
df = pd.read_csv(f'{path}data_s_{s[0]}.csv', names=
                 ['m', 'x', 'y', 'z', 'vx', 'vy', 'vz'], delimiter=' ')
data_frames = [pd.read_csv(f'{path}data_s_{i}.csv', 
                           names=['m', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'ax', 
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
        else: 
            pts[i].set_color('dimgrey')
            pts[i].set_sizes([6])
            pts[i].set_alpha(0.4)
        

fig = plt.figure(1, figsize=[12., 10.])
ax = fig.add_subplot(projection='3d')
ax.set(xlim3d=(-5.,5.), xlabel='x')
ax.set(ylim3d=(-5.,5.), ylabel='y')
ax.set(zlim3d=(-5.,5.))
plt.axis('off')
ax.view_init(elev=30., azim=45.)

pts = [ax.scatter([],[],[]) for _ in df['x']]
 
plot_pts(data_frames[0], pts)
   
# change interval for a slow/faster animation
fig_ani = ani.FuncAnimation(fig, plot_pts, frames=data_frames, fargs=(pts,), interval=50)

# change this to where you would like to save
fig_ani.save("C:/Users/steve/Downloads/trojan_final.gif")