import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd

names = ['m', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'm2', 'x2', 'y2', 'z2', 'vx2', 'vy2', 'vz2']
delim = ' '
files = ['./data/5_ecc_leap_out.csv', './data/5_ecc_rk4_out.csv', './data/9_ecc_leap_out.csv',
         './data/9_ecc_rk4_out.csv']

fig_5 = plt.figure('E=0.5', figsize=[15.,10.])
fig_9 = plt.figure('E=0.9', figsize=[15.,10.])

gs = GridSpec(2, 2)

axes_5 = []
axes_9 = []

for i in range(4):
    axes_5.append(fig_5.add_subplot(gs[i//2,i%2]))
    axes_9.append(fig_9.add_subplot(gs[i//2,i%2]))

fig_5.suptitle('N-Body Simulation with $e=0.5$')
fig_9.suptitle('N-Body Simulation with $e=0.9$')
    
i = 0
for f in files:
    df = pd.read_csv(f, names=names, delimiter=delim)

    x = df['x'].to_numpy()
    y = df['y'].to_numpy()
    z = df['z'].to_numpy()
    x2 = df['x2'].to_numpy()
    y2 = df['y2'].to_numpy()
    z2 = df['z2'].to_numpy()
    m = df['m'].to_numpy()
    m2 = df['m2'].to_numpy()
    vx = df['vx'].to_numpy()
    vy = df['vy'].to_numpy()
    vz = df['vz'].to_numpy()
    vx2 = df['vx2'].to_numpy()
    vy2 = df['vy2'].to_numpy()
    vz2 = df['vz2'].to_numpy()

    r = np.sqrt((x-x2)**2 + (y-y2)**2 + (z-z2)**2)
    v = np.sqrt(vx**2 + vy**2 + vz**2)
    v2 = np.sqrt(vx2**2 + vy2**2 + vz2**2)
    vr = np.dot(2*v, r)/r

    t = np.arange(1., np.size(x)+1., 1.)
    if i<2:
        t *= 0.05*9
    else:
        t *= 0.003*150

    E = 0.5*v**2 + 0.5*v2**2 + 1/r
    if i<2:
        axes_5[(i%2)*2].scatter(vr, r, s=8.)
        axes_5[(i%2)*2 + 1].scatter(t, E, s=8.)
        axes_5[(i%2)*2 + 1].set_yscale('log')
        if i==0:
            axes_5[0].set_title('Phase Diagram of Leapfrog')
            axes_5[0].set_xlabel('$v_r$ (30 km/s)')
            axes_5[0].set_ylabel('$r$ (AU)')
            axes_5[1].set_title('Energy versus Time of Leapfrog')
            axes_5[1].set_xlabel('$t$ (yr)')
            axes_5[1].set_ylabel('$\log E$')
        else:
            axes_5[2].set_title('Phase Diagram of RK4')
            axes_5[2].set_xlabel('$v_r$ (30 km/s)')
            axes_5[2].set_ylabel('$r$ (AU)')
            axes_5[3].set_title('Energy versus Time of RK4')
            axes_5[3].set_xlabel('$t$ (yr)')
            axes_5[3].set_ylabel('$\log E$')
    else:
        axes_9[(i%2)*2].scatter(vr, r, s=8.)
        axes_9[(i%2)*2 + 1].scatter(t, E, s=8.)
        axes_9[(i%2)*2 + 1].set_yscale('log')
        if i==2:
            axes_9[0].set_title('Phase Diagram of Leapfrog')
            axes_9[0].set_xlabel('$v_r$ (30 km/s)')
            axes_9[0].set_ylabel('$r$ (AU)')
            axes_9[1].set_title('Energy versus Time of Leapfrog')
            axes_9[1].set_xlabel('$t$ (yr)')
            axes_9[1].set_ylabel('$\log E$')
        else:
            axes_9[2].set_title('Phase Diagram of RK4')
            axes_9[2].set_xlabel('$v_r$ (30 km/s)')
            axes_9[2].set_ylabel('$r$ (AU)')
            axes_9[3].set_title('Energy versus Time of RK4')
            axes_9[3].set_xlabel('$t$ (yr)')
            axes_9[3].set_ylabel('$\log E$')

    i+=1

plt.tight_layout()
plt.subplots_adjust(hspace=0.25)
plt.show()

fig_5.savefig('./data/5_ecc_plots.png')
fig_9.savefig('./data/9_ecc_plots.png')
