from numpy import *
import matplotlib.pyplot as plt
from msc_gcv import *

#-------------------------------------------------------------------------------
# INPUT

#parameters to vary and the values to vary them over
pvar = {
    'kb': logspace(-7, -5, 10)/yrsec,
    'tauc': linspace(40, 150, 10),
    'a': linspace(1.1, 2.5, 10)
}

#-------------------------------------------------------------------------------
# MAIN

fig, axs = plt.subplots(1,len(pvar.keys()))

for i,k in enumerate(pvar):
    #fresh parameter set
    p = Param()
    #fresh plot
    for val in pvar[k]:
        #change the parameter value
        p[k] = val
        #solve the system
        traj = trajectory(p)
        #plot the solution
        axs[i].plot(traj['t']/kyrsec, traj['zm'], 'k', linewidth=0.5)
    axs[i].set_title(k)
    yl = axs[i].get_ylim()
    if yl[0] < -1e3:
        axs[i].set_ylim(-1e3, yl[1])

fig.tight_layout()
plt.show()
