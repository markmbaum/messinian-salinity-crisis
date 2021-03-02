from os.path import join
from numpy import *
import matplotlib.pyplot as plt

import sys
sys.path.append(join('..', '..', 'python'))
from msc_gcv import *

#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')

#-------------------------------------------------------------------------------
# INPUT

#output folder
dirout = join('..', 'out')
#plot folder
dirplot = join('..', '..', '..', 'plots')
#number of values for each parameter range
N = 20
#number of parameters
L = 6
#parameter columns
cols = ['a', 'Cw', 'kb', 'L', 'tauc', 'U']
#proper names
lab = dict(zip(cols, ['$a$', '$C_w$', '$\log_{10}(k_b)$', '$L$', r'$\tau_c$', 'U']))
#units
unit = dict(zip(cols, ['-', '-', 'm/yr/Pa$^{-a}$', 'km', 'mm/yr']))

#-------------------------------------------------------------------------------
# FUNCTIONS

def disp(k, x):
    if k == 'kb':
        return( log10(x) )
    elif k == 'L':
        return( x/1e3 )
    else:
        return(x)

def sim(k, x):
    if k == 'kb':
        return( x/yrsec )
    elif k == 'U':
        return( x/mmyr )
    else:
        return(x)

def isim(k, x):
    if k == 'kb':
        return( x*yrsec )
    elif k == 'U':
        return( x*mmyr )
    else:
        return(x)

def savefig(fig, fn):
    fig.savefig(fn)
    print('figure saved: %s' % fn)

#-------------------------------------------------------------------------------
# MAIN

#read the parameter columns
p = {col:fromfile(join(dirout, col), dtype='float32') for col in cols}
#read the classification column
c = fromfile(join(dirout, 'classification'), dtype='int8')

#oscillatory results
m = (c == 0)

#plot the oscillatory solution with the lowest τc and kb values
#appear not to actually oscillate, slipped through checks (prob Newton failure)
for k in ('tauc', 'kb'):
    po = p[k][m].min()
    idx = nonzero(m & (p[k] == po))[0][0]
    param = Param()
    for x in p:
        param[x] = sim(x, p[x][idx])
    plot(param, tint=1e3)
    show()

#histogram of oscillatory results for each param
fig, axs = plt.subplots(2, 3)
axs = list(axs[0,:]) + list(axs[1,:])
for k,ax in zip(p,axs):
    #get bin edges
    x = array(sorted(unique(disp(k, p[k]))))
    h = diff(x)[0]/2
    bins = x - h
    bins = append(bins, [bins[-1] + 2*h])
    #plot histogram
    ax.hist(disp(k, p[k][m]), bins, color='grey')
    ax.set_xlabel(lab[k])
    ax.set_yticks([])
    #plot the reference value
    yl = ax.get_ylim()
    x = Param()[k]
    x = disp(k, isim(k, x))
    ax.plot([x]*2, yl, '--', color='C3')
fig.tight_layout()
savefig(fig, join(dirplot, 'oscillatory_hists'))

fig, axs = plt.subplots(10, 1, figsize=(7.5, 8))
param = Param()
#get a group of random oscillatory solutions
idx = nonzero(m)[0]
idx = idx[random.randint(0, len(idx), len(axs))]
for i in range(len(axs)):
    ax = axs[i]
    #pick a random oscillatory solution and set parameter values
    for k in p:
        param[k] = sim(k, p[k][idx[i]])
    #integrate the system and plot results
    traj = trajectory(param)
    ax.plot(traj['t']/kyrsec, traj['zm'], color='cornflowerblue')
    ax.plot(traj['t']/kyrsec, traj['zs'], color='k')
    if i < len(axs) - 1:
        ax.set_xticklabels([])
    if i == len(axs)//2:
        ax.set_ylabel('Level (m)')
axs[-1].set_xlabel('Time (kyr)')
fig.tight_layout()
savefig(fig, join(dirplot, 'oscillatory_examples'))
