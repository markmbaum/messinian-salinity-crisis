import sys
from numpy import *
from os.path import join
import matplotlib.pyplot as plt

sys.path.append(join('..', '..', 'python'))
from msc_gcv import *

#-------------------------------------------------------------------------------
#FUNCTIONS

read = lambda p, fn: fromfile(join('..', 'out', 'msc_' + p + '_' + fn))

#-------------------------------------------------------------------------------
#MAIN

#prognostics
traj = {
    'zs': read('gcv', '0'),
    'zm': read('gcv', '1'),
    'zo': fzo(read('gcv', '1')),
    't': read('gcv', 't')
}

#plot
plot(traj)

show()
