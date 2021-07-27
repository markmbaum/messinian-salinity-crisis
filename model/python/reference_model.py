from msc_gcv import *
from os.path import join
import matplotlib.pyplot as plt

traj = trajectory(Param())
for k in 'zs', 'zm', 'zo':
    print('%s = %g m' % (k, traj[k][-1]))
plot(traj)
p = join('..', '..', 'plots', 'reference_simulation')
for fmt in ['png', 'ps', 'svg', 'pdf']:
    fn = p + '.' + fmt
    print('saving figure:', fn)
    plt.savefig(fn, format=fmt)
show()
