from msc_gcv import *
import matplotlib.pyplot as plt

traj = trajectory(Param())
for k in 'zs', 'zm', 'zo':
    print('%s = %g m' % (k, traj[k][-1]))
plot(traj)
show()
