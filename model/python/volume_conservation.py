from numpy import *
from scipy.integrate import quad
import matplotlib.pyplot as plt

from msc_gcv import *

N = 1000
zm = linspace(-1, -2e3, N)
Vm = zeros((N,))
Vo = zeros((N,))

for i in range(N):
    Vm[i] = quad(fAm, zm[i], 0)[0]
    Vo[i] = Ao*fzo(zm[i])

aerr = abs(Vo - Vm)
rerr = aerr/abs(Vo)

fig, (axa, axb) = plt.subplots(1,2)
axa.semilogx(aerr, zm, 'k.')
axa.set_ylabel('$z_m$')
axa.set_xlabel('Absolute Error\n$|V_o - V_m|$')
axb.semilogx(rerr, zm, 'k.')
axb.set_xlabel('Relative Error\n$|(V_o - V_m)/V_o|$')
axb.set_yticklabels([])
fig.tight_layout()
plt.show()
