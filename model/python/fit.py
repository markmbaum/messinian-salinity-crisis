from numpy import *
from os.path import join
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

#-------------------------------------------------------------------------------
# INPUT

#file with hypsometry points
fnhyp = join('..', '..', 'data', 'hypsometry.csv')

#output file
fnout = join('..', '..', 'plots', 'mediterranean_area_fit')

#-------------------------------------------------------------------------------
# FUNCTIONS

#double exponential to fit
fA = lambda z, c1, α1, c2, α2: c1*exp(z/α1) + c2*exp(z/α2)

#-------------------------------------------------------------------------------
# MAIN

#read the csv
hyp = loadtxt(fnhyp, delimiter=',', skiprows=2)
z, hyp = hyp[:,0], hyp[:,1]

#fit the model
popt, _ = curve_fit(fA, z, hyp, [1e12, 1e3, 1e11, 1e2])
c1, α1, c2, α2 = popt

#print the results
print("Mediterranean area as a function of water level:")
print("  A(z) = %g*exp(z/%g) + %g*exp(z/%g)" % tuple(popt))

#plot the fit and residual
fig, (axa, axb) = plt.subplots(1,2)
zz = linspace(z.min(), z.max(), 250)
axa.plot(hyp/1e12, z, 'k.', label='Mediterranean Hypsometry')
axa.plot(fA(zz, *popt)/1e12, zz, label='Fit')
axa.set_xlabel('Mediterranean Area (10$^{12}$ m$^2$)')
axa.set_ylabel('Mediterranean Level (m)')
axa.legend(fontsize=8)
axb.plot(100*(hyp - fA(z, *popt))/hyp, z, 'r')
axb.set_xlabel('Residual (%)')
axb.set_yticklabels([])
fig.tight_layout()
p = join('..', '..', 'plots', 'mediterranean_area_fit')
for fmt in ['png', 'ps', 'svg', 'pdf']:
    fn = p + '.' + fmt
    print('saving figure:', fn)
    plt.savefig(fn, format=fmt)
plt.show()
