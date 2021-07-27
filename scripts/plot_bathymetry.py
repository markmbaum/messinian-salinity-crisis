from os.path import join
from copy import deepcopy
from numpy import *
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib_scalebar.scalebar import ScaleBar

#-------------------------------------------------------------------------------
# INPUT

#name of bathymetry file
fn = join('..', 'data', 'EMODnet_bathymetry.asc')

#colorscale choice
cmap_name = 'Blues_r'

#color of land
land_color = 'grey'

#output file path
fnout = join('..', 'plots', 'bathymetry')

#earth's radius [m]
R = 6.371e6

#saved figure formats
fmts = ['png', 'ps', 'svg', 'pdf']

#-------------------------------------------------------------------------------
# FUNCTIONS

def tonum(x):
    try:
        x = float(x)
    except TypeError:
        return(x)
    else:
        if isclose(x, int(x)):
            return(int(x))
        else:
            return(x)

#-------------------------------------------------------------------------------
# MAIN

#read properties from the first 6 rows of the file
with open(fn, 'r') as ifile:
    prop = [ifile.readline().strip().split() for i in range(6)]
for i in range(len(prop)):
    prop[i][1] = tonum(prop[i][1])
prop = dict(prop)

#generate the x and y coordinates
x = linspace(
    prop['XLLCORNER'],
    prop['XLLCORNER'] + prop['CELLSIZE']*prop['NCOLS'],
    prop['NCOLS'] + 1)
y = linspace(
    prop['YLLCORNER'] + prop['CELLSIZE']*prop['NROWS'],
    prop['YLLCORNER'],
    prop['NROWS'] + 1)

#find approximate horizontal distance represented by a pixel
dx = (x[1] - x[0])*(pi/180)*R

#read the bathymetry data
bath = genfromtxt(fn, delimiter=' ', skip_header=6)
#flip the sign
bath *= -1
#mask the positive values
bath[bath > 0] = nan

#generate the colorscale
cmap = get_cmap(cmap_name).copy()
cmap.set_bad(land_color)

#make the plot
fig, ax = plt.subplots(1, 1, figsize=(8,3))
r = ax.pcolormesh(x, y, bath, cmap=cmap)
ax.contour((x[1:] + x[:-1])/2, (y[1:] + y[:-1])/2, bath,
        levels=arange(-1400, 0, 150),
        colors='k',
        linestyles='-',
        linewidths=0.5)
cb = plt.colorbar(r, ax=ax)
cb.set_label('Depth (m)')
cb.set_ticks(cb.get_ticks())
cb.set_ticklabels([int(abs(i)) for i in cb.get_ticks()])
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
#label the sill
ax.annotate('Camarinal Sill', (-5.75,35.92), (-5.6,35.8),
        va='top', ha='center',
        arrowprops=dict(arrowstyle='simple', fc='k'))
#scale bar
sb = ScaleBar(dx, 'km',
    length_fraction=0.5,
    location='upper right',
    border_pad=0.4,
    box_alpha=0.85)
ax.add_artist(sb)

fig.tight_layout()
for fmt in fmts:
    fn = fnout + '.' + fmt
    print('saving figure:', fn)
    fig.savefig(fn, format=fmt)
plt.show()
