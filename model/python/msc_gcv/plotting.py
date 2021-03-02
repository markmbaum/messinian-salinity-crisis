import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.pyplot import show

from .physical_system import *

def _handle_X(X, tint, zs0, zm0, method, tol):

    if type(X) is Param:
        traj = trajectory(X, tint, zs0, zm0, method, tol)
    else:
        traj = X

    return(traj)

def _handle_ax(ax):

    if ax is None:
        fig, ax = plt.subplots(1,1)
        return(fig, ax)
    else:
        return(None, ax)

def plot(X, ax=None, legend=True, tint=1e2, zs0=-60.0, zm0=0.0, method='LSODA', tol=1e-7):
    """Plots a trajectory. It's important to remember that the `tint` keyword variable has units of kyr, for convenience.

    :param X: either a trajectory dictionary or a Param object
    :param ax: an existing Axes object to plot in, if desired
    :param legend: whether to show a legend
    :param tint: integration time [kyr]
    :param zs0: initial sill level [m]
    :param zm0: initial Mediterranean level [m]
    :param method: ODE integration method
    :param tol: default tolerance for ODE integrator

    :return: Axes object"""

    traj = _handle_X(X, tint, zs0, zm0, method, tol)

    fig, ax = _handle_ax(ax)

    ax.plot(traj['t']/kyrsec, traj['zs'], color='k', label='$z_s$')
    ax.plot(traj['t']/kyrsec, traj['zm'], color='cornflowerblue', label='$z_m$')
    ax.plot(traj['t']/kyrsec, traj['zo'], color='mediumblue', label='$z_o$')
    if legend:
        ax.legend()
    ax.set_xlabel('Time (kyr)')
    ax.set_ylabel('Level (m)')

    return(ax)

def plot_phase(X, ax=None, legend=True, cmap='nipy_spectral', tint=1e2, zs0=-60.0, zm0=0.0, method='LSODA', tol=1e-7):
    """Plots a trajectory as a phase portrait. It's important to remember that the `tint` keyword variable has units of kyr, for convenience.

    :param X: either a trajectory dictionary or a Param object
    :param ax: an existing Axes object to plot in, if desired
    :param legend: whether to show a legend
    :param cmap: colormap name
    :param tint: integration time [kyr]
    :param zs0: initial sill level [m]
    :param zm0: initial Mediterranean level [m]
    :param method: ODE integration method
    :param tol: default tolerance for ODE integrator

    :return: Axes object"""

    traj = _handle_X(X, tint, zs0, zm0, method, tol)

    fig, ax = _handle_ax(ax)

    x, y = traj['zm'], traj['zs']
    segs = [[(x[i],y[i]),(x[i+1],y[i+1])] for i in range(len(x)-1)]
    lc = LineCollection(segs, cmap=cmap)
    lc.set_array(traj['t'][:-1]/kyrsec)
    ax.add_collection(lc)
    xpad = 0.02*(x.max() - x.min())
    ypad = 0.02*(y.max() - y.min())
    ax.set_xlim(x.min() - xpad, x.max() + xpad)
    ax.set_ylim(y.min() - ypad, y.max() + ypad)
    cb = plt.colorbar(lc, ax=ax)
    cb.set_label('Time (kyr)')
    ax.set_xlabel('$z_m$')
    ax.set_ylabel('$z_s$')

    return(ax)
