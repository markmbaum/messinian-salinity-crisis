"""
This is the brief documentation for the brief python package `msc_gcv`. It's the simplified, consolidated, Python implementation of a prior Messinian Salinity Crisis model:

`Garcia-Castellanos, D. & Villaseñor, A. Messinian salinity crisis regulated by competing tectonics and erosion at the Gibraltar arc. Nature 480, 359–363 (2011) <http://www.nature.com/articles/nature10651>`_

The C++ implementation of the model is in the `model/c++` directory.

To make the model available, the folder containing the package must be in your Python path. The simplest way to achieve this is by appending the proper path to your `sys.path`. For example, if the package is on your computer at the path `/Users/me/model/python`, you simply have to

.. code-block:: python

    import sys
    sys.path.append('/Users/me/model/python')

and you should then be able to run `from msc_gcv import *`

Then, to run the model,

1. construct a :class:`Param` object
2. modify any of the parameters, which are attributes of the :class:`Param` object
3. integrate the system using the :func:`trajectory` function
4. plot the results with the :func:`plot` or :func:`plot_phase` functions

For example, the following six lines of code will use `scipy.solve_ivp <https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html>`_'s `LSODA` wrapper to integrate the model with default parameter values and :math:`a=1.75` for 100 kyr and produce the plot below:

.. code-block:: python

    from msc_gcv import *

    p = Param()
    p.a = 1.75
    traj = trajectory(p)
    plot(p)
    show()

.. image:: trajectory.png

Or, in the simplest demonstration, a simulation with parameters left at the default reference values can be run and plotted with

.. code-block:: python

    from msc_gcv import *

    plot(Param())
    show()

All the classes and functions in the model are documented below!

--------

"""

from .mediterranean_area import Ao, fAm, fzo
from .physical_system import *
from .plotting import plot, plot_phase, show
