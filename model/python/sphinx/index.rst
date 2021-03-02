.. msc_gcv documentation master file, created by
   sphinx-quickstart on Fri Jan  8 15:32:14 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

`msc_gcv`
=========

.. automodule:: msc_gcv

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Physical Constants
==================

.. data:: yrsec = 31557600.0

    seconds in a year

.. data:: kyrsec = 31557600000.0

    seconds in 1000 years

.. data:: mmyr = 31557600000.0

    converts m/s to mm/yr

.. data:: ρ = 1e3

    density of water [kg/m^3]

.. data:: g = 9.8

    gravitational acceleration [m/s^2]

Parameter Class
===============

.. autoclass:: Param

Model Functions/Equations
=========================

.. autofunction:: fAm
.. autofunction:: fzo
.. autofunction:: fT
.. autofunction:: fS
.. autofunction:: zsdot
.. autofunction:: zmdot
.. autofunction:: trajectory

Plotting Functions
==================

.. autofunction:: plot
.. autofunction:: plot_phase
