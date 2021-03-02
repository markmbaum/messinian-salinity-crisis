from numpy import exp

#parameters of fit
c1 = 2.068e12
α1 = 2754
c2 = 4.035e11
α2 = 127.5

#: world ocean area (not including Mediterranean), https://www.ngdc.noaa.gov/mgg/global/etopo1_ocean_volumes.html
Ao = 360.0e12

fAm = lambda zm: c1*exp(zm/α1) + c2*exp(zm/α2)
fAm.__doc__ = """Mediterranean area as a function of water level. Curve was fit to digitized version of Figure 2 in:

- Meijer, P. & Krijgsman, W. A quantitative analysis of the desiccation and re-filling of the Mediterranean during the Messinian Salinity Crisis. Earth and Planetary Science Letters 240, 510–520 (2005).

:param zm float: level of Mediterranean Sea [m]
:return: surface area of Mediterranean Sea [m :math:`^2`]"""

fzo = lambda zm: (c1*α1/Ao)*(1 - exp(zm/α1)) + (c2*α2/Ao)*(1 - exp(zm/α2))
fzo.__doc__ = """World ocean level, assuming change in Mediterranean level is wholly compensated by change in world ocean level.

:param zm float: level of Mediterranean [m]
:return: level of world ocean [m]"""
