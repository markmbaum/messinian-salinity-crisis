from scipy.integrate import solve_ivp

from .mediterranean_area import fAm, fzo

__all__ = [
    'yrsec', 'kyrsec', 'mmyr',
    'ρ', 'g',
    'Param',
    'fT', 'fS',
    'zsdot', 'zmdot',
    'trajectory'
]

#-------------------------------------------------------------------------------
# constants and physical parameters

#seconds in a year
yrsec = 31557600.0
#seconds in 1000 years
kyrsec = 1e3*yrsec
#convert m/s to mm/yr
mmyr = 1e3*yrsec

#density of water [kg/m^3]
ρ = 1e3
#gravitational acceleration [m/s^2]
g = 9.8

#-------------------------------------------------------------------------------
# container for parameters

class Param:
    """Container for parameter values that can be accessed like an object or like a dictionary, just for convenience. Constructor takes no arguments. Parameters have consistent units with meters and seconds. Initialized with reference values (run `print(Param())` to see default values). Parameters/attributes are:

    1. `kb`, erodability constant [m/s/Pa^a]
    2. `tauc`, critical shear stress [Pa]
    3. `Cw`, constant for Turowski channel width formula [?]
    4. `U`, sill uplift rate [m/s]
    5. `a`, erosion exponent [-]
    6. `L`, stream length for slope [m]
    7. `n`, roughness coefficient [s/m^(1/3)]
    8. `P`, precipitation rate [m/s]
    9. `E`, evaporation rate [m/s]
    10. `R`, stream discharge rate (excluding sill/strait) [m^3/s]"""

    def __init__(self):

        #erodability [m/yr/Pa^-a]
        self.kb = 8e-6/yrsec
        #critical shear stress [Pa]
        self.tauc = 50
        #constant for Turowski channel width formula [?]
        self.Cw = 6 #value not clear, but appears to be 6 in DEMO.csh
        #uplift rate [m/s]
        self.U = 4.9/mmyr
        #erosion exponent [-]
        self.a = 1.5
        #stream length [m]
        self.L = 100e3
        #roughness coefficient [s/m^(1/3)]
        self.n = 0.05
        #precipitation rate [m/s]
        self.P = 0.6/yrsec
        #evaporation rate [m/s]
        self.E = 1.2/yrsec
        #input from various rivers [m^3/s]
        self.R = 4500.0 + 12000.0

    def __getitem__(self, key):
        #allow dict-like access
        return(self.__dict__[key])

    def __setitem__(self, key, value):
        #allow dict-like access
        self.__dict__[key] = value

    def __str__(self):
        #print all values in rows
        return('\n'.join(['%s = %g' % (k, self[k]) for k in self.__dict__]))

#-------------------------------------------------------------------------------
# functions evaluating the system of ODEs

def fT(kb, tauc, Cw, U):
    r"""Computes a convenience number from a bundle of parameters in the channel width formula

    :math:`T \equiv C_w (\tau_c + U/k_b)/(\rho g)`

    :param float kb: erodability [m/yr/Pa^-a]
    :param float tauc: critical shear stress [Pa]
    :param float Cw: constant for Turowski channel width formula [?]
    :param float U: uplift rate [m/s]

    :return: :math:`T`"""

    return( Cw*((tauc + U/kb)/(ρ*g))**(-3.0/13.0) )

def fS(zs, zm, p):
    """Computes the sill/channel slope

    :math:`S = (z_o - z_m)/L`

    :param float zs: sill level [m]
    :param float zm: Mediterranean level [m]
    :param float p: Param object

    :return: :math:`S`"""

    #NOTE:
    #The slope is not allowed to be negative here for numerical reasons.
    #If the model behaves physically, the slope should never be negative because
    #the flow velocity over the sill goes to zero at exactly zo = zm. The basin
    #should never fill above the ocean, not even by a tiny amount. In some cases,
    #however, a tiny negative slope appears somewhere in the time stepping
    #process. This is almost certainly just a numerical inaccuracy that occurs
    #when the Mediterranean is filling extremely rapidly. However, if S is ever
    #a negative number, even a vanishingly small one, it blows up the
    #calculation of Q upon the evaluation of S^(13/14). When a NaN appears
    #there, it carries through all subsequent calculations and the model fails.
    #So, even though it should not be necessary on a physical basis, negative
    #slope values are forced to zero for numerical feasibility. This issue
    #cropped up in lots of different ODE solving algorithms, explicity and
    #implicit, so it is not specific to the algorithm used here.
    S = max((fzo(zm) - zm)/p.L, 0)

    return(S)

def zsdot(zs, zm, p):
    r"""Computes the time derivative of the sill level

    :param float zs: sill level [m]
    :param float zm: Mediterranean level [m]
    :param float p: Param object

    :return: :math:`dz_s/dt`"""

    #ocean level
    zo = fzo(zm)
    #channel slope
    S = fS(zs, zm, p)
    #shear stress
    τ = ρ*g*(zo - zs)*S

    return( p.U - p.kb*max(τ - p.tauc, 0)**p.a )

def zmdot(zs, zm, p):
    """Computes the time derivative of the Mediterranean level

    :param float zs: sill level [m]
    :param float zm: Mediterranean level [m]
    :param float p: Param object

    :return: :math:`dz_m/dt`"""

    #ocean level
    zo = fzo(zm)
    #jumble of numbers
    T = fT(p.kb, p.tauc, p.Cw, p.U)
    #channel slope
    S = fS(zs, zm, p)
    #water discharge over the sill
    Q = (T**(13.0/7)/p.n)*(max(zo - zs, 0)**(65.0/21))*(S**(13.0/14))

    return( p.P - p.E + (p.R + Q)/fAm(zm) )

def trajectory(p, tint=1e2, zs0=-60.0, zm0=0.0, method='LSODA', tol=1e-7):
    """Integrates the model, using the set of parameters specified in a Param object. It's important to remember that the `tint` keyword variable has units of kyr, for convenience.

    :param p: Param object
    :param tint: integration time [kyr]
    :param zs0: initial sill level
    :param zm0: initial Mediterranean level
    :param method: integration algorithm for scipy.integrate.solve_ivp
    :param tol: tolerance value for scipy.integrate.solve_ivp

    :return: dictionary with `t`, `zs`, `zm`, and `zo` keys"""

    #ODE system in proper format for scipy
    fun = lambda t, z, p: (zsdot(z[0], z[1], p), zmdot(z[0], z[1], p))
    #solve the system
    sol = solve_ivp(fun, [0.0,tint*kyrsec], [zs0,zm0],
            method=method,
            args=(p,),
            rtol=tol,
            atol=tol)
    #unpack
    traj = dict(t=sol.t, zs=sol.y[0,:], zm=sol.y[1,:], zo=fzo(sol.y[1,:]))

    return(traj)
