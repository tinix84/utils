"""
Model script for a transmission line with 
arbitrary complex load termination.  The full lossy
equation for input Z is used to solve for per-length
R, L, C, and G (no skin effect).
Calculations from http://en.wikipedia.org/wiki/Transmission_line
"""

import numpy as np

j = 1j

# List of parameter dictionaries with names, initial values,
# and min/max bounds. Set 'vary': False to hold a param constant.
PARAMS = [
    {"name": "Lpu", "init":   700e-9, "vary": True, "min": None, "max": None},
    {"name": "Rpu", "init":     1.00, "vary": True, "min": None, "max": None},
    {"name": "Cpu", "init": 50.0e-12, "vary": True, "min": None, "max": None},
    {"name": "Gpu", "init":  40.0e-6, "vary": True, "min": None, "max": None},
]

def model(w, params, **kws):
    """
    Calculate impedance using equations here for all frequencies w.
    :param w: radian frequency array
    :param params: list of component values to apply to the model equations
    :param kws: dict of optional args (eg load, fsf, zsf)
    :param kws: dict of optional args (eg load, fsf, zsf)
    :param kws: dict of optional args (eg load, fsf, zsf)
    :return: complex impedance array corresponding to freqs w
    """
    # Length of cable in "length units" (meters)
    LENGTH = 3.251
    
    # Get per-unit-length values from parameter list
    Lpu = params['Lpu']
    Rpu = params['Rpu']
    Cpu = params['Cpu']
    Gpu = params['Gpu']
    load = kws['load']
    
    # Pre-calculate terms
    s = Rpu + j*w*Lpu
    p = Gpu + j*w*Cpu
    # Characteristic impedance and propagation constant
    Z0 = np.sqrt(s/p)
    gamma = np.sqrt(s*p)
    # Pre-calculate
    # tanh calc below returns function error for large complex values 
    # (known numpy bug in this version).  Use the expansion below to
    # avoid this.
    # t = np.tanh(gamma*LENGTH)
    x = gamma*LENGTH
    t = (1.0 - np.exp(-2.0*x)) / (1.0 + np.exp(-2.0*x)) 
    # Numerator & denominator of Zin formula
    num = load + Z0*t
    den = Z0 + load*t
    # Input impedance of a lossy line with arbitrary load
    Zin = Z0*num/den
    return Zin
