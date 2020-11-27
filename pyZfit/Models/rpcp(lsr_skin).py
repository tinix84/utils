"""
Model script for R || C || (L + Rdc + skin effect)
"""

import numpy as np

j = 1j

# List of parameter dictionaries with names, initial values,
# and min/max bounds. Set 'vary': False to hold a param constant.
PARAMS = [
    {"name":   "C", "init": 10.0e-12, "vary": True, "min": 1.00e-12, "max":  1.00e12},
    {"name":   "L", "init":  50.0e-6, "vary": True, "min":      0.0, "max":  1.00e12},
    {"name": "Rdc", "init":   100e-3, "vary": True, "min":      0.0, "max":  1.00e12},
    {"name":  "sf", "init":  1.00e-3, "vary": True, "min":      0.0, "max":  1.00e12},
    {"name":   "R", "init":    100e3, "vary": True, "min":      0.0, "max":  1.00e12},
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
    # Extract individual component values from params list
    C = params['C']
    L = params['L']
    Rdc = params['Rdc']
    sf = params['sf']
    R = params['R']
    # This is the definition of the model impedance which we want to
    # fit to the data points.  Modify it to represent the circuit you
    # want to fit to the data.
    Zlr = j*w*L + Rdc + sf*np.sqrt(w)
    Y = j*w*C + 1/Zlr + 1/R
    Z = 1/Y
    return Z
