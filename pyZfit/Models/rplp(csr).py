"""
Model script for R || L || (C + Rc)
"""

import numpy as np

j = 1j

# List of parameter dictionaries with names, initial values,
# and min/max bounds. Set 'vary': False to hold a param constant.
PARAMS = [
    {"name":  "R", "init":   10.0e3, "vary": True, "min":      0.0, "max":  1.00e12},
    {"name":  "C", "init": 10.0e-12, "vary": True, "min":      0.0, "max":  1.00e12},
    {"name": "Rc", "init":   1.00e3, "vary": True, "min":      0.0, "max":  1.00e12},
    {"name":  "L", "init":  40.0e-6, "vary": True, "min":      0.0, "max":  1.00e12},
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
    R = params['R']
    C = params['C']
    Rc = params['Rc']
    L = params['L']
    # This is the definition of the model impedance which we want to
    # fit to the data points.  Modify it to represent the circuit you
    # want to fit to the data.
    Zc = 1/(j*w*C) + Rc
    Y = 1/R + 1/Zc + 1/(j*w*L)
    Z = 1/Y
    return Z
