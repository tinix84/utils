"""
Model script for L + Rdc plus skin effect
"""

import numpy as np

j = 1j

# List of parameter dictionaries with names, initial values,
# and min/max bounds. Set 'vary': False to hold a param constant.
PARAMS = [
    {"name":   "L", "init":  50.0e-6, "vary": True, "min":      0.0, "max":  1.00e12},
    {"name": "Rdc", "init":   100e-3, "vary": True, "min":      0.0, "max":  1.00e12},
    {"name":   "s", "init":  1.00e-6, "vary": True, "min":      0.0, "max":  1.00e12},
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
    # Extract individual component values from params list, in the same
    # order as defined in PARAMS above
    L = params['L']
    Rdc = params['Rdc']
    s = params['s']
    # This is the definition of the model impedance which we want to
    # fit to the data points.  Modify it to represent the circuit you
    # want to fit to the data.
    Z = Rdc + s*np.sqrt(w) + j*w*L
    return Z
