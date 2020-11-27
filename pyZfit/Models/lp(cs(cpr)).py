"""
Model script for L || (C2 + (C1 || R1))
Good for distributed C in transformer winding.
"""

import numpy as np

j = 1j

# List of parameter dictionaries with names, initial values,
# and min/max bounds. Set 'vary': False to hold a param constant.
PARAMS = [
    {"name": "C1", "init": 20.0e-12, "vary": True, "min":      0.0, "max":  1.00e12},
    {"name": "R1", "init":      100, "vary": True, "min":      0.0, "max":  1.00e12},
    {"name": "C2", "init": 50.0e-12, "vary": True, "min":      0.0, "max":  1.00e12},
    {"name":  "L", "init":   150e-6, "vary": True, "min":      0.0, "max":  1.00e12},
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
    C1 = params['C1']
    R1 = params['R1']
    C2 = params['C2']
    L = params['L']
    # This is the definition of the model impedance which we want to
    # fit to the data points.  Modify it to represent the circuit you
    # want to fit to the data.
    Y1 = j*w*C1 + 1/R1
    Z1 = 1/Y1 + 1/(j*w*C2)
    Y = 1/Z1 + 1/(j*w*L)
    Z = 1 / Y
    return Z
