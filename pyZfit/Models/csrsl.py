"""
Model script for C + R + L
"""

import numpy as np

# List of parameter dictionaries with names, initial values,
# and min/max bounds. Set 'vary': False to hold a param constant.
PARAMS = [
    {"name": "R", "init":   5.00e1, "vary": True, "min":      0.0, "max":  1.00e12},
    {"name": "C", "init": 30.0e-12, "vary": True, "min": 1.00e-12, "max":  1.00e12},
    {"name": "L", "init":  3.00e-9, "vary": True, "min": 1.00e-12, "max":  1.00e12},
]

j = 1j

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
    L = params['L']
    # This is the definition of the model impedance which we want to
    # fit to the data points.  Modify it to represent the circuit you
    # want to fit to the data.
    return R + 1.0/(j*w*C) + j*w*L
