"""
Model script for R || C || (L + Rl)
"""

import numpy as np

j = 1j

# List of parameter dictionaries with names, initial values,
# and min/max bounds. Set 'vary': False to hold a param constant.
PARAMS = [
    {"name":  "R", "init":   2.00e2, "vary": True, "min": 1e-12, "max":  None},
    {"name":  "C", "init": 20.0e-12, "vary": True, "min": 1e-12, "max":  None},
    {"name":  "L", "init":     1e-6, "vary": True, "min": 1e-12, "max":  None},
    {"name": "Rl", "init":   100e-3, "vary": True, "min": 1e-12, "max":  None},
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
    L = params['L']
    Rl = params['Rl']
    # This is the definition of the model impedance which we want to
    # fit to the data points.  Modify it to represent the circuit you
    # want to fit to the data.
    Zl = j*w*L + Rl
    Zc = 1/(j*w*C)
    Y = 1/R + 1/Zc + 1/Zl
    Z = 1/Y
    return Z
