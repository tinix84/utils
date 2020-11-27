"""
Model script for R || (C + L +R)
"""

import numpy as np

j = 1j

# List of parameter dictionaries with names, initial values,
# and min/max bounds. Set 'vary': False to hold a param constant.
PARAMS = [
    {"name": "Rp", "init":      200, "vary": True, "min": 1e-12, "max": None},
    {"name":  "C", "init": 10.0e-12, "vary": True, "min": 1e-12, "max": None},
    {"name":  "L", "init":     4e-6, "vary": True, "min": 1e-12, "max": None},
    {"name": "Rs", "init":   100e-3, "vary": True, "min": 1e-12, "max": None},
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
    Rp = params['Rp']
    C  = params['C']
    L  = params['L']
    Rs = params['Rs']
    # This is the definition of the model impedance which we want to
    # fit to the data points.  Modify it to represent the circuit you
    # want to fit to the data.
    Zcl = j*w*L + 1/(j*w*C) + Rs
    Y = 1/Rp + 1/Zcl
    Z = 1/Y
    return Z
