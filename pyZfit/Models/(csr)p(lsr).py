"""
Model script for (C + Rc) || (L + Rl)
"""

import numpy as np

j = 1j

# List of parameter dictionaries with names, initial values,
# and min/max bounds. Set 'vary': False to hold a param constant.
PARAMS = [
    {"name": "Rc", "init":  1.00e-3, "vary": True, "min":      0.0, "max":  1.00e12},
    {"name":  "C", "init": 10.0e-12, "vary": True, "min": 1.00e-12, "max":  1.00e12},
    {"name": "Rl", "init":  1.00e-3, "vary": True, "min": 1.00e-12, "max":  1.00e12},
    {"name":  "L", "init":  1.00e-6, "vary": True, "min":      0.0, "max":  1.00e12},
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
    Rc = params['Rc']
    C  = params['C']
    Rl = params['Rl']
    L  = params['L']
    # This is the definition of the model impedance which we want to
    # fit to the data points.  Modify it to represent the circuit you
    # want to fit to the data.
    Zlr = (j * w * L) + Rl
    Zcr = 1/(j * w * C) + Rc
    Y = 1/Zcr + 1/Zlr
    Z = 1 / Y
    return Z
