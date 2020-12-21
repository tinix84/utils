"""
Model script for Rs + (C || L || Rp)
"""

import numpy as np

j = 1j

# List of parameter dictionaries with names, initial values,
# and min/max bounds. Set 'vary': False to hold a param constant.
PARAMS = [
    {"name": "Rs", "init":     1.00, "vary": True, "min": 1e-12, "max":  1.00e12},
    {"name":  "C", "init": 10.0e-12, "vary": True, "min": 1e-12, "max":  1.00},
    {"name":  "L", "init":   150e-6, "vary": True, "min": 1e-12, "max":  1.00},
    {"name": "Rp", "init":   1.00e3, "vary": True, "min": 1e-12, "max":  1.00e12},
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
    Rs = params['Rs']
    C = params['C']
    L = params['L']
    Rp = params['Rp']
    # This is the definition of the model impedance which we want to
    # fit to the data points.  Modify it to represent the circuit you
    # want to fit to the data.
    Zl = j*w*L
    Zc = 1/(j*w*C)
    Y = 1/Zc + 1/Zl + 1/Rp
    Z = 1/Y + Rs
    return Z
