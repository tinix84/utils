"""
Model script for Ls + Rs + (Cp || Rp)
"""

import numpy as np

j = 1j

# List of parameter dictionaries with names, initial values,
# and min/max bounds. Set 'vary': False to hold a param constant.
PARAMS = [
    {"name": "Ls", "init":   200e-9, "vary": True, "min":      0.0, "max":  1.00e12},
    {"name": "Rs", "init":   500e-3, "vary": True, "min":      0.0, "max":  1.00e12},
    {"name": "Cp", "init":  1.00e-9, "vary": True, "min": 1.00e-12, "max":  1.00e12},
    {"name": "Rp", "init":      200, "vary": True, "min":      0.0, "max":  1.00e12},
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
    Ls = params['Ls']
    Rs = params['Rs']
    Cp = params['Cp']
    Rp = params['Rp']
    # This is the definition of the model impedance which we want to
    # fit to the data points.  Modify it to represent the circuit you
    # want to fit to the data.
    Ycp = j*w*Cp + 1/Rp
    Z = 1/Ycp + j*w*Ls + Rs
    return Z
