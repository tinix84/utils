"""
Model script for Rs + Ls + (Cp || Lp)
"""

import numpy as np

j = 1j

# List of parameter dictionaries with names, initial values,
# and min/max bounds. Set 'vary': False to hold a param constant.
PARAMS = [
    {"name": "Rs", "init":  50.0e-3, "vary": True, "min":  40.0e-3, "max":  1.00e12},
    {"name": "Ls", "init":   100e-9, "vary": True, "min":      0.0, "max":  1.00e12},
    {"name": "Cp", "init": 10.0e-12, "vary": True, "min":      0.0, "max":  1.00e12},
    {"name": "Lp", "init":  10.0e-9, "vary": True, "min":      0.0, "max":  1.00e12},
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
    Ls = params['Ls']
    Cp = params['Cp']
    Lp = params['Lp']
    # This is the definition of the model impedance which we want to
    # fit to the data points.  Modify it to represent the circuit you
    # want to fit to the data.
    Y = j*w*Cp + 1/(j*w*Lp)
    Z = 1/Y + Rs + j*w*Ls
    return Z
