"""
Model script for (Lp1 || Rp) + (Lp2 || Rp2)
"""

import numpy as np

j = 1j

# List of parameter dictionaries with names, initial values,
# and min/max bounds. Set 'vary': False to hold a param constant.
PARAMS = [
    {"name": "Rp1", "init":      100, "vary": True, "min":      0.0, "max":  1.00e12},
    {"name": "Lp1", "init":  50.0e-9, "vary": True, "min":      0.0, "max":  1.00e12},
    {"name": "Rp2", "init":      100, "vary": True, "min":      0.0, "max":  1.00e12},
    {"name": "Lp2", "init":  50.0e-9, "vary": True, "min":      0.0, "max":  1.00e12},
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
    Rp1 = params['Rp1']
    Lp1 = params['Lp1']
    Rp2 = params['Rp2']
    Lp2 = params['Lp2']
    # This is the definition of the model impedance which we want to
    # fit to the data points.  Modify it to represent the circuit you
    # want to fit to the data.
    Y1 = 1/Rp1 + 1/(j*w*Lp1)
    Y2 = 1/Rp2 + 1/(j*w*Lp2)
    Z = 1/Y1 + 1/Y2
    return Z
