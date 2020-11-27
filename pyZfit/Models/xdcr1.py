"""
Model script for standard transducer model: (C1 + R1) || (C2 + R2 + L2)
"""

import numpy as np

j = 1j

# List of parameter dictionaries with names, initial values,
# and min/max bounds. Set 'vary': False to hold a param constant.
PARAMS = [
    {"name": "R1", "init":     10.0, "vary": True, "min":  1e-12, "max": None},
    {"name": "C1", "init":  200e-12, "vary": True, "min":  1e-12, "max": None},
    {"name": "R2", "init":      400, "vary": True, "min":  1e-12, "max": None},
    {"name": "C2", "init": 20.0e-12, "vary": True, "min":  1e-12, "max": None},
    {"name": "L2", "init":   100e-6, "vary": True, "min":  1e-12, "max": None},
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
    R1 = params['R1']
    C1 = params['C1']
    R2 = params['R2']
    C2 = params['C2']
    L2 = params['L2']
    # This is the definition of the model impedance which we want to
    # fit to the data points.  Modify it to represent the circuit you
    # want to fit to the data.
    # (series R1, C1) || (series R2, C2, L2) 
    Z1 = R1 + 1/(j*w*C1)
    Z2 = R2 + 1/(j*w*C2) + j*w*L2
    Z = 1/(1/Z1 + 1/Z2)
    return Z
