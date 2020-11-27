"""
Model script for L + (C || R)
"""

import numpy as np

j = 1j

# List of parameter dictionaries with names, initial values,
# and min/max bounds. Set 'vary': False to hold a param constant.
PARAMS = [
    {"name": "Ls", "init":   100e-9, "vary": True, "min": 1e-12, "max":  None},
    {"name": "Cp", "init":    10e-9, "vary": True, "min": 1e-12, "max":  None},
    {"name": "Rp", "init":      100, "vary": True, "min":     1, "max":  None},
]

def model(w, params, **kwargs):
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
    Cp = params['Cp']
    Rp = params['Rp']
    
    # This is the definition of the model impedance which we want to
    # fit to the data points.  Modify it to represent the circuit you
    # want to fit to the data.
    Ycp = j*w*Cp + 1/Rp
    Z = 1/Ycp + j*w*Ls
    return Z
