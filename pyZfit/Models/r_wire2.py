"""
Model script for wire resistance with skin effect.
This model matches data from 
http://chemandy.com/calculators/round-wire-ac-resistance-calculator.htm
using a piecewise-linear fit which matches Rdc at low frequencies and
sqrt(w) at high frequencies.  This is slightly faster than r_wire, but
not much.
"""

import numpy as np

# List of parameter dictionaries with names, initial values,
# and min/max bounds. Set 'vary': False to hold a param constant.
PARAMS = [
    {"name": "Rdc", "init":    100e-3, "vary": True, "min": 1e-12, "max":  None},
    {"name":  "w0", "init":    6.0e6, "vary": True, "min": 1e3, "max":  None},
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
    Rdc = params['Rdc']
    w0 = params['w0']
    fsf = kws['fsf']
    zsf = kws['zsf']
    
    w0 = w0 / fsf
    # This is the definition of the model impedance which we want to
    # fit to the data points.  Modify it to represent the circuit you
    # want to fit to the data.
    _Rdc = np.full(w.size, float(Rdc))
    Rhf = np.sqrt(w/w0)/zsf
    Z = np.maximum(_Rdc, Rhf)
    return Z
