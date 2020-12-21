"""
Model script for transformer.
This model includes:
> Primary and secondary finite inductances
> Coupling factor k
> Winding shunt capacitance Cp
> Winding shunt resistance Rp
"""

import numpy as np

j = 1j

# List of parameter dictionaries with names, initial values,
# and min/max bounds. Set 'vary': False to hold a param constant.
PARAMS = [
    {"name": "Lp", "init": 100e-6, "vary": True, "min":  1e-6, "max":  None},
    {"name": "Ls", "init": 100e-6, "vary": True, "min":  1e-6, "max":  None},
    {"name":  "k", "init":     .9, "vary": True, "min":  1e-12, "max":     1},
    {"name": "Cp", "init": 10e-12, "vary": True, "min": 1e-12, "max":  None},
    {"name": "Rp", "init":   10e3, "vary": True, "min":   100, "max":  None},
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
    Lp = params['Lp']
    Ls = params['Ls']
    k  = params['k']
    Cp = params['Cp']
    Rp = params['Rp']
    load = kws['load']
    # This is the definition of the model impedance which we want to
    # fit to the data points.  Modify it to represent the circuit you
    # want to fit to the data.
    #
    # Zp equation below is Eqn (3-21) for the effect of a coupled 
    # secondary on primary impedance from Terman, 1955, pp 57-59.  
    # It is a good transformer model to use when inductances and
    # coupling factor k are sought for use in a SPICE simulation.
    # Zsc is Eqn (3-20) = (w*M)**2/Zs, where Zs is the total secondary
    # impedance (winding plus load) and M = k * sqrt(Lp * Ls).
    Zs = j*w*Ls + load                  # secondary winding plus load          
    Zsc = (w**2 * k**2 * Lp * Ls)/Zs    # Eqn (3-20), coupled secondary
                                        # impedance
    Zp = j*w*Lp + Zsc                   # Eqn (3-21), primary with
                                        # coupled secondary
    Y = 1/Zp + j*w*Cp + 1/Rp
    Zp = 1/Y
    return Zp
