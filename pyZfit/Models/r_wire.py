"""
Model script for wire resistance with skin effect.
This model matches data from 
http://chemandy.com/calculators/round-wire-ac-resistance-calculator.htm
closely, for several different wire sizes AWG 40 - 10, frequency 1 kHz - 100 
MHz, and arbitrary lengths.  The website formula is claimed to be very 
accurate, for straight isolated wire. Proximity effect in coils would greatly 
alter the resistance, however. It is not known whether the high-frequency 
resistance would be proportional to the square root of frequency.
 
Length affects Rdc but not knee frequency or shape.  AWG affects knee 
frequency, and shape somewhat.  Temp affects Rdc and weakly affects knee 
shape ('a').  'a' varies with AWG from about 2 for AWG 40 to about 1.2 for 
AWG 10.  Using a fixed 'a' = 2.0 fits most wires very well. 
"""

import numpy as np

# List of parameter dictionaries with names, initial values,
# and min/max bounds. Set 'vary': False to hold a param constant.
PARAMS = [
    {"name": "Rdc", "init":   100e-3, "vary": True, "min": 1e-6, "max":  None},
    {"name":  "w0", "init":    200e3, "vary": True, "min":  1e3, "max":  None},
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
    
    w0 = w0 / fsf
    # This is the definition of the model impedance which we want to
    # fit to the data points.  Modify it to represent the circuit you
    # want to fit to the data.
    #
    # This formula returns in Rdc at low frequencies and resistance 
    # proportional to (w**.485) at high frequencies.  w0 and 'a' adjust the
    # placement and sharpness of the knee.
    # Fixed 'a' of 2.0 fits a wide variety of wires very well.
    a = 2.0
    Z = Rdc*pow(1 + (w/w0)**a, (.485/a))
    return Z
