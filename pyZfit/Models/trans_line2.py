"""
Model script for a segmented transmission line with 
arbitrary complex load termination.  Each segment 
is a T with 1/2 LR in series in the top branches and 
CG in shunt between them.
"""

import numpy as np

j = 1j

# List of parameter dictionaries with names, initial values,
# and min/max bounds. Set 'vary': False to hold a param constant.
PARAMS = [
    {"name": "Lpm", "init":   500e-9, "vary": True, "min": 1e-12, "max": None},
    {"name": "Rpm", "init":     1e-0, "vary": True, "min": 1e-12, "max": None},
    {"name": "Cpm", "init": 30.0e-12, "vary": True, "min": 1e-12, "max": None},
    {"name": "Gpm", "init": 100.0e-6, "vary": True, "min": 1e-12, "max": None},
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
    # Length of cable in meters
    LENGTH = 3.251
    # Integer number of segments to split cable into
    SEGMENTS = 3

    def segment_z(Zl, w, Lps, Rps, Cps, Gps):
        """
        Return impedance array of a cable segment
        :param Zl: load impedance to this segment
        :param w: radian frequency array
        :param Lps: series inductance per segment
        :param Rps: series resistance per segment
        :param Cps: shunt capacitance per segment
        """
        # Half of series leg
        Zleg = 0.5*(j*w*Lps + Rps)
        # Termination Z in series with half of series elements
        Zt = Zl + Zleg
        # Now in parallel with segment C and G
        Yt = 1.0/Zt + j*w*Cps + Gps     # Bc is positive
        # In series with remainder of series elements
        Zt = 1.0/Yt + Zleg
        return Zt

    # Get per-segment values from parameter list
    Lps = params['Lpm'] * LENGTH / SEGMENTS
    Rps = params['Rpm'] * LENGTH / SEGMENTS
    Cps = params['Cpm'] * LENGTH / SEGMENTS
    Gps = params['Gpm'] * LENGTH / SEGMENTS
    load = kws['load']
    # Initialize Z array
    Zmod = load
    # Iterate cable segments, supplying previous Zmod as load to each iteration
    for s in range(SEGMENTS):
        Zmod = segment_z(Zmod, w, Lps, Rps, Cps, Gps)
    # Return impedance
    return Zmod
