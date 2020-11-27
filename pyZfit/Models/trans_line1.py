"""
Model script for a segmented transmission line with skin effect
See LTspice schematic for one segment in ./Data/LineSegment.asc
"""

import numpy as np

j = 1j

# List of parameter dictionaries with names, initial values,
# and min/max bounds. Set 'vary': False to hold a param constant.
PARAMS = [
    {"name":  "Lpm", "init":     1e-6, "vary": True, "min": 1e-12, "max":  None},
    {"name": "Lspm", "init":   200e-9, "vary": True, "min": 1e-12, "max":  None},
    {"name":  "Rpm", "init":   426e-3, "vary": True, "min": 1e-12, "max":  None},
    {"name":  "Cpm", "init": 30.0e-12, "vary": True, "min": 1e-12, "max":  None},
    {"name":  "Gpm", "init":  100e-12, "vary": True, "min": 1e-12, "max":  None},
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
    # Skin effect parallel L/R branches
    BRANCHES = 3

    def segment_z(Zl, w, Lps, Lsps, Rps, Cps, Gps):
        """
        Return impedance array of a cable segment
        See LTspice schem ".\Data\LineSegment.asc"
        :param Zl: load impedance to this segment
        :param w: radian frequency array
        :param Lps: series inductance per segment
        :param Lsps: skin inductance per segment
        :param Rps: series resistance per segment
        :param Cps: shunt capacitance per segment
        :param Gps: shunt conductance per segment
        """
        dL = Lps / 2.0
        dLs = Lsps / 2.0
        dR = Rps / 2.0
        dC = Cps
        dG = Gps
        # Branch factors
        bf = [10**(x/2.0) for x in range(BRANCHES)]
        Ys = np.zeros(len(w), dtype='complex')
        # For each skin effect branch:
        for factor in bf:
            Lb = dLs / factor
            Xb = np.array(1j * w * Lb)
            Rb = dR * factor
            Zb = Rb + Xb
            Ys += 1.0 / Zb
        # Total skin effect impedance
        Zskin = 1.0 / Ys
        # "Termination Z" is load in series with skin effect and incremental L
        Zt = Zl + Zskin + 1j*w*dL
        # "Shunt Z" is parallel combination of Zt and incremental C and G
        Zs = 1.0/(1.0/Zt + 1j*w*dC + dG)
        # Total segment Z is Zs in series with skin effect and incremental L
        return Zs + Zskin + 1j*w*dL

    # Get per-segment values from parameter list
    Lps = params['Lpm']  * LENGTH / SEGMENTS
    Lsps = params['Lspm'] * LENGTH / SEGMENTS
    Rps = params['Rpm'] * LENGTH / SEGMENTS
    Cps = params['Cpm'] * LENGTH / SEGMENTS
    Gps = params['Gpm'] * LENGTH / SEGMENTS
    load = kws['load']
    # Initialize Z array 
    Zmod = load
    # Iterate cable segments, supplying previous Zmod as load to each iteration
    for s in range(SEGMENTS):
        Zmod = segment_z(Zmod, w, Lps, Lsps, Rps, Cps, Gps)
    # Return impedance
    return Zmod
