"""
Model script for a 3-stage LC filter with stray coupling between windings.
See LTspice schematic in ./Data/CoupledFilter.asc
"""

import numpy as np
import numpy.linalg as la

j = 1j

# List of parameter dictionaries with names, initial values,
# and min/max bounds. Set 'vary': False to hold a param constant.
PARAMS = [
    {"name": "Rs1", "init": 20e-3  , "vary": True, "min": 1e-12 , "max": None},
    {"name": "Rp1", "init": 200    , "vary": True,  "min": 1e-12, "max": None},
    {"name": "Rp2", "init": 200    , "vary": True,  "min": 1e-12, "max": None},
    {"name": "Rs2", "init": 60e-3  , "vary": True, "min": 1e-12 , "max": None},
    {"name": "L1" , "init": 820e-9 , "vary": True,  "min": 1e-12, "max": None},
    {"name": "L2" , "init": 820e-9 , "vary": True,  "min": 1e-12, "max": None},
    {"name": "L3" , "init": 50e-9  , "vary": True,  "min": 1e-12, "max": None},
    {"name": "M1" , "init": 1e-12  , "vary": False, "min": 1e-12 , "max": None},
    {"name": "M2" , "init": 1e-12  , "vary": False, "min": 1e-12 , "max": None},
    {"name": "C1" , "init": 1e-9   , "vary": True, "min": 1e-12 , "max": None},
    {"name": "C2" , "init": 470e-12, "vary": True,  "min": 1e-12, "max": None},
    {"name": "C3" , "init": 240e-12, "vary": True,  "min": 1e-12, "max": None},
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
    Rs1 = params["Rs1"]
    Rp1 = params["Rp1"]
    Rp2 = params["Rp2"]
    Rs2 = params["Rs2"]
    L1  = params["L1" ]
    L2  = params["L2" ]
    L3  = params["L3" ]
    M1  = params["M1" ]
    M2  = params["M2" ]
    C1  = params["C1" ]
    C2  = params["C2" ]
    C3  = params["C3" ]
    load = kws['load']

    # Calculate the network impedance from a matrix describing node voltages and coupled
    # inductor currents (see "MatrixGen.py").  The network is driven with input voltage =
    # 1, so network impedance is just 1 / input current.

    # A and B matrices from "MatrixGen.py":
    def get_A(s, Zl):
        # Evaluate A matrix for a single point in s and Zl
        A = np.array([
            [1/Rs1+1/Rp1+1/(L1*s), -1/Rp1-1/(L1*s)           , 0                   , 0                        , 0              , 0                 , 0     , -M1/L1, 0     ],
            [1/Rp1+1/(L1*s)      , -C1*s-1/Rs1-1/Rp1-1/(L1*s), 1/Rs1               , 0                        , 0              , 0                 , 0     , -M1/L1, 0     ],
            [0                   , -1/Rs1                    , 1/Rs1+1/Rp2+1/(L2*s), -1/Rp2-1/(L2*s)          , 0              , 0                 , -M1/L2, 0     , -M2/L2],
            [0                   , 0                         , -1/Rp2-1/(L2*s)     , C2*s+1/Rs2+1/Rp2+1/(L2*s), -1/Rs2         , 0                 , M1/L2 , 0     , M2/L2 ],
            [0                   , 0                         , 0                   , 1/Rs2                    , -1/Rs2-1/(L3*s), 1/(L3*s)          , 0     , M2/L3 , 0     ],
            [0                   , 0                         , 0                   , 0                        , -1/(L3*s)      , C3*s+1/Zl+1/(L3*s), 0     , M2/L3 , 0     ],
            [-1                  , 1                         , 0                   , 0                        , 0              , 0                 , L1*s  , M1*s  , 0     ],
            [0                   , 0                         , -1                  , 1                        , 0              , 0                 , M1*s  , L2*s  , M2*s  ],
            [0                   , 0                         , 0                   , 0                        , -1             , 1                 , 0     , M2*s  , L3*s  ]
        ])
        return A

    # B matrix is constant over s and Zl
    B_row = np.array([1/Rs1, 0, 0, 0, 0, 0, 0, 0, 0])
    # Convert single row to column
    B = B_row[:,None]

    # Solution vector = [v1,v2,v3,v4,v5,v6,i3,i7,i10] (see "MatrixGen.py")
    soln = np.array([la.solve(get_A(j*_w, _Zl), B) for _w, _Zl in zip(w, load)])
    # List comprehension like this results in transposed matrix, so extract
    # *column* 0 and flatten it from lists of one number to only the number
    v1 = soln[:,0].flatten()
    # Input current is difference between source (1.0) and v1, through Rs1
    i_in = (1.0-v1)/Rs1
    # Input impedance is source voltage (1.0) over input current
    Z = 1.0/i_in

    return Z
