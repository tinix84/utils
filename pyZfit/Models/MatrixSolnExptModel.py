"""
Model script for matrix solution exp't, with coupled windings.
"""

import numpy as np
import numpy.linalg as la

j = 1j

# List of parameter dictionaries with names, initial values,
# and min/max bounds. Set 'vary': False to hold a param constant.
PARAMS = [
    {"name": "L1", "init": 500e-6 , "vary": True, "min": 1e-12, "max": None},
    {"name": "R1", "init":  20    , "vary": True, "min": 1e-12, "max": None},
    {"name": "L2", "init": 10e-6  , "vary": True, "min": 1e-12, "max": None},
    {"name": "C1", "init": 10e-9  , "vary": True, "min": 1e-12, "max": None},
    {"name": "M" , "init": 4e-6   , "vary": True, "min": 1e-12, "max": None},
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
    L1  = params["L1"]
    R1  = params["R1"]
    L2  = params["L2"]
    C1  = params["C1"]
    M   = params["M"]
    # Calculate the network impedance from a matrix describing node voltages and coupled
    # inductor currents (see "MatrixSolnExpt.ipynb").  The network is driven with input voltage =
    # 1, so network impedance is just 1 / input current.

    # A and B matrices from "MatrixSolnExpt.ipynb":
    def get_A(s):
        # Evaluate A matrix for a single point in s and Zl
        A = np.array([
            [-1/(s*L1) - 1/(s*L2) - 1/R1, 1/(s*L2)        , M/L2 , -M/L1],
            [1/(s*L2)                   , -1/(s*L2) - s*C1, -M/L2, 0    ],
            [1                          , 0               , s*L1 , s*M  ],
            [-1                         , 1               , s*M  , s*L2 ]
        ])
        return A

    def get_B(s):
        # Evaluate B matrix for single s point
        B = np.array([
            [-1/(s*L1)],
            [0],
            [1],
            [0]
        ])
        return B

    # Solution vector = [v2,v3,i1,i2] (see "MatrixSolnExpt.ipynb")
    soln = np.array([la.solve(get_A(j*_w), get_B(j*_w)) for _w in w])
    # List comprehension like this results in transposed matrix, so extract
    # *column* 2 and flatten it from lists of one number to only the number
    i1 = soln[:,2].flatten()
    # Input impedance is source voltage (1.0) over input current
    Z = 1.0/i1

    return Z
