"""
Model script for output LC filter of buck conveter
"""

import numpy as np


# List of parameter dictionaries with names, initial values,
# and min/max bounds. Set 'vary': False to hold a param constant.
PARAMS = [
    {"name":   "Ls", "init": 1e-6, "vary": True, "min":  100e-9, "max":  10e-3},
    {"name":  "rLs", "init": 1e-3, "vary": True, "min":   100, "max":  1},
    {"name":   "Cp", "init": 100e-9, "vary": True, "min": 10e-9, "max":  100e-6},
    {"name":   "Rp", "init":   1, "vary": True, "min":   1e-3, "max":  10},
    {"name":   "w0", "init":   100e3, "vary": True, "min":  1e3, "max":  10e6},
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
    Ls  = params['Ls']
    rLs   = params['rLs']
    Cp  = params['Cp']
    Rp  = params['Rp']
    load = params['load']
    
    
    ZLs = (1j*w*Ls + rLs)
    ZCp = 1/1j/w/Cp + Rp
                                            
    Zext2int = ZLs * ZCp / (ZLs + ZCp)
    Zint2ext = ZLs + ZCp*load/(ZCp+load)
    attenuation = ZCp / (ZLs + ZCp)
                                            
    return attenuation, Zint2ext, Zext2int


if __name__ == "__main__":
    import numpy as np
    import matplotlib.pyplot as plt

    # Data for plotting
    f = np.linspace(10, 20.0e6, 10000)
    w = 2 * np.pi * f
    params = {'Ls': 1e-6, 'rLs': 1e-3, 'Cp': 100e-6, 'Rp': 10, 'load': 10e6}
    att, zi2e, ze2i = model(w, params)

    # Create figure
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

    # log y axis
    ax1.semilogx(w, 20*np.log10(np.abs(att)))
    ax1.set(title='att')
    ax1.grid()

    # log x axis
    ax2.semilogx(w, 20*np.log10(np.abs(zi2e)))
    ax2.set(title='zi2e')
    ax2.grid()

    # log x axis
    ax3.semilogx(w, 20*np.log10(np.abs(ze2i)))
    ax3.set(title='ze2i')
    ax3.grid()

    fig.tight_layout()
    plt.show()
