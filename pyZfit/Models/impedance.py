"""
Model script for R || C || L
"""
# List of parameter dictionaries with names, initial values,
# and min/max bounds. Set 'vary': False to hold a param constant.
PARAMS = [
    {"name": "R", "init":   1.00e3, "vary": True, "min": 1e-12, "max": None},
    {"name": "C", "init": 40.0e-12, "vary": True, "min": 1e-12, "max": None},
    {"name": "L", "init":   150e-6, "vary": True, "min": 1e-12, "max": None},
]

import numpy as np
import lmfit

# Define one or the other:
LIMITS = "lmfit"
from scipy.stats import gmean   # geometric mean


class ImpedanceModel():
    def __init__(self):
        pass

    @staticmethod
    def model(w, params, **kwargs):
        pass

    def _fit(self, freq, mag, phase, load=None, do_norm_denorm=True, fit_method="leastsq", PARAM_FILE="params.csv"):
        # Radian frequency
        freq = np.array(freq)
        mag = np.array(mag)
        phase = np.array(phase)
        w = 2*np.pi*freq

        # Multiple data segments.
        # Create null weighting array, same size as m but full of 1's
        weight = mag.copy()
        weight.fill(1.0)

        # Complex impedance target
        z = mag*np.exp(1j*phase)

        # Instantiate clean class for lmfit fitter
        params = lmfit.Parameters()
        params.clear()
        # Make working copy of PARAMS list from model
        param_list = list(PARAMS)

        # Init list of name/value tuples
        values = []
        
        # Get selected fitting method
        method = fit_method

        # Adjust min and max if necessary
        for param in param_list:
            param["min"] = self._min_max_set(
                param["min"], method, param["init"] / 1e2)
            param["max"] = self._min_max_set(
                param["max"], method, param["init"] * 1e2)

        if do_norm_denorm:
            # Normalize component, frequency, and impedance values
            # Determine frequency and Z scaling factors from initial values
            fsf, zsf = self._find_sf(param_list)
            # Normalize each component value, min, and max
            for param in param_list:
                type = param["name"][0].upper()
                norm = self._normalize(
                    param["init"], param["min"], param["max"], type, fsf, zsf)
                param["init"] = norm["init"]
                param["min"] = norm["min"]
                param["max"] = norm["max"]
            # Normalize frequency, target Z, and load
            w = w / fsf
            z = z / zsf
            # TODO: check what does gui in case of None
            load = load / zsf
        else:
            fsf, zsf = 1.0, 1.0

        # Add modified params to lmfit Parameter class
        # .add converts min/max of None to -/+inf
        for param in param_list:
            params.add(param["name"], value=param["init"],
                       vary=param["vary"], min=param["min"], max=param["max"])

        # Perform weighted model optimization.
        # Errors will be caught and displayed by zfit_excepthook() in main window.
        kw_args = {"load": load, "fsf": fsf, "zsf": zsf}
        result = lmfit.minimize(self._fcn2min, params, args=(
            w, z, weight), kws=kw_args, method=method)

        # Don't use params class after minimize -- some values are scrambled or changed.

        # Populate values[] with modeling results, denormalized if necessary
        for param in param_list:
            name = param["name"]
            val = result.params[name].value
            if do_norm_denorm:
                comp_type = name[0]
                val = self._denormalize(val, comp_type, fsf, zsf)
            v = (name, val)
            values.append(v)

        if do_norm_denorm:
            # Denormalize frequency, target Z, and load
            w = w * fsf
            z = z * zsf
            # TODO: check what does gui in case of None
            load = load * zsf

        # Write denormalized modeling results to file
        with open(PARAM_FILE, mode='w', encoding='utf-8') as f:
            print('name, value', file=f)
            for param in values:
                print('{}, {}'.format(param[0], param[1]), file=f)
                print('{}, {}'.format(param[0], param[1]))

        # Convert list of tuples to a single dict for the model, to be compatible
        # with the way minimize() uses the model
        values_d = {param[0]: param[1] for param in values}

        # Get complex impedance of model using modeled or locked parameters
        # Use denormalized values
        zfit = Rpcpl.model(w=w, params=values_d)

        # Break into magnitude and degree phase
        Zfit_abs = np.abs(zfit)
        Zdeg_fit = np.angle(zfit, deg=True)

        status = "Number of function calls: " + str(result.nfev)
        if result.aborted:
            status = "Process aborted"

        print(status)
        lmfit.report_fit(result)


    def _fcn2min(self, params, w, Z, weight, log_mag=True, **kwargs):
        """
        This is the function to minimize.  It is the difference between model
        and target (aka residuals) with modeling and penalty weights applied.
        :param params:
        :param w: radian frequency array
        :param Z: complex impedances corresponding to frequencies w
        :param weight: array of weights corresponding to frequencies w
        :param **kwargs: keyword arguments
        :return: must return array for leastsq method, optional for others.
        """

        Zmodel = Rpcpl.model(w=w, params=params, **kwargs)

        if log_mag:
            diff = np.log10(Z) - np.log10(Zmodel)
        else:
            diff = Z - Zmodel
        diff *= weight
        # Flatten complex impedance into re/im adjacent floats
        residuals = diff.view('double')
        return residuals


    # def _min_real_imag(self, params, log_mag=False, **kwargs):
    #     """
    #     This is the function to minimize.  It is the difference between model
    #     and target (aka residuals) with modeling and penalty weights applied.
    #     :param params:
    #     :param w: radian frequency array
    #     :param Z: complex impedances corresponding to frequencies w
    #     :param weight: array of weights corresponding to frequencies w
    #     :param **kwargs: keyword arguments
    #     :return: must return array for leastsq method, optional for others.
    #     """
    #     Zmodel = Rpcpl.model(freq=self.freq, params=params, **kwargs)

    #     if log_mag:
    #         diff = np.log10(self.Zxy) - np.log10(Zmodel)
    #     else:
    #         diff = self.Zxy - Zmodel
    #     diff *= self.weight
    #     # Flatten complex impedance into re/im adjacent floats
    #     residuals = (np.real(diff), np.imag(diff))
    #     return residuals

    def _min_abs_deg(self, params, log_mag=True, **kwargs):
        raise NotImplementedError

    def _min_max_set(self, min_max, method, scaled_val):
        """
        Set min or max for the minimization function being used, and
        whether limits are handled by lmfit or zfit.
        :param min_max:
        :param method:
        :param scaled_val:
        :return:
        """
        if method == "differential_evolution":
            # Diff Evo requires finite min & max values
            return scaled_val if min_max == None else min_max
        # elif LIMITS == "zfit":
        #     # lmfit doesn't do the limiting
        #     return None
        elif LIMITS == "lmfit":
            # lmfit gets spec'd limit
            return min_max

    def _find_sf(self, model_list):
        """
        Find scale factors for frequency (FSF) and impedance (ZSF).
        Ignore component if it doesn't begin with R, G, C, L, or M
        :param model_list: PARAM list from model script
        :return: FSF and ZSF
        """
        R, G, C, L = [], [], [], []
        for comp in model_list:
            type = comp['name'][0].upper()
            if type == 'R':
                R.append(comp['init'])
            elif type == 'G':
                G.append(comp['init'])
            elif type == 'C':
                C.append(comp['init'])
            elif type == 'L' or type == 'M':
                L.append(comp['init'])

        # Get zsf depending on whether Rs and/or Gs are present
        if R == []:
            if G == []:
                # No Rs or Gs
                zsf = 1.0
            else:
                # Gs but no Rs
                zsf = 1.0 / gmean(G)
        else:
            if G == []:
                # Rs but no Gs
                zsf = gmean(R)
            else:
                # Both Rs and Gs.  Use the geometric mean of the
                # preferred zsf for each
                Rzsf = gmean(R)
                Gzsf = 1.0 / gmean(G)
                zsf = gmean([Rzsf, Gzsf])

        # Get fsf depending on whether Ls and/or Cs are present
        if L == []:
            if C == []:
                # No Ls or Cs
                fsf = 1.0
            else:
                # No Ls, but Cs
                fsf = 1.0 / (gmean(C) * zsf)
        else:
            if C == []:
                # No Cs, but Ls
                fsf = zsf / gmean(L)
            else:
                # Both Ls and Cs.  Use the geometric mean of the
                # preferred fsf for each
                Lfsf = zsf / gmean(L)
                Cfsf = 1.0 / (gmean(C) * zsf)
                fsf = gmean([Lfsf, Cfsf])

        return fsf, zsf

    def _normalize(self, val, min, max, comp_type, fsf, zsf):
        # Scale val, min, and max by the appropriate factor depending
        # on its type and return the scaled values.  Default scale is 1.0.
        scale = {
            'R': 1.0/zsf,
            'G': zsf,
            'C': fsf*zsf,
            'L': fsf/zsf,
            'M': fsf/zsf
        }
        s = scale.get(comp_type, 1.0)
        min_ret = None if min is None else min * s
        max_ret = None if max is None else max * s
        return {'init': val*s, 'min': min_ret, 'max': max_ret}

    def _denormalize(self, val, comp_type, fsf, zsf):
        # Scale val by the appropriate factor depending on its type
        # and return the scaled value.  Default scale is 1.0.
        scale = {
            'R': zsf,
            'G': 1.0/zsf,
            'C': 1.0/(fsf*zsf),
            'L': zsf/fsf,
            'M': zsf/fsf
        }
        s = scale.get(comp_type, 1.0)
        return val * s