import numpy as np
from importlib import import_module, reload

from lmfit import minimize, Parameters, Minimizer, report_fit
from scipy.stats import gmean   # geometric mean
# import csv
# import pandas as pd

# indices for weight, magnitude, and phase lists
M, P, W = 0, 1, 2
# Set max for model weighting.  Minimum is 1.0
MAX_MODEL_WEIGHT = 100.0
# Boolean for whether to allow negative param results
ALLOW_NEG = True
# File name for parameter save/restore
PARAM_FILE = "params.csv"
# List of fitting methods available to lmfit
# (Dogleg and Newton CG are not included since they require a Jacobian
# to be supplied)
METHODS = [
    ("Levenberg-Marquardt", "leastsq"),
    ("Least-Squares minimization, using Trust Region Reflective method", "least_squares"),
    ("differential evolution", "differential_evolution"),
    ("brute force method", "brute"),
    ("basinhopping", "basinhopping"),
    ("Adaptive Memory Programming for Global Optimization", "ampgo"),
    ("Nelder-Mead", "nelder"),
    ("L-BFGS-B", "lbfgsb"),
    ("Powell", "powell"),
    ("Conjugate-Gradient", "cg"),
    ("Newton-CG", "newton"),
    ("Cobyla", "cobyla"),
    ("BFGS", "bfgs"),
    ("Truncated Newton", "tnc"),
    ("Newton-CG trust-region", "trust-ncg"),
    ("nearly exact trust-region", "trust-exact"),
    ("Newton GLTR trust-region", "trust-krylov"),
    ("trust-region for constrained optimization", "trust-constr"),
    ("Dog-leg trust-region", "dogleg"),
    ("Sequential Linear Squares Programming", "slsqp"),
    ("Maximum likelihood via Monte-Carlo Markov Chain", "emcee"),
    ("Simplicial Homology Global Optimization", "shgo"),
    ("Dual Annealing optimization", "dual_annealing")
]
# Define one or the other:
LIMITS = "lmfit"


class DoModel:
    """
    Encapsulate all operations for modeling.  This class is instantiated
    only once, and the instance is used to interface to the remaining
    program and to perform the actual modeling.
    """

    def __init__(self):
        self.amw = None
        self.ya = None
        self.ya_list = None
        self.range = None
        self.range_list = None
        self.print_results = None
        self.draw_formatted = None
        self.exc_handler = None
        self.model = None

    # Local fitting functions ============================

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
        if log_mag:
            diff = np.log10(Z) - np.log10(self.model.model(w, params, **kwargs))
        else:
            diff = Z - self.model.model(w, params, **kwargs)
        diff *= weight
        # Flatten complex impedance into re/im adjacent floats
        residuals = diff.view('double')
        if LIMITS == "zfit":
            # Get mean-square float of entire difference array for weight scaling
            mean_sq = (residuals**2).mean()
            # Append penalties to residuals
            bounded = np.hstack(
                (residuals, self._bound_penalties(params, mean_sq)))
            return bounded
        else:
            return residuals

    def _bound_penalties(self, params, weight):
        """
        This function is only used when zfit_constants.LIMITS == 'zfit'.
        Return a list of numbers each of which increases rapidly when the min or max for
        a param is approached.  This represents boundary penalties as a parameter goes
        out of bounds.
        :param params:
        :param weight:
        :return: a list of N elements where N is the total number of min or max bounds
        in the set of params.  Each element is a number which increases rapidly when
        the min or max bound is approached.  Append penalty elements if params are out of bounds.
         """
        penalties = []
        # Custom limiting is done here:
        # This is an exponent which controls the abruptness of penalty increase as
        # a min or max limit is approached:
        PENALTY_WALL = 6
        full_penalty = 1e4 * weight
        for p in self.model.PARAMS:
            name = p['name']
            val = params[name].value
            # max and min must be >0 or None.
            # Use min/max limits from model.PARAMS here and set them to None
            # for lmfit (except for diff evo method).
            if p['max'] == None:
                max_pen = 0
            else:
                max_pen = full_penalty if val >= p['max'] else \
                    full_penalty * np.power(val/p['max'], PENALTY_WALL)
            if p['min'] == None:
                min_pen = 0
            else:
                min_pen = full_penalty if val <= p['min'] else \
                    full_penalty * np.power(p['min']/val, PENALTY_WALL)
            penalty = np.maximum(np.abs(max_pen), np.abs(min_pen))
            penalties.append(penalty)
        return penalties

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
        elif LIMITS == "zfit":
            # lmfit doesn't do the limiting
            return None
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

    def do_model_cli(self, range, ya, model, method_nr=0, do_norm_denorm=True):
        # Clear status and params label boxes

        # Import the model script, or reload it if already imported
        self.model = model
        self.model = reload(self.model)

        # Append data to form full set for modeling
        mag = np.array(ya.inputData[M])
        phase = np.array(np.radians(ya.inputData[P]))
        freq = np.array(range.xa["Hz"])
        load = np.array(range.load_array)

        # Radian frequency
        w = 2*np.pi*freq

        # Multiple data segments.
        # Create null weighting array, same size as m but full of 1's
        weight = mag.copy()
        weight.fill(1.0)

        # Complex impedance target
        z = mag*np.exp(1j*phase)

        # Instantiate clean class for lmfit fitter
        params = Parameters()
        params.clear()

        # Init list of name/value tuples
        values = []

        # Get selected fitting method
        method = METHODS[method_nr][1]

        # Do actual modeling.
        # Make working copy of PARAMS list from model
        param_list = list(self.model.PARAMS)

        # Adjust min and max if necessary
        for phase in param_list:
            phase["min"] = self._min_max_set(
                phase["min"], method, phase["init"] / 1e2)
            phase["max"] = self._min_max_set(
                phase["max"], method, phase["init"] * 1e2)

        if do_norm_denorm:
            # Normalize component, frequency, and impedance values
            # Determine frequency and Z scaling factors from initial values
            fsf, zsf = self._find_sf(param_list)
            # Normalize each component value, min, and max
            for phase in param_list:
                type = phase["name"][0].upper()
                norm = self._normalize(
                    phase["init"], phase["min"], phase["max"], type, fsf, zsf)
                phase["init"] = norm["init"]
                phase["min"] = norm["min"]
                phase["max"] = norm["max"]
            # Normalize frequency, target Z, and load
            w = w / fsf
            z = z / zsf
            # TODO: check what does gui in case of None
            #load = load / zsf
        else:
            fsf, zsf = 1.0, 1.0

        # Add modified params to lmfit Parameter class
        # .add converts min/max of None to -/+inf
        for phase in param_list:
            params.add(phase["name"], value=phase["init"],
                       vary=phase["vary"], min=phase["min"], max=phase["max"])

        # Perform weighted model optimization.
        # Errors will be caught and displayed by zfit_excepthook() in main window.
        kw_args = {"load": load, "fsf": fsf, "zsf": zsf}
        result = minimize(self._fcn2min, params, args=(
            w, z, weight), kws=kw_args, method=method)

        # Don't use params class after minimize -- some values are scrambled or changed.

        # Populate values[] with modeling results, denormalized if necessary
        for phase in param_list:
            name = phase["name"]
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
            #load = load * zsf

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
        kw_args = {"load": load, "fsf": 1.0, "zsf": 1.0}
        zfit = self.model.model(w, values_d, **kw_args)

        # Break into magnitude and degree phase
        magfit = np.abs(zfit)
        phasefit = np.angle(zfit, deg=True)

        # TODO: understand why a, b, i are 0!!!
        # # Split into segments as required, add to data list
        # a, b, i = 0, 0, 0
        # ya.modeledData[M] = magfit[a:b]
        # ya.modeledData[P] = phasefit[a:b]

        ya.modeledData[M] = magfit
        ya.modeledData[P] = phasefit

        self.ya = ya

        status = "Number of function calls: " + str(result.nfev)
        if result.aborted:
            status = "Process aborted"

        print(status)
        report_fit(result)


class Range:
    """
    This class is a holder for global data pertaining to the x or
    abscissa range.
    """
    def __init__(self):
        self.xa = {"Hz": None, "omega": None}
        self.load_array = None
        self.segment_str = ""
        self.segment_index = 0


class YAxes:
    """
    This class contains lists of items which differ between the mag, phase,
    and weight axes, and an index to select between them: M, P, or W.
    Items initialized with None are filled at run time_survey.
    """
    def __init__(self, index):
        # Variables containing the state of the axes
        self.index = index
        # Lists indexed by self.index
        self.ax = [None, None, None]
        self.label = ["Magnitude", "Phase", ""]
        self.axis_label = ["MAGNITUDE", "PHASE", ""]
        self.inputData = [None, None]   # mag, pha
        self.modeledData = [None, None]  # mag, pha
        self.drawnData = [None, None, None]     # drawn mag, pha, weight curves


if __name__ == "__main__":
    pass
