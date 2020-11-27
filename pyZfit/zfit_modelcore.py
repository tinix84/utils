import numpy as np
from importlib import import_module, reload
from zfit_constants import *
from lmfit import minimize, Parameters
from scipy.stats import gmean   # geometric mean
import csv


class CliInterface():
    def __init__(self):
        self.checkBoxLogMag
        

class DoModel:
    """
    Encapsulate all operations for modeling.  This class is instantiated
    only once, and the instance is used to interface to the remaining
    program and to perform the actual modeling.
    """

    def __init__(self, cli=False):
        if cli is False:
            # Interface to the main program which instantiates this class
            self.amw = None
            self.ya = None
            self.ya_list = None
            self.range = None
            self.range_list = None
            self.print_results = None
            self.draw_formatted = None
            self.exc_handler = None
        else:
            # Interface to the main program which instantiates this class
            self.amw = CliInterface()
            self.ya = None
            self.ya_list = None
            self.range = None
            self.range_list = None
            self.print_results = None
            self.draw_formatted = None
            self.exc_handler = None

    # Local fitting functions ============================

    def _fcn2min(self, params, w, Z, weight, **kwargs):
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
        if self.amw.checkBoxLogMag.isChecked():
            diff = np.log10(
                Z) - np.log10(self.model.model(w, params, **kwargs))
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

    def _prog_bar_update(self, x1, x2, x3, *x4, **x5):
        """
        Update the progress bar in the main app.
        (Could look for abort here too.)
        :param xn: Dummy args to match actual sent by minimize() callback
        :return: nothing
        """
        self.amw.prog_bar_tick += 1
        if self.amw.prog_bar_tick >= 100:
            self.amw.prog_bar_tick = 0
            next = (self.amw.progressBar.value() + 1) % 100
            self.amw.progressBar.setValue(next)

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

    def do_model(self):
        # Clear status and params label boxes
        self.amw.labelParams.setText("")
        self.amw.labelParams.repaint()
        self.amw.labelStatus.setText("Modeling...")
        self.amw.labelStatus.repaint()

        # Import the model script, or reload it if already imported
        self.model = import_module("Models." + self.amw.lineEditModel.text())
        self.model = reload(self.model)

        # Clear any previous modeling, M and P axes
        for line in self.ya.ax[M].get_lines() + self.ya.ax[P].get_lines():
            if line.get_label() == "modeledZPlot":
                line.remove()

        # Create local concatenated arrays for freq, mag, phase, and load.
        # Overwrites data in ya class
        m = np.array([])
        p = np.array([])
        f = np.array([])
        l = np.array([])
        for y, r in zip(self.ya_list, self.range_list):
            # Copy data from lists to working objects
            self.ya.data_unbundle(y)
            self.range.data_unbundle(r)
            # Append data to form full set for modeling
            m = np.append(m, self.ya.inputData[M])
            p = np.append(p, np.radians(self.ya.inputData[P]))
            f = np.append(f, self.range.xa["Hz"])
            l = np.append(l, self.range.load_array)

        # Radian frequency
        w = 2*np.pi*f

        # Use drawn curves if they exist
        if self.ya.drawnData[M] is not None:
            # Drawn data exists for magnitude, use it instead
            m = self.ya.drawnData[M]
        if self.ya.drawnData[P] is not None:
            # Drawn data exists for phase, use it instead
            p = np.radians(self.ya.drawnData[P])

        if len(self.ya_list) > 1:
            # Multiple data segments.
            # Create null weighting array, same size as m but full of 1's
            weight = m.copy()
            weight.fill(1.0)
        else:
            weight = self.ya.drawnData[W]

        # Complex impedance target
        z = m*np.exp(1j*p)

        # Instantiate clean class for lmfit fitter
        params = Parameters()
        params.clear()

        # Init list of name/value tuples
        values = []

        # Get selected fitting method
        method = METHODS[self.amw.comboBoxMethod.currentIndex()][1]

        if self.amw.checkBoxLocked.isChecked():
            # Read last saved or edited params data (denormalized)
            with open(PARAM_FILE, mode='r', encoding='utf-8', newline='') as f:
                reader = csv.reader(f)
                next(f)      # skip header line
                for line in reader:
                    v = (line[0], float(line[1]))
                    # Build a list of name/value tuples
                    values.append(v)
        else:
            # Do actual modeling.
            # Make working copy of PARAMS list from model
            param_list = list(self.model.PARAMS)

            # Adjust min and max if necessary
            for p in param_list:
                p["min"] = self._min_max_set(p["min"], method, p["init"] / 1e2)
                p["max"] = self._min_max_set(p["max"], method, p["init"] * 1e2)

            if self.amw.do_norm_denorm:
                # Normalize component, frequency, and impedance values
                # Determine frequency and Z scaling factors from initial values
                fsf, zsf = self._find_sf(param_list)
                # Normalize each component value, min, and max
                for p in param_list:
                    type = p["name"][0].upper()
                    norm = self._normalize(
                        p["init"], p["min"], p["max"], type, fsf, zsf)
                    p["init"] = norm["init"]
                    p["min"] = norm["min"]
                    p["max"] = norm["max"]
                # Normalize frequency, target Z, and load
                w = w / fsf
                z = z / zsf
                l = l / zsf
            else:
                fsf, zsf = 1.0, 1.0

            # Add modified params to lmfit Parameter class
            # .add converts min/max of None to -/+inf
            for p in param_list:
                params.add(p["name"], value=p["init"],
                           vary=p["vary"], min=p["min"], max=p["max"])

            # Perform weighted model optimization.
            # Errors will be caught and displayed by zfit_excepthook() in main window.
            kw_args = {"load": l, "fsf": fsf, "zsf": zsf}
            result = minimize(self._fcn2min, params, args=(w, z, weight), kws=kw_args, method=method,
                              iter_cb=self._prog_bar_update)

            # Don't use params class after minimize -- some values are scrambled or changed.

            # Populate values[] with modeling results, denormalized if necessary
            for p in param_list:
                name = p["name"]
                val = result.params[name].value
                if self.amw.do_norm_denorm:
                    comp_type = name[0]
                    val = self._denormalize(val, comp_type, fsf, zsf)
                v = (name, val)
                values.append(v)

            if self.amw.do_norm_denorm:
                # Denormalize frequency, target Z, and load
                w = w * fsf
                z = z * zsf
                l = l * zsf

            # Write denormalized modeling results to file
            with open(PARAM_FILE, mode='w', encoding='utf-8') as f:
                print('name, value', file=f)
                for p in values:
                    print('{}, {}'.format(p[0], p[1]), file=f)

            self.amw.progressBar.setValue(0)

        # Convert list of tuples to a single dict for the model, to be compatible
        # with the way minimize() uses the model
        values_d = {p[0]: p[1] for p in values}

        # Get complex impedance of model using modeled or locked parameters
        # Use denormalized values
        kw_args = {"load": l, "fsf": 1.0, "zsf": 1.0}
        zfit = self.model.model(w, values_d, **kw_args)

        # Break into magnitude and degree phase
        magfit = np.abs(zfit)
        phasefit = np.angle(zfit, deg=True)
        # Split into segments as required, add to data list
        a, b, i = 0, 0, 0
        for y in self.ya_list:
            self.ya.data_unbundle(y)
            seg_length = len(self.ya.inputData[M])
            b += seg_length
            self.ya.modeledData[M] = magfit[a:b]
            self.ya.modeledData[P] = phasefit[a:b]
            a += seg_length
            self.ya_list[i] = self.ya.data_bundle()
            i += 1

        # Refresh working classes with currently indexed segment data
        self.ya.data_unbundle(self.ya_list[self.range.segment_index])
        self.range.data_unbundle(self.range_list[self.range.segment_index])

        # Add to plot
        self.ya.ax[M].plot(self.range.xa["Hz"], self.ya.modeledData[M], self.ya.modeledLinePlot[M],
                           ls=self.ya.modeledLineStyle[M], lw=1, label="modeledZPlot")
        self.ya.ax[P].plot(self.range.xa["Hz"], self.ya.modeledData[P], self.ya.modeledLinePlot[P],
                           ls=self.ya.modeledLineStyle[P], lw=1, label="modeledZPlot")

        # Update results text box
        self.print_results(values)

        if self.amw.checkBoxLocked.isChecked():
            self.amw.labelStatus.setText("")
        else:
            # Append "written to" text to results box
            outstr = self.amw.labelParams.text()
            outstr += "<br>Written to<br>" + PARAM_FILE
            self.amw.labelParams.setText(outstr)
            # Print optimization info
            status = "Number of function calls: " + str(result.nfev) + "<br>"
            if result.aborted:
                status = RICH_TEXT_RED + "Process aborted:<br>"
            #status += result.lmdif_message
            self.amw.labelStatus.setText(status)

        self.draw_formatted()
