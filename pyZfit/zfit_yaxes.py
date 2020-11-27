from zfit_constants import *
import copy
from matplotlib.ticker import EngFormatter


class YAxes:
    """
    This class contains lists of items which differ between the mag, phase,
    and weight axes, and an index to select between them: M, P, or W.
    Items initialized with None are filled at run time_survey.
    """

    def __init__(self, index):
        # Variables containing the state of the axes
        self.index = index
        self.fig = None
        # Lists indexed by self.index
        self.ax = [None, None, None]
        self.label = ["Magnitude", "Phase", ""]
        self.axis_label = ["MAGNITUDE", "PHASE", ""]
        self.dataColor = [M_COLOR, P_COLOR, DW_COLOR]
        self.drawnLinePlot = [DM_COLOR, DP_COLOR, DW_COLOR]
        self.modeledLinePlot = [MM_COLOR, MP_COLOR, ""]
        self.modeledLineStyle = [MM_STYLE, MP_STYLE, ""]
        self.inputData = [None, None]   # mag, pha
        self.modeledData = [None, None] # mag, pha
        self.drawnData = [None, None, None]     # drawn mag, pha, weight curves

    def select(self, index):
        # Mouse events will only register on the top plot, so
        # set current axes for active one, also set z order and make the
        # canvas (patch) of the top plot not visible so it doesn't block the other.
        # Set only the lowest patch visible.  Set index to match label.
        self.index = index
        if self.index == M:
            # Magnitude
            self.fig.sca(self.ax[M])
            self.ax[M].set_zorder(1)
            self.ax[P].set_zorder(.5)
            self.ax[W].set_zorder(.1)
            self.ax[M].patch.set_visible(False)
            self.ax[P].patch.set_visible(False)
            self.ax[W].patch.set_visible(True)
        elif self.index == P:
            # Phase
            self.fig.sca(self.ax[P])
            self.ax[M].set_zorder(.5)
            self.ax[P].set_zorder(1)
            self.ax[W].set_zorder(.1)
            self.ax[M].patch.set_visible(False)
            self.ax[P].patch.set_visible(False)
            self.ax[W].patch.set_visible(True)
        elif self.index == W:
            # Weight
            self.fig.sca(self.ax[W])
            self.ax[M].set_zorder(.5)
            self.ax[P].set_zorder(.1)
            self.ax[W].set_zorder(1)
            self.ax[M].patch.set_visible(False)
            self.ax[P].patch.set_visible(True)
            self.ax[W].patch.set_visible(False)

    def data_bundle(self):
        # Return a list of copies of the current data
        iMag = copy.copy(self.inputData[M])
        iPha = copy.copy(self.inputData[P])
        mMag = copy.copy(self.modeledData[M])
        mPha = copy.copy(self.modeledData[P])
        return [iMag, iPha, mMag, mPha]

    def data_unbundle(self, bundle):
        # Complement of function above
        self.inputData[M] = bundle[0]
        self.inputData[P] = bundle[1]
        self.modeledData[M] = bundle[2]
        self.modeledData[P] = bundle[3]
