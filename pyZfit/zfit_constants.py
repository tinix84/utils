# indices for weight, magnitude, and phase lists
M, P, W = 0, 1, 2
# Set max for model weighting.  Minimum is 1.0
MAX_MODEL_WEIGHT = 100.0
# Rich text control for red font
RICH_TEXT_RED = "<FONT COLOR='red'>"
# Color definitions
# (see http://www.w3schools.com/tags/ref_colorpicker.asp )
M_COLOR = "#6f93ff"  # magnitude data
P_COLOR = "#ffcc00"  # phase data
DM_COLOR = "#3366ff" # drawn magnitude
DP_COLOR = "#b38f00" # drawn phase
DW_COLOR = "#33cc33" # drawn weight
MM_COLOR = "000000"  # modeled magnitude lines
MP_COLOR = "000000"  # modeled phase lines
# Line styles for modeled data
MM_STYLE = "-"
MP_STYLE = "--"
# Boolean for whether to allow negative param results
ALLOW_NEG = True
# File name for parameter save/restore
PARAM_FILE = "params.csv"
# Path and file name for preferred text editor for editing model scripts
# TODO: make this usable from os.startfile() -- not successful so far
EDITOR = "notepad++.exe"
EDIT_PATH = "C:\\Program Files\\Notepad++\\"
# List of fitting methods available to lmfit
# (Dogleg and Newton CG are not included since they require a Jacobian
# to be supplied)
METHODS = [
    ("Levenberg-Marquardt",         "leastsq"),
    ("Nelder-Mead",                 "nelder"),
    ("L-BFGS-B",                    "lbfgsb"),
    ("Powell",                      "powell"),
    ("Conjugate Gradient",          "cg"),
    ("COBYLA",                      "cobyla"),
    ("Truncated Newton",            "tnc"),
    ("Sequential Linear Squares",   "slsqp"),
    ("Differential Evolution",      "differential_evolution")
]
# Define one or the other:
LIMITS = "lmfit"
# LIMITS = "zfit"