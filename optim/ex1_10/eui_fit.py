from bokeh.plotting import figure
from bokeh.io import push_notebook, show, output_notebook
from bokeh.layouts import row

import numpy as np
import pandas as pd

from pymoo.algorithms.nsga2 import NSGA2
from pymoo.model.problem import Problem
from pymoo.optimize import minimize
from pymoo.visualization.scatter import Scatter
from pymoo.factory import get_sampling, get_crossover, get_mutation
from pymoo.visualization.pcp import PCP

class Struct(object):
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

def input_design_parameters():
    # design information-------------------------------------------------------
    D = Struct()             # D = structure of design paramters
    D.kpf = 0.7              # packing factor
    D.rowmc = 4683.7         # mass density of core material (kg/m**3)
    D.Bmx = 0.6171           # flux limit of core material (kg/m**3)
    D.rowwc = 8890     # mass density of wire conductor (kg/m**3)
    D.sigmawc = 5.959e7  # conductivity of wire conductor (1/(Ohm*m))
    D.mu0 = np.pi*4e-7       # permeability of free space (H/m)
    D.irt = 10               # rated current (A)
    D.Lmn = 1e-3             # minimim inductance (H)
    D.Jmx = 7.5e6            # maximum current density (A/m**2)
    D.Pmx = 1  # maximum power dissipation at rated current (W)
    D.Mmx = 1  # maximum mass (kg)
    
    return D

class MyProblem(Problem):
    def __init__(self):
        super().__init__()
        self.D = input_design_parameters()

    def _evaluate(self, x, out, *args, **kwargs):
        f, g = eui_fit(x, self.D, None)
        out["F"] = f
        out["G"] = g

def eui_modesign():
    """
    eui_modesign performs a multi-objective optimal design of a ui core 
                  inductor based on an elementary analysis.  Used in 
                  section 10.1 of "Power Magnetic Devices: A Multi-Objective
                  Deign Approach" by S.D. Sudhoff
    """

    # setup design space
    #               N   ds   ws   wc   lc   g
    GAP = Struct()
    GAP.gd_min = np.array([1,  1e-3, 1e-3, 1e-3, 1e-3, 1e-5])
    GAP.gd_max = np.array([1e3, 1e-1, 1e-1, 1e-1, 1e-1, 1e-2])


    # setup genetic algorithm parameters--------------------------------------
    nobj=2                 # number of objectives
    ngen=100              # number of generations
    npop = 100  # population size
    
    problem = MyProblem()
    problem.n_var = len(GAP.gd_min)
    problem.n_obj = nobj
    problem.n_constr = 2
    problem.xl = GAP.gd_min
    problem.xu = GAP.gd_max
    problem.elementwise_evaluation = True

    algorithm = NSGA2(
        pop_size=npop,
        eliminate_duplicates=True
    )
    
    # conduct the optimization-------------------------------------------------
    res = minimize(problem, algorithm, ("n_gen", ngen), verbose=True)

    # save results-------------------------------------------------------------
    return res

def calc_inductor(x, D):
    # break paramameter vector into components
    # Internal:
    # N = number of turns
    # ds = slot depth(m)
    # ws = slot width(m)
    # wc = core width(m)
    # lc = core length(m)
    # g = air gap(m)
    # M = mass(kg)
    # Prt = loss at rated power(W)
    # L = inductance(H)
    # Brt = flux denity at rated current(T)
    # Jrt = current density at rated current(A/m ** 2)

    N = round(x[0]) # force number of turns to be integer
    # N = x[0]
    ds = x[1]
    ws = x[2]
    wc = x[3]
    lc = x[4]
    g = x[5]

    # compute mass
    M = 2.0*(2.0*wc+ws+ds)*lc*wc*D.rowmc + \
        (2*lc+2*wc+np.pi*ds)*ds*ws*D.kpf*D.rowwc
    # compute loss at rated current
    Prt = (2*lc+2*wc+np.pi*ds)*(N*D.irt) ** 2/(ds*ws*D.kpf*D.sigmawc)
    # compute inductance
    L = D.mu0*lc*wc*N ** 2/(2*g)
    # compute the flux density
    Brt = D.mu0*N*D.irt/(2*g)
    # current density
    Jrt = N * D.irt / (ws * ds * D.kpf)
    
    return M, Prt, L, Brt, Jrt

def eui_fit(x, D, fn):
    # eui_fit is a fitness funtion for an elementary design of a ui core
    #         inductor, as described in section 1.10 of "Power
    #         Magnetic Devices: A Multi-Objective Design Approach"
    #         by S.D. Sudhoff
    #
    # Call:
    # f = eui_fit(x, D)
    # f = eui_fit(x, D, fn)
    #
    # Input:
    # x = parameter vector
    #  x(1) = desired number of turns(as real number)
    #  x(2) = slot depth(m)
    #  x(3) = slot width(m)
    #  x(4) = core width(m)
    #  x(5) = core length(m)
    #  x(6) = air gap(m)

    M, Prt, L, Brt, Jrt = calc_inductor(x, D)
    
    # optimization fitness functions
    f1 = M
    f2 = Prt

    # optimization constrains
    # energy metric constraint
    g1 = D.Lmn - L
    # flux density constraint
    g2 = Brt - D.Bmx
    # current density constraint
    g3 = Jrt - D.Jmx
    # power loss constraint
    g4 = Prt - D.Pmx
    # mass constraint
    g5 = M-D.Mmx


    # vectorization
    f = np.array([f1, f2])
    g = np.array([g1, g2, g3, g4, g5])  
    return f, g

def draw_inductor():
    # xc = x-coordinates of objects used to draw inductor(m)
    # yc = y-coordinates of objects used to draw inductor(m)
    # dgray = color code for dark gray
    # lgray = color code of light gray
    # # draw cross section of inductor
    # figure(fn)
    # dgray = [0.5 0.5 0.5]
    # # color code for dark grey
    # lgray = [0.9 0.9 0.9]
    # xc = [-(ws/2+wc)(ws/2+wc)(ws/2+wc) - (ws/2+wc}
    # yc = [0         0 - (ds+wc) - (ds+wc};
    # fill(xc, yc, dgray)
    # hold on
    # yc = [(wc+g)(wc+g)      g            g];
    # fill(xc, yc, dgray)
    # xc = [-ws/2 ws/2 ws/2 - ws/2];
    # yc = [0    0 - ds - ds];
    # fill(xc, yc, lgray)
    # yc = [0    0 - ds - ds]-(wc+ds);
    # fill(xc, yc, lgray)
    # hold off
    # axis equal
    pass

def print_design(x, D):
    """Print in the console the design parameters

    Args:
        x ([type]): optimization variables x=[N, ds, ws, wc, lc, g]
    """

    N = round(x[0])
    ds = x[1]
    ws = x[2]
    wc = x[3]
    lc = x[4]
    g = x[5]

    # compute mass
    M = 2.0*(2.0*wc+ws+ds)*lc*wc*D.rowmc + \
        (2*lc+2*wc+np.pi*ds)*ds*ws*D.kpf*D.rowwc
    # compute loss at rated current
    Prt = (2*lc+2*wc+np.pi*ds)*(N*D.irt) ** 2/(ds*ws*D.kpf*D.sigmawc)
    # compute inductance
    L = D.mu0*lc*wc*N ** 2/(2*g)
    # compute the flux density
    Brt = D.mu0*N*D.irt/(2*g)
    # current density
    Jrt = N*D.irt/(ws*ds*D.kpf)
    print('Design Data')
    print(f'Turns = {N}')
    print(f'Slot depth (m) = {ds}')
    print(f'Slot width (m) = {ws}')
    print(f'Core width (m) = {wc}')
    print(f'Core length (m) = {lc}')
    print(f'Air gap (m) = {g}')
    print(' ')
    print('Design Metrics')
    print(f'Mass (kg) = {M}')
    print(f'Loss at rated current (W) = {Prt}')
    print(' ')
    print('Constrained Quantities')
    print(f'Inductance (H) = {L}')
    print(f'Flux Density at Rated Current (T) = {Brt}')
    print(f'Current Density Rated Current (A/m**2) = {Jrt}')


def lte(x, xmx):
    """LTE      less-than-or-equal to function
        c = lte(x, xmx)

        Args:
            x ([type]): a quantitity
            xmn ([type]): minimum allowed value

        Returns:
            [type]: = c = constraint variable - 1 if x <= xmx, 0 < c < 1 if x > xmx
        """
    if (x <= xmx):
        c = 1
    else:
        c = 1/(1+x-xmx)

    return c

def gte(x, xmn):
    """GTE        greater-than-or-equal to function
    c = gte(x, xmn)

    Args:
        x ([type]): a quantitity
        xmn ([type]): minimum allowed value

    Returns:
        [type]: = constraint variable - 1 if x >= xmn, 0 <c<1 if x<xmn
    """
    if (x >= xmn):
        c = 1
    else:
        c = 1/(1+xmn-x)

    return c

def get_design(x, D):
    des = dict()
    """Print in the console the design parameters

    Args:
        x ([type]): optimization variables x=[N, ds, ws, wc, lc, g]
    """

    des['N'] = N = round(x[0])
    des['ds'] = ds = x[1]
    des['ws'] = ws = x[2]
    des['wc'] = wc = x[3]
    des['lc'] = lc = x[4]
    des['g'] = g = x[5]

    # compute mass
    des['M'] = 2.0*(2.0*wc+ws+ds)*lc*wc*D.rowmc + \
        (2*lc+2*wc+np.pi*ds)*ds*ws*D.kpf*D.rowwc
    # compute loss at rated current
    des['Prt'] = (2*lc+2*wc+np.pi*ds)*(N*D.irt) ** 2/(ds*ws*D.kpf*D.sigmawc)
    # compute inductance
    des['L'] = D.mu0*lc*wc*N ** 2/(2*g)
    # compute the flux density
    des['Brt'] = D.mu0*N*D.irt/(2*g)
    # current density
    des['Jrt'] = N*D.irt/(ws*ds*D.kpf)

    return pd.Series(des)

if __name__ == "__main__":

    from bokeh.models import ColumnDataSource, NumeralTickFormatter, HoverTool

    res = eui_modesign()

    TOOLS = "hover,pan,wheel_zoom,zoom_in,zoom_out,box_zoom,undo,redo,reset,tap,save,box_select,poly_select,lasso_select,"

    TOOLTIPS = [
        ("index", "$index"),
        ('L', "@L"),
        ('ds', "@ds"),
        ('ws', "@ws"),
        ('lc', "@lc"),
        ('wc', "@wc"),
        ('B', "@Brt"),
        ('Nturns', "@N"),
        ('lgap', "@g")
    ]

    sol = pd.DataFrame()
    for v in res.X:
        row = get_design(v, input_design_parameters())
        sol = sol.append(row, ignore_index=True)

    source = ColumnDataSource(sol)

    p = figure(tools=TOOLS, tooltips=TOOLTIPS, plot_width=600, plot_height=450,
               title="inductor optimization")

    p.circle(x='M', y='Prt', size=5, source=source)
    #p.scatter(x=sol.M, y=sol.Prt)
    p.xaxis.axis_label = 'mass'
    p.yaxis.axis_label = 'losses'
    # Add the HoverTool to the figure
    #p.add_tools(HoverTool(tooltips=tooltips))
    show(p)

    print_design(res.X[26], input_design_parameters())
