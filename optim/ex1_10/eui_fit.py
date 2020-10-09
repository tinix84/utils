from bokeh.plotting import figure
from bokeh.io import push_notebook, show, output_notebook
from bokeh.layouts import row

import numpy as np

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
    # mass density of wire conductor (kg/m**3)
    D.rowwc = 8890
    # conductivity of wire conductor (1/(Ohm*m))
    D.sigmawc = 5.959e7
    D.mu0 = np.pi*4e-7       # permeability of free space (H/m)
    D.irt = 10               # rated current (A)
    D.Lmn = 1e-3             # minimim inductance (H)
    D.Jmx = 7.5e6            # maximum current density (A/m**2)
    # maximum power dissipation at rated current (W)
    D.Pmx = 1
    D.Mmx = 1  # maximum mass (kg)
    
    return D


class MyProblem(Problem):
    def __init__(self):
        super().__init__()



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

    # set up genetic algorithm parameters--------------------------------------
    nobj=2                 # number of objectives
    oopt=0                 # objective to optimize
    ngen=100              # number of generations
    npop=100              # population size
    
    # GAP=gapdefault(nobj,oopt,npop,ngen)

    GAP=Struct()
    # set up search space------------------------------------------------------
    #               N   ds   ws   wc   lc   g            
    GAP.gd_min=np.array([1,  1e-3, 1e-3, 1e-3, 1e-3, 1e-5])
    GAP.gd_max=np.array([1e3, 1e-1, 1e-1, 1e-1, 1e-1, 1e-2])
    GAP.gd_type=np.array([3, 3, 3, 3, 3, 3])
    GAP.gd_cid=np.array([1, 1, 1, 1, 1, 1])

    problem = MyProblem()
    problem.n_var = len(GAP.gd_min)
    problem.n_obj = nobj
    problem.n_constr = 2
    problem.xl = GAP.gd_min
    problem.xu = GAP.gd_max
    problem.elementwise_evaluation=True

    algorithm = NSGA2(
        pop_size=npop,
        eliminate_duplicates=True
    )

    res = minimize(problem, algorithm, ("n_gen", ngen), verbose=True)

    # # conduct the optimization-------------------------------------------------
    # # [fP,GAS,bi,bf]= gaoptimize(@eui_fit,GAP,D)
    # eui_fit(x=GAP, D=D, fn=None)

    # save results-------------------------------------------------------------
    # savemoresults
    return res

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
    #  D = structure of design paramters
    #   D.kpf = packing factor
    #   D.rowmc = mass density of core material(kg/m ** 3)
    #   D.Bmx = flux limit of core material(kg/m ** 3)
    #   D.rowwc = mass density of wire conductor(kg/m ** 3)
    #   D.sigmawc = conductivity of wire conductor(1/(Ohm*m))
    #   D.mu0 = permeability of free space(H/m)
    #   D.irt = rated current(A)
    #   D.Lmn = minimim inductance(H)
    #   D.Jmx = maximum current density(A/m ** 2)
    #   D.Pmx = maximum power dissipation at rated current(W)
    #   D.Mmx = maximum mass(kg)
    # fn = figure number for design report
    #
    # Output:
    # f = fitness
    #  f(1) = recprical of mass(if all constriants met)(1/kg)
    #  f(2) = recprical of power loss at rated current(if all constraints
                                                    # met)(1/W)
    #
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
    # xc = x-coordinates of objects used to draw inductor(m)
    # yc = y-coordinates of objects used to draw inductor(m)
    # dgray = color code for dark gray
    # lgray = color code of light gray
    #

    # break paramameter vector into components
    # N = round(x(1))
    # force number of turns to be integer
    N =  x[0]
    ds = x[1]
    ws = x[2]
    wc = x[3]
    lc = x[4]
    g =  x[5]

    # compute mass
    M = 2.0*(2.0*wc+ws+ds)*lc*wc*D.rowmc+(2*lc+2*wc+np.pi*ds)*ds*ws*D.kpf*D.rowwc
    # compute loss at rated current
    Prt = (2*lc+2*wc+np.pi*ds)*(N*D.irt) ** 2/(ds*ws*D.kpf*D.sigmawc)
    # compute inductance
    L = D.mu0*lc*wc*N ** 2/(2*g)
    # compute the flux density
    Brt = D.mu0*N*D.irt/(2*g)
    # current density
    Jrt = N*D.irt/(ws*ds*D.kpf)
    # energy metric constraint
    c1 = gte(L, D.Lmn)
    # flux density constraint
    c2 = lte(Brt, D.Bmx)
    # current density constraint
    c3 = lte(Jrt, D.Jmx)
    # power loss constraint
    c4 = lte(Prt, D.Pmx)
    # mass constraint
    c5 = lte(M, D.Mmx)
    # aggrecate constraint
    c = (c1+c2+c3+c4+c5)/5

    # # fitness
    # if (c < 1):
    #     f = 1e-3*(c-1)*np.array([1, 1])
    # else:
    #     f = np.array([M, Prt])
    f = np.array([M, Prt])
    g = np.array([D.Lmn - L, Brt - D.Bmx,
                    Jrt - D.Jmx, Prt - D.Pmx,
                    M-D.Mmx])  
    return f, g

def draw_inductor():
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

    N = x[0]
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


if __name__ == "__main__":


    res = eui_modesign()
    
    TOOLS = "hover,pan,wheel_zoom,zoom_in,zoom_out,box_zoom,undo,redo,reset,tap,save,box_select,poly_select,lasso_select,"
    
    TOOLTIPS = [
        ("index", "$index"),
        ("(x,y)", "($x, $y)"),
        ("foo", "@res.F")]
    
    p = figure(tools=TOOLS, tooltips=None, plot_width=600, plot_height=450,
               title="inductor optimization")
    p.scatter(x=res.F[:, 0], y=res.F[:, 1])
    p.xaxis.axis_label = 'mass'
    p.yaxis.axis_label = 'losses'
    show(p)

    print_design(res.X[26], D)
