#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 17:51:58 2019

@author: tinivella
"""

import numpy as np
import plotly.figure_factory as ff
import plotly.graph_objects as go
from plotly.offline import plot
from plotly.subplots import make_subplots

# import matplotlib.pyplot as plt

# from scipy import interpolate
# from scipy import signal
# from scipy import integrate

# import mat4py
# import itertools

# from fft import calc_fourier_coefficients, calc_time_series, get_N_highest_peaks, calc_square_wave_DFT, BodePlotMatplotlib, TimeWaveformsMatplotlib

# from bokeh.plotting import figure, output_file, show, gridplot, output_notebook
# from bokeh.models import LinearAxis, HoverTool, Range1d
# from bokeh.palettes import Spectral6 as palette

# output_notebook()

# # colors has a list of colors which can be used in plots
# colors = itertools.cycle(palette)
# _tools_to_show = 'box_zoom,pan,save,hover,reset,tap,wheel_zoom'


class LogLogPlotDualAxisPlotly():
    def __init__(self, title=None):
        # Create figure with secondary y-axis
        self.fig = make_subplots(specs=[[{"secondary_y": True}]])
        self.fig.layout.template = 'plotly_dark'
        self.fig.layout.width = 1000
        self.fig.layout.height = 500
        # Add figure title
        self.fig.update_layout(title_text=title)
        self.fig.update_xaxes(title_text="FREQUENCY [Hz]")
        # Set y-axes titles
        self.fig.update_yaxes(title_text="'MAG [Ω]", secondary_y=False)
        self.fig.update_layout(xaxis_type="log", yaxis_type="log")
        self.fig.update_yaxes(title_text="PHASE [deg]", secondary_y=True)

    def add_trace(self, f, zabs=None, zdeg=None, label="None"):
        if zabs is not None and zdeg is not None:
            # Add traces
            self.fig.add_trace(
                go.Scatter(x=f, y=zabs, name=f"abs({label})"), secondary_y=False)
            self.fig.add_trace(
                go.Scatter(x=f, y=zdeg, name=f"deg({label})"), secondary_y=True)
        elif zabs is not None and zdeg is None:
            self.fig.add_trace(
                go.Scatter(x=f, y=zabs, name=f"abs({label})"), secondary_y=False)
        elif zabs is None and zdeg is not None:
            self.fig.add_trace(
                go.Scatter(x=f, y=zdeg, name=f"deg({label})"), secondary_y=True)
        else:
            raise AttributeError

    def show(self):
        self.fig.show()


# class LogLogPlotDualAxisBokeh():
#     def __init__(self, title=None):
#         # Create figure with secondary y-axis
#         self.p = figure(title=title, x_axis_type="log",
#                         y_axis_type="log",
#                         background_fill_color="#fafafa")

#         # # Add figure title
#         # self.fig.update_layout(title_text=title)
#         # self.fig.update_xaxes(title_text="FREQUENCY [Hz]")
#         # # Set y-axes titles
#         # self.fig.update_yaxes(title_text="'MAG [Ω]", secondary_y=False)
#         # self.fig.update_layout(xaxis_type="log", yaxis_type="log")
#         # self.fig.update_yaxes(title_text="PHASE [deg]", secondary_y=True)

#     def add_trace(self, f, zabs=None, zdeg=None, label="None"):
#         if zabs is not None and zdeg is not None:
#             # Add traces
#             self.p.line(x=f, y=zabs, legend_label=f"abs({label})")
#             self.p.extra_y_ranges={"foo": Range1d(start = -180, end = 180)}
#             self.p.line(x=f, y=zdeg, color = "red", y_range_name = "foo")
#             self.p.add_layout(LinearAxis(y_range_name="foo"), 'right')
#         elif zabs is not None and zdeg is None:
#             pass

#         elif zabs is None and zdeg is not None:
#             pass

#         else:
#             raise AttributeError

#     def show(self):
#         show(self.p)

