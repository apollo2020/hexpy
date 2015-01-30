# created by: David Mangold, mangoldd
# last updated: 20140630
# purpose: tools for working with regular hexagons
# input parameters: NA, module
# dev notes:

import os
from geom import *

import sys
if sys.version_info[0] < 3:
    import Tkinter as Tk
else:
    import tkinter as Tk

import matplotlib
matplotlib.use('TkAgg', warn=False)
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure


class Plotter(object):
    """Tools for visually showing hexpy.geom instance geometries."""

    def __init__(self):
        """Setup Plotter to use TkInter."""

        self.fig = Figure(figsize=(5, 4), dpi=100)
        self.fig_sub = self.fig.add_subplot(1, 1, 1)

    def build_canvas(self):
        """Create Tk canvas for current figure."""

    def plot_all(self, hex_geom_list=[]):
        """Plot all object in hex_geom_list using self.plot."""

        for hex_geom_obj in hex_geom_list:
            self.plot(hex_geom_obj)

    def plot(self, hex_geom_obj=None):
        """Plot hexpy.geom instances using pyplot."""

        if type(hex_geom_obj) is Extent:

            x_vals = [v.x for v in hex_geom_obj.vertices] + [hex_geom_obj.vertices[0].x]
            y_vals = [v.y for v in hex_geom_obj.vertices] + [hex_geom_obj.vertices[0].y]
            self.fig_sub.plot(x_vals, y_vals)
            return

        elif type(hex_geom_obj) is HexPoint:

            x_vals = hex_geom_obj.x
            y_vals = hex_geom_obj.y
            self.fig_sub.plot(x_vals, y_vals, 'ro', markersize=10)
            return

        elif type(hex_geom_obj) is HexCell:

            x_vals = [v.x for v in hex_geom_obj.vertices] + [hex_geom_obj.vertices[0].x]
            y_vals = [v.y for v in hex_geom_obj.vertices] + [hex_geom_obj.vertices[0].y]
            self.fig_sub.plot(x_vals, y_vals)
            return

        elif type(hex_geom_obj) is HexGroup:

            for cell in hex_geom_obj.cells:

                x_vals = [v.x for v in cell.vertices] + [cell.vertices[0].x]
                y_vals = [v.y for v in cell.vertices] + [cell.vertices[0].y]
                self.fig_sub.plot(x_vals, y_vals)
            return

        else:

            print("%s is not valid. Use a HexPoint, HexCell, or HexGroup instance instead." % hex_geom_obj)
            return

    def show(self):
        """Show the current plot."""

        root = Tk.Tk()
        root.title("HexPy Plotter")

        canvas = FigureCanvasTkAgg(self.fig, master=root)
        canvas.show()
        canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

        toolbar = NavigationToolbar2TkAgg(canvas, root)
        toolbar.update()
        canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

        Tk.mainloop()

    def save(self, path, name):
        """Save the current plot as an image."""

        root = Tk.Tk()
        FigureCanvasTkAgg(self.fig, master=root)
        filepath = os.path.join(path, name)
        self.fig.savefig(filepath)
        root.destroy()

    def clear(self):
        """Clear all figure contents and create a new default sub-plot."""

        self.fig.clf()
        self.fig_sub = self.fig.add_subplot(1, 1, 1)
