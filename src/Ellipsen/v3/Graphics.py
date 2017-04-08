########################################################################################################################
#                                                                                                                      #
#                                    Software Projekt: Frechet Distanz                                                 #
#                                    Teilgebiet: Ellipsen-Alg. einer Zelle                                             #
#                                    Erstellt: WS 16/17 FU Berlin                                                      #
#                                                                                                                      #
#                                    Team: Josephine Mertens, Jana Kirschner,                                          #
#                                    Alexander Korzech, Fabian Kovacs, Alexander                                       #
#                                    Timme, Kilian Kraatz & Anton Begehr                                               #
#                                                                                                                      #
########################################################################################################################

# -*- coding: utf-8 -*-

from Algorithm import *

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from numbers import Number
from matplotlib.widgets import Slider

class ActivePlot:

    sample = []
    plot_3d = True
    ax_3d = 0
    ax_2d = 0
    axis_2d = 0
    bounds_l = 0
    pad_1 = 0
    curr_val = 0

    def __init__(self):
        self.plt = plt
    
    def vectors_to_xy(self, vectors: [Vector]) -> ( [float], [float]):  # converts array of vectors to x- & y-coordinate arrays
        x = []
        y = []
    
        for vector in vectors:
            if isinstance(vector, Vector):
                x.append(vector.x)
                y.append(vector.y)
            else:
                print("Error: not a Vector: " + str(vector))
    
        return x, y
    
    
    def xy_to_vectors(self, xs: [float], ys: [float]) -> [Vector]:  # converts x- & y-coordinate arrays to arrays of vectors
        vectors = []
    
        for i in range(min(len(xs), len(ys))):
            x = xs[i]
            y = ys[i]
            if isinstance(x, Number) and isinstance(y, Number):
                vectors.append(Vector(x, y))
            else:
                print("Error: either is not a valid float: x:" + str(x) + " y:" + str(y))
    
        return vectors
        
    def update(self,val, plot_heatmap, plot_3d, sample):
        self.curr_val = val
        
        self.plotTraversals(val, plot_heatmap, plot_3d, sample)
        self.plt.show()
        
    def plotTraversals(self, val, plot_heatmap, plot_3d, sample):
        
        # plot 3d
        if plot_3d:
           heatmap = sample["heatmap"]
           x = heatmap[0]
           y = heatmap[1]
           z = heatmap[2]            
           self.ax_3d.set_zlim([self.bounds_l[0] - self.pad_1, self.bounds_l[1] + self.pad_1])
            
           surf = self.ax_3d.plot_surface(x, y, z, cmap=cm.coolwarm, rstride=1, cstride=1,
                                                linewidth=0, antialiased=False)
        
        if plot_heatmap:
            heatmap = sample["heatmap"]
            x = heatmap[0]
            y = heatmap[1]
            z = heatmap[2]
            surf = self.ax_2d.pcolor(x, y, z, cmap=cm.coolwarm)
        
        
    
    
    def sample_to_matplotlib(self, sample, plot_borders: bool = True, plot_ellipsis: bool = True, plot_heatmap: bool = True,
                            plot_traversals: bool = True, plot_axis: bool = False, plot_l_lines: bool = False,
                            plot_3d: bool = True, show_legend: bool = True, show_colorbar: bool = True,
                            plot_input: bool = True, show_labels: bool = True, plot_critical_traversals: bool = False,
                            plot_cross_sections: bool = False):
        # plot sample with matplotlib
        
        self.sample = sample
        self.plot_3d = plot_3d
        padding = 0.04
        curr_val = 0
    
        # plot cross-sections
        if plot_cross_sections:
            self.bounds_l = sample["bounds-l"]
            fig_both = plt.figure(figsize=plt.figaspect(0.5))
            ax_hor = fig_both.add_subplot(2, 1, 1, aspect=1, ylim=self.bounds_l, xlabel="p", ylabel="ε")
            ax_ver = fig_both.add_subplot(2, 1, 2, aspect=1, ylim=self.bounds_l, xlabel="q", ylabel="ε")
            cross_sections_hor = sample["cross-sections-hor"]
            cross_sections_ver = sample["cross-sections-ver"]
            for i_q, cross_section_hor in cross_sections_hor:
                ax_hor.plot(*self.vectors_to_xy(cross_section_hor), label="q = " + str(i_q))
            for i_p, cross_section_ver in cross_sections_ver:
                ax_ver.plot(*self.vectors_to_xy(cross_section_ver), label="p = " + str(i_p))
            ax_hor.legend()
            ax_ver.legend()
    
        axcolor = 'lightgoldenrodyellow'	
        if self.plot_3d and plot_input:
            fig = plt.figure(figsize=plt.figaspect(0.5))
            fig_in = plt.figure(figsize=plt.figaspect(0.5))
            self.ax_2d = fig.add_subplot(1, 2, 1)
            self.ax_3d = fig.add_subplot(1, 2, 2, projection='3d')
            ax_in = fig_in.add_subplot(1, 1, 1)
            ax_slider = fig.add_axes([0.01, 0.05, 1, 0.03], axisbg=axcolor)
            curr_trav = Slider(ax_slider, '', 0.1, 9.0, valinit=0)
        elif self.plot_3d:
            fig = plt.figure(figsize=plt.figaspect(0.5))
            self.ax_2d = fig.add_subplot(1, 2, 1)
            self.ax_3d = fig.add_subplot(1, 2, 2, projection='3d')
            ax_slider = fig.add_axes([0.01, 0.05, 1, 0.03], axisbg=axcolor)
            curr_trav = Slider(ax_slider, '', 0.1, 9.0, valinit=0)
        elif plot_input:
            fig = plt.figure(figsize=plt.figaspect(0.5))
            fig_in = plt.figure(figsize=plt.figaspect(0.5))
            self.ax_2d = fig.add_subplot(1, 1, 1)
            ax_in = fig_in.add_subplot(1, 1, 1)
            ax_slider = fig.add_axes([0.01, 0.05, 1, 0.03], axisbg=axcolor)
            curr_trav = Slider(ax_slider, '', 0.1, 9.0, valinit=0)
        else:
            fig = plt.figure(figsize=plt.figaspect(0.5))
            self.ax_2d = fig.add_subplot(1, 1, 1)
            ax_slider = fig.add_axes([0.01, 0.05, 1, 0.03], axisbg=axcolor)
            curr_trav = Slider(ax_slider, '', 0.1, 9.0, valinit=0)
    
        # set aspect ratio
        self.ax_2d.set_aspect('equal')
        if self.plot_3d:
            self.ax_3d.set_aspect('equal')
        if plot_input:
            ax_in.set_aspect('equal')
        
    
        # show axis labels
        if show_labels:
            self.ax_2d.set_xlabel("Path P")
            self.ax_2d.set_ylabel("Path Q")
    
            if self.plot_3d:
                self.ax_3d.set_xlabel("Path P")
                self.ax_3d.set_ylabel("Path Q")
                self.ax_3d.set_zlabel("Length l")
    
            if plot_input:
                ax_in.set_xlabel("X")
                ax_in.set_ylabel("Y")
    
        # plot input #not in heatmap
        if plot_input:
            paths = sample["input"]
            xa, ya = np.array(self.vectors_to_xy(paths[0]))
            xb, yb = np.array(self.vectors_to_xy(paths[1]))
    
            ax_in.quiver(xa[:-1], ya[:-1], xa[1:] - xa[:-1], ya[1:] - ya[:-1], color="b", scale_units='xy', angles='xy', scale=1)
            ax_in.quiver(xb[:-1], yb[:-1], xb[1:] - xb[:-1], yb[1:] - yb[:-1], color="c", scale_units='xy', angles='xy', scale=1)
    
            ax_in.plot([], [], "b", label="Path P")
            ax_in.plot([], [], "c", label="Path Q")
    
            xlim = [min(xa.min(), xb.min()), max(xa.max(), xb.max())]
            xd = xlim[1] - xlim[0]
            xpad = xd * padding
            ylim = [min(ya.min(), yb.min()), max(ya.max(), yb.max())]
            yd = ylim[1] - ylim[0]
            ypad = yd * padding
            ax_in.axis([xlim[0] - xpad, xlim[1] + xpad, ylim[0] - ypad, ylim[1] + ypad])
    
            ax_in.legend()
            
            #traversals not in heatmap
            if plot_traversals:
                in_traversal = sample["traversal"]["in-traversal"]
                in_traversal_l = sample["traversal"]["in-traversal-l"]
                for ps in in_traversal:
                    x, y = self.vectors_to_xy(ps)
                    ax_in.plot(x, y, "k", linewidth=0.5)
                for ps in in_traversal_l:
                    x, y = self.vectors_to_xy(ps)
                    ax_in.plot(x, y, "r", linewidth=1.5)
    
        # set padding
        # 2d
        l_a, l_b = sample["size"]
        pad_x = padding * l_a
        pad_y = padding * l_b
        self.axis_2d = [-pad_x, l_a + pad_x, -pad_y, l_b + pad_y]
        self.ax_2d.axis(self.axis_2d)
        self.bounds_l = sample["bounds-l"]
        # 3d
        if self.plot_3d:
            dl = self.bounds_l[1] - self.bounds_l[0]
            self.pad_l = padding * dl
            self.ax_3d.axis(self.axis_2d)
            self.ax_3d.set_zlim([self.bounds_l[0] - self.pad_1, self.bounds_l[1] + self.pad_1])
    
        #plot 3d
        if self.plot_3d:
           heatmap = sample["heatmap"]
           x = heatmap[0]
           y = heatmap[1]
           z = heatmap[2]
           surf = self.ax_3d.plot_surface(x, y, z, cmap=cm.coolwarm, rstride=1, cstride=1,
                                     linewidth=0, antialiased=False)
    
        # plot borders
        if plot_borders:
            for border in sample["borders-v"]:
                x, y = self.vectors_to_xy(border[1])
                self.ax_2d.plot(x, y, "", color="0.5")
            for border in sample["borders-h"]:
                x, y = self.vectors_to_xy(border[1])
                self.ax_2d.plot(x, y, "", color="0.5")
    
        # plot cells
        for cell in sample["cells"]:
            data = cell[1]
    
            # plot ellipses
            if plot_ellipsis:
                for ellipsis_to_plot in data["ellipses"]:
                    l = ellipsis_to_plot[0]
                    if len(ellipsis_to_plot[1]) > 0:
                        x, y = self.vectors_to_xy(ellipsis_to_plot[1])
                        if ellipsis_to_plot[0] == 0:
                            self.ax_2d.plot(x, y, "k.")
                        else:
                            self.ax_2d.plot(x, y, "k", linewidth=0.8)
            # plot axis
            if plot_axis:
                for axis in data["axis"]:
                    x, y = self.vectors_to_xy(axis[1])
                    self.ax_2d.plot(x, y, "y:")
            # plot l-lines
            if plot_l_lines:
                for l_line in data["l-lines"]:
                    x, y = self.vectors_to_xy(l_line[1])
                    self.ax_2d.plot(x, y, "c--")
    
        # plot heatmap
        if plot_heatmap:
            heatmap = sample["heatmap"]
            x = heatmap[0]
            y = heatmap[1]
            z = heatmap[2]
            surf = self.ax_2d.pcolor(x, y, z, cmap=cm.coolwarm)
            
    
    
        # plot colorbar
        if (plot_heatmap or self.plot_3d) and show_colorbar:
            fig.colorbar(surf, shrink=0.95, aspect=5)
        ####################################################################################################################
        # plot critical traversals
        if plot_critical_traversals:
            for tra in sample["critical-traversals"]:
                p1 = tra.a
                p2 = tra.b
                self.ax_2d.plot([p1.x], [p1.y], "b.")
                self.ax_2d.plot([p2.x], [p2.y], "b.")
                x, y = self.vectors_to_xy(tra.points)
                self.ax_2d.plot(x, y, "b--", linewidth=1.0)
    
        # plot traversals in 2d
        if plot_traversals:
            for traversal in sample["traversals"]:
                for i in range(len(traversal.points)):
                    if about_equal(traversal.epsilon, traversal.epsilons[i]):
                        point = traversal.points[i]
                        self.ax_2d.plot([point.x], [point.y], "ro")
                x, y = self.vectors_to_xy(traversal.points)
                self.ax_2d.plot(x, y, "r--", label="traversal: l=" + str(traversal.epsilon), linewidth=1.5)
    
        # plot traversal in 3d
        if self.plot_3d and plot_traversals:
            x, y, z = sample["traversal"]["traversal-3d"]
            x_l, y_l, z_l = sample["traversal"]["traversal-3d-l"]
    
            self.ax_3d.plot(x, y, z, "r", label="traversal: l=" + str(sample["traversal"]["epsilon"]), linewidth=0.5)
    
            for i in range(len(x_l)):
                self.ax_3d.plot([x_l[i]] * 2, [y_l[i]] * 2, self.bounds_l, "k", linewidth=0.5)
        #####################################################################################################################
        # show legend in 3d
        if show_legend and self.plot_3d:
            self.ax_3d.legend()
            
            
        curr_trav.on_changed(self.update(self, plot_heatmap, plot_3d, sample))            
        
        plt.show()
	
class DiscreteSlider(Slider):
    """A matplotlib slider widget with discrete steps."""
    def __init__(self, *args, **kwargs):
        """Identical to Slider.__init__, except for the "increment" kwarg.
        "increment" specifies the step size that the slider will be discritized
        to."""
        self.inc = kwargs.pop('increment', 0.5)
        Slider.__init__(self, *args, **kwargs)

    def set_val(self, val):
        discrete_val = int(val)
        # We can't just call Slider.set_val(self, discrete_val), because this 
        # will prevent the slider from updating properly (it will get stuck at
        # the first step and not "slide"). Instead, we'll keep track of the
        # the continuous value as self.val and pass in the discrete value to
        # everything else.
        xy = self.poly.xy
        xy[2] = discrete_val, 1
        xy[3] = discrete_val, 0
        self.poly.xy = xy
        self.valtext.set_text(self.valfmt % discrete_val)
        if self.drawon: 
            self.ax.figure.canvas.draw()
        self.val = val
        if not self.eventson: 
            return
        for cid, func in self.observers.iteritems():
            func(discrete_val)
