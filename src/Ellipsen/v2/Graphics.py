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

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm


def vectors_to_xy(vectors: [Vector]) -> ([float], [float]):
    x = []
    y = []

    for vector in vectors:
        if isinstance(vector, Vector):
            x.append(vector.x)
            y.append(vector.y)
        else:
            print("Error: not a Vector: " + str(vector))

    return x, y


def sample_to_matplotlib(sample, plot_borders: bool = True, plot_ellipsis: bool = True, plot_heatmap: bool = True,
                         plot_traversals: bool = True, plot_axis: bool = False, plot_l_lines: bool = False,
                         plot_3d: bool = True, show_legend: bool = False, show_colorbar: bool = True):

    # plot sample with matplotlib
    if plot_3d:
        fig = plt.figure(figsize=plt.figaspect(0.5))
        ax_2d = fig.add_subplot(1, 2, 1)
        ax_3d = fig.add_subplot(1, 2, 2, projection='3d')
    else:
        fig = plt.figure()
        ax_2d = fig.add_subplot(1, 1, 1)

    # plot 3d
    if plot_3d:
        heatmap = sample["heatmap"]
        x = heatmap[0]
        y = heatmap[1]
        z = heatmap[2]
        surf = ax_3d.plot_surface(x, y, z, cmap=cm.coolwarm, rstride=1, cstride=1,
                                  linewidth=0, antialiased=False)

    # plot borders
    if plot_borders:
        for border in sample["borders-v"]:
            x, y = vectors_to_xy(border[1])
            ax_2d.plot(x, y, "", color="0.5")
        for border in sample["borders-h"]:
            x, y = vectors_to_xy(border[1])
            ax_2d.plot(x, y, "", color="0.5")

    # plot cells
    for cell in sample["cells"]:
        name = cell[0]
        data = cell[1]

        # plot ellipses
        if plot_ellipsis:
            for ellipsis in data["ellipses"]:
                if len(ellipsis[1]) > 0:
                    x, y = vectors_to_xy(ellipsis[1])
                    num = len(x)
                    if num == 1:
                        ax_2d.plot(x, y, "k.", linewidth=0.8)
                    else:
                        ax_2d.plot(x, y, "k", linewidth=0.8)
        # plot axis
        if plot_axis:
            for axis in data["axis"]:
                x, y = vectors_to_xy(axis[1])
                ax_2d.plot(x, y, "y:")
        # plot l-lines
        if plot_l_lines:
            for l_line in data["l-lines"]:
                x, y = vectors_to_xy(l_line[1])
                ax_2d.plot(x, y, "c--")

    # plot heatmap
    if plot_heatmap:
        heatmap = sample["heatmap"]
        x = heatmap[0]
        y = heatmap[1]
        z = heatmap[2]
        surf = ax_2d.pcolor(x, y, z, cmap=cm.coolwarm)

    # plot colorbar
    if (plot_heatmap or plot_3d) and show_colorbar:
        fig.colorbar(surf, shrink=0.5, aspect=5)

    # set padding
    l_a, l_b = sample["size"]
    pad = 0.04
    pad_x = pad * l_a
    pad_y = pad * l_b
    ax_2d.axis([-pad_x, l_a + pad_x, -pad_y, l_b + pad_y])

    # plot traversals
    if plot_traversals:
        for traversal in sample["traversals"]:
            x, y = vectors_to_xy(traversal[2])
            ax_2d.plot(x, y, "r--", label=traversal[0], linewidth=1.5)

    # show legend
    if show_legend:
        ax_2d.legend()

    plt.show()
