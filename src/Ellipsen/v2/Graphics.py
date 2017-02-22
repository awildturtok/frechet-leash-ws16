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
from matplotlib import cm
from numbers import Number
import numpy as np


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


def xy_to_vectors(xs: [float], ys: [float]) -> [Vector]:
    vectors = []

    for i in range(min(len(xs), len(ys))):
        x = xs[i]
        y = ys[i]
        if isinstance(x, Number) and isinstance(y, Number):
            vectors.append(Vector(x, y))
        else:
            print("Error: either isn not a valid float: x:" + str(x) + " y:" + str(y))

    return vectors



def sample_to_matplotlib(sample, plot_borders: bool = True, plot_ellipsis: bool = True, plot_heatmap: bool = True,
                         plot_traversals: bool = True, plot_axis: bool = False, plot_l_lines: bool = False,
                         plot_3d: bool = True, show_legend: bool = False, show_colorbar: bool = True,
                         plot_input: bool = True):

    # plot sample with matplotlib
    """
    Display sampled data with matplotlib.
    :param sample: Hashmap containing the input data:
                    - Input: The pair of curves to produce the sample
                    - in-traversal: Route/Alignment of the lexicographic frechet path
                    - heatmap: The resulting heightmap from the resulting ellipses
                    - cells: The grid cells of the heightmap
    :param plot_borders:
    :param plot_ellipsis:
    :param plot_heatmap:
    :param plot_traversals:
    :param plot_axis:
    :param plot_l_lines:
    :param plot_3d:
    :param show_legend:
    :param show_colorbar:
    :param plot_input:
    """
    if plot_3d and plot_input:
        fig = plt.figure(figsize=plt.figaspect(0.5))
        fig_in = plt.figure(figsize=plt.figaspect(0.5))
        ax_2d = fig.add_subplot(1, 2, 1)
        ax_3d = fig.add_subplot(1, 2, 2, projection='3d')
        ax_in = fig_in.add_subplot(1, 1, 1)
    elif plot_3d:
        fig = plt.figure(figsize=plt.figaspect(0.5))
        ax_2d = fig.add_subplot(1, 2, 1)
        ax_3d = fig.add_subplot(1, 2, 2, projection='3d')
    elif plot_input:
        fig = plt.figure(figsize=plt.figaspect(0.5))
        fig_in = plt.figure(figsize=plt.figaspect(0.5))
        ax_2d = fig.add_subplot(1, 1, 1)
        ax_in = fig_in.add_subplot(1, 1, 1)
    else:
        fig = plt.figure(figsize=plt.figaspect(0.5))
        ax_2d = fig.add_subplot(1, 1, 1)

    # plot input
    if plot_input:
        paths = sample["input"]
        xa, ya = np.array(vectors_to_xy(paths[0]))
        xb, yb = np.array(vectors_to_xy(paths[1]))

        ax_in.quiver(xa[:-1], ya[:-1], xa[1:] - xa[:-1], ya[1:] - ya[:-1], color="b", scale_units='xy', angles='xy', scale=1)
        ax_in.quiver(xb[:-1], yb[:-1], xb[1:] - xb[:-1], yb[1:] - yb[:-1], color="c", scale_units='xy', angles='xy', scale=1)

        ax_in.plot([], [], "b", label="Path A")
        ax_in.plot([], [], "c", label="Path B")
        ax_in.legend()

        if plot_traversals:
            in_traversal = sample["in-traversal"]
            for ps in in_traversal:
                x, y = vectors_to_xy(ps)
                ax_in.plot(x, y, "k", linewidth=0.5)

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
        fig.colorbar(surf, shrink=0.95, aspect=5)

    # set padding
    l_a, l_b = sample["size"]
    pad = 0.04
    pad_x = pad * l_a
    pad_y = pad * l_b
    ax_2d.axis([-pad_x, l_a + pad_x, -pad_y, l_b + pad_y])

    # plot traversals
    if plot_traversals:
        for tra in sample["traversals"]:
            for i in range(len(tra[2])):
                if tra[0] <= tra[1][i] + 1e-13:
                    p = tra[2][i]
                    ax_2d.plot([p.x], [p.y], "r.")
            x, y = vectors_to_xy(tra[2])
            ax_2d.plot(x, y, "r--", label="traversal: l=" + str(tra[0]), linewidth=1.5)
            if plot_3d:
                ax_3d.plot(x, y, tra[1], "r", label="traversal: l=" + str(tra[0]), linewidth=0.5)

    # show legend
    if show_legend:
        #ax_2d.legend()
        ax_3d.legend()

    plt.show()

