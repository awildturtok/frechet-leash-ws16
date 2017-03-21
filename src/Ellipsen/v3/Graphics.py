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


def vectors_to_xy(vectors: [Vector]) -> ([float], [float]):  # converts array of vectors to x- & y-coordinate arrays
    x = []
    y = []

    for vector in vectors:
        if isinstance(vector, Vector):
            x.append(vector.x)
            y.append(vector.y)
        else:
            print("Error: not a Vector: " + str(vector))

    return x, y


def xy_to_vectors(xs: [float], ys: [float]) -> [Vector]:  # converts x- & y-coordinate arrays to arrays of vectors
    vectors = []

    for i in range(min(len(xs), len(ys))):
        x = xs[i]
        y = ys[i]
        if isinstance(x, Number) and isinstance(y, Number):
            vectors.append(Vector(x, y))
        else:
            print("Error: either is not a valid float: x:" + str(x) + " y:" + str(y))

    return vectors


def sample_to_matplotlib(sample, plot_borders: bool = True, plot_ellipsis: bool = True, plot_heatmap: bool = True,
                         plot_traversals: bool = True, plot_axis: bool = False, plot_l_lines: bool = False,
                         plot_3d: bool = True, show_legend: bool = True, show_colorbar: bool = True,
                         plot_input: bool = True, show_labels: bool = True, plot_critical_traversals: bool = False):
    # plot sample with matplotlib

    padding = 0.04

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

    # set aspect ratio
    ax_2d.set_aspect('equal')
    #if plot_3d:  # Todo
        #ax_3d.set_aspect('equal')
    if plot_input:
        ax_in.set_aspect('equal')

    # show axis labels
    if show_labels:
        ax_2d.set_xlabel("Path A")
        ax_2d.set_ylabel("Path B")

        if plot_3d:
            ax_3d.set_xlabel("Path A")
            ax_3d.set_ylabel("Path B")
            ax_3d.set_zlabel("Length l")

        if plot_input:
            ax_in.set_xlabel("X")
            ax_in.set_ylabel("Y")

    # plot input
    if plot_input:
        paths = sample["input"]
        xa, ya = np.array(vectors_to_xy(paths[0]))
        xb, yb = np.array(vectors_to_xy(paths[1]))

        ax_in.quiver(xa[:-1], ya[:-1], xa[1:] - xa[:-1], ya[1:] - ya[:-1], color="b", scale_units='xy', angles='xy', scale=1)
        ax_in.quiver(xb[:-1], yb[:-1], xb[1:] - xb[:-1], yb[1:] - yb[:-1], color="c", scale_units='xy', angles='xy', scale=1)

        ax_in.plot([], [], "b", label="Path A")
        ax_in.plot([], [], "c", label="Path B")

        xlim = [min(xa.min(), xb.min()), max(xa.max(), xb.max())]
        xd = xlim[1] - xlim[0]
        xpad = xd * padding
        ylim = [min(ya.min(), yb.min()), max(ya.max(), yb.max())]
        yd = ylim[1] - ylim[0]
        ypad = yd * padding
        ax_in.axis([xlim[0] - xpad, xlim[1] + xpad, ylim[0] - ypad, ylim[1] + ypad])

        ax_in.legend()

        if plot_traversals:
            in_traversal = sample["traversal"]["in-traversal"]
            in_traversal_l = sample["traversal"]["in-traversal-l"]
            for ps in in_traversal:
                x, y = vectors_to_xy(ps)
                ax_in.plot(x, y, "k", linewidth=0.5)
            for ps in in_traversal_l:
                x, y = vectors_to_xy(ps)
                ax_in.plot(x, y, "r", linewidth=1.5)

    # set padding
    # 2d
    l_a, l_b = sample["size"]
    pad_x = padding * l_a
    pad_y = padding * l_b
    axis_2d = [-pad_x, l_a + pad_x, -pad_y, l_b + pad_y]
    ax_2d.axis(axis_2d)
    bounds_l = sample["bounds-l"]
    # 3d
    if plot_3d:
        dl = bounds_l[1] - bounds_l[0]
        pad_l = padding * dl
        ax_3d.axis(axis_2d)
        ax_3d.set_zlim([bounds_l[0] - pad_l, bounds_l[1] + pad_l])

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
                        ax_2d.plot(x, y, "k.")
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

    # plot critical traversals
    if plot_critical_traversals:
        for tra in sample["critical-traversals"]:
            p1 = tra[3][0]
            p2 = tra[3][-1]
            ax_2d.plot([p1.x], [p1.y], "b.")
            ax_2d.plot([p2.x], [p2.y], "b.")
            x, y = vectors_to_xy(tra[3])
            ax_2d.plot(x, y, "b--", linewidth=1.0)

    # plot traversals
    if plot_traversals:

        for tra in sample["traversals"]:
            for i in range(len(tra[2])):
                if tra[0] <= tra[2][i] + 1e-13:
                    p = tra[3][i]
                    ax_2d.plot([p.x], [p.y], "ro")
            x, y = vectors_to_xy(tra[3])
            ax_2d.plot(x, y, "r--", label="traversal: l=" + str(tra[0]), linewidth=1.5)

    if plot_3d and plot_traversals:
        x, y, z = sample["traversal"]["traversal-3d"]
        x_l, y_l, z_l = sample["traversal"]["traversal-3d-l"]

        ax_3d.plot(x, y, z, "r", label="traversal: l=" + str(tra[0]), linewidth=0.5)

        for i in range(len(x_l)):
            ax_3d.plot([x_l[i]] * 2, [y_l[i]] * 2, bounds_l, "k", linewidth=0.5)

    # show legend
    if show_legend and plot_3d:
        ax_3d.legend()

    plt.show()
