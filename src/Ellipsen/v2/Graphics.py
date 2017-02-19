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
                         plot_3d: bool = True, show_legend: bool = False, show_colorbar: bool = True,
                         plot_input: bool = True):

    # plot sample with matplotlib
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
        xa, ya = vectors_to_xy(paths[0])
        xb, yb = vectors_to_xy(paths[1])

        ax_in.plot(xa, ya, "b", label="Path A", linewidth=2.5)
        ax_in.plot(xb, yb, "c", label="Path B", linewidth=2.5)
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
        fig.colorbar(surf, shrink=0.95, aspect=5)

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
            if plot_3d:
                ax_3d.plot(x, y, traversal[1], "r", label=traversal[0], linewidth=0.5)

    # show legend
    if show_legend:
        ax_2d.legend()
        ax_3d.legend()

    plt.show()


'''ap1 = Vector(0, 0)
ap2 = Vector(1.2, -0.2)
ap3 = Vector(1, 1)
ap4 = Vector(2.2, 0.8)
ap5 = Vector(2, 2)
ap6 = Vector(3.2, 1.8)
ap7 = Vector(3, 3)

bp1 = Vector(0, 0)
bp2 = Vector(-0.2, 1.2)
bp3 = Vector(1, 1)
bp4 = Vector(0.8, 2.2)
bp5 = Vector(2, 2)
bp6 = Vector(1.8, 3.2)
bp7 = Vector(3, 3)'''


# Aktueller Test:

'''ap1 = Vector(0, 0)
ap2 = Vector(20, -3)
ap3 = Vector(-10, 4)
ap4 = Vector(3, 4)
ap5 = Vector(-10, -2)


bp1 = Vector(0, 5)
bp2 = Vector(10, 3)
bp3 = Vector(4, -12)
bp4 = Vector(16, 5)
bp5 = Vector(-12, 0)

patha = [ap1, ap2, ap3, ap4, ap5]  # , ap6, ap7]
pathb = [bp1, bp2, bp3, bp4, bp5]  # , bp6, bp7]

input1 = CellMatrix(patha, pathb)
print(input1)
sample1 = input1.sample_l(9, 100)
#print(sample1)
sample_heatmap1 = input1.sample_heatmap_a(100)
sample1["heatmap"] = sample_heatmap1
#print(sample_heatmap1)
traversal = sample1["traversals"][0]
sample_traversal = input1.sample_traversal(traversal, 100)
sample1["in-traversal"] = sample_traversal
print(sample_traversal)
sample_to_matplotlib(sample1, plot_3d=True, show_legend=False)'''
