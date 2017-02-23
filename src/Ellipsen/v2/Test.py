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

from Graphics import *

import numpy as np


n_l = 9  # number of height lines
n_p = 50  # points per ellipsis


# path_a:
a_xs = [10, 0, 10, 0]
a_ys = [-8, -8, 2, 2]

# path_b:
b_xs = [10, 0]
b_ys = [0, 0]

'''# path_a:
a_xs = np.random.random_sample(6)
a_ys = np.random.random_sample(6)

# path_b:
b_xs = np.random.random_sample(6)
b_ys = np.random.random_sample(6)'''


path_a = xy_to_vectors(a_xs, a_ys)
path_b = xy_to_vectors(b_xs, b_ys)

input1 = CellMatrix(path_a, path_b)
print(input1)

sample1 = input1.sample_l(n_l, n_p)

sample_to_matplotlib(sample1)
