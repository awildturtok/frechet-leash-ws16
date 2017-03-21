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


n_l = 7  # number of height lines
n_p = 50  # points per ellipsis

'''# path_a:
a_xs = [0, 1]
a_ys = [0, 1]

# path_b:
b_xs = [1, 0]
b_ys = [0, 0]'''

'''# path_a:
a_xs = [-8, 12, 12, 0, -10]
a_ys = [-2, -2, 2, 10, 5]

# path_b:
b_xs = [-10, 10, -10]
b_ys = [0, 0, 4]'''

# new type of events ex. 1
# path_a:
a_xs = [0, 1.5, 1, 2]
a_ys = [0, 0, 1+math.sqrt(0.5), 1+math.sqrt(0.5)]

# path_b:
b_xs = [0, 2, 1]
b_ys = [math.sqrt(0.5), math.sqrt(0.5), 0]

'''# new type of event ex. 2
# path_a:
a_xs = [4, 1.7, 2, 0]
a_ys = [1-math.sqrt(0.5), 1-math.sqrt(0.5), 0, 0]

# path_b:
b_xs = [4, 0]
b_ys = [1, 1]'''

'''# path_a:
a_xs = np.random.random_sample(10)
a_ys = np.random.random_sample(10)

# path_b:
b_xs = np.random.random_sample(10)
b_ys = np.random.random_sample(10)'''


path_a = xy_to_vectors(a_xs, a_ys)
path_b = xy_to_vectors(b_xs, b_ys)

input1 = CellMatrix(path_a, path_b, traverse=0)
print(input1)

#sample1 = input1.sample_l(n_l, n_p)
sample1 = input1.sample([1], n_p)

sample_to_matplotlib(sample1, plot_3d=True, plot_traversals=False, plot_critical_traversals=False, plot_l_lines=False)
