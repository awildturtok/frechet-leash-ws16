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

'''# new type of events ex. 1
# path_a:
a_xs = [0, 1.5, 1, 2]
a_ys = [0, 0, 1+math.sqrt(0.5), 1+math.sqrt(0.5)]
# path_b:
b_xs = [0, 2, 1]
b_ys = [math.sqrt(0.5), math.sqrt(0.5), 0]'''

'''# new type of event ex. 2
# path_a:
a_xs = [4, 1.7, 2, 0]
a_ys = [1-math.sqrt(0.5), 1-math.sqrt(0.5), 0, 0]
# path_b:
b_xs = [4, 0]
b_ys = [1, 1]'''

'''# old type of event acute:
# path_a:
a_xs = [0, 10]
a_ys = [0, 0]
# path_b:
b_xs = [1, 4, 7, 10]
b_ys = [0, 3, -3, 0]'''

'''# old type of event not acute:
# path_a:
a_xs = [0, 10]
a_ys = [0, 0]
# path_b:
b_xs = [1, 7, 4, 10]
b_ys = [0, 3, -3, 0]'''

'''# critical events test 1
# path_a:
a_xs = [5, 10, 10, 0, 0, 10]
a_ys = [10, 10, 0, 0, 5, 5]
# path_b:
b_xs = [-10, -5, -5, -15, -15, -5]
b_ys = [10, 10, 5, 5, 0, 0]'''

'''# critical events test 2
# path_a:
a_xs = [0, 20]
a_ys = [0, 0]
# path_b:
b_xs = [10, 15, 20]
b_ys = [0, 10, 0]'''

'''# critical events test 3
# path_a:
a_xs = [0, 10, 3, 0, 10]
a_ys = [0, 0, 8, 10, 10]
# path_b:
b_xs = [0, 10]
b_ys = [3, 3]'''

# random test 1
# path_a:
a_xs = np.random.random_sample(8)
a_ys = np.random.random_sample(8)
# path_b:
b_xs = np.random.random_sample(7)
b_ys = np.random.random_sample(7)

'''# random test save 1
# path_a:
a_xs = [0.4772922092138784, 0.01530327285512556, 0.18623288333654608, 0.9571042784382312, 0.7531384890787559, 0.5586263445891971, 0.4087144781918828, 0.742737156990017]
a_ys = [0.9717921811331545, 0.2352575436671138, 0.5101121131543056, 0.04746020008338181, 0.45603230757057844, 0.44729409272770604, 0.4203528406682898, 0.27839889123856565]
# path_b:
b_xs = [0.46355916536642316, 0.7848048532306606, 0.8984278039692681, 0.07598172633778322, 0.2466328327882239, 0.04643325875918691, 0.14385640913682562]
b_ys = [0.2252646424367046, 0.3992611940661893, 0.1551268264371849, 0.7049183480244169, 0.7245003824862535, 0.5009416536371276, 0.2867299714276945]'''

'''# random test save 2
# path_a:
a_xs = [0.47624109293124584, 0.2837557133906694, 0.7135247957391838, 0.17614720063872424, 0.9887899709217016, 0.09772925938927324, 0.8221920560820938, 0.2368133627598178]
a_ys = [0.6917492186536913, 0.6919972280719934, 0.025707124294570338, 0.896297605839744, 0.6505998191096224, 0.5765972821537704, 0.0035238865886881854, 0.6612671752727163]
# path_b:
b_xs = [0.41887808669385074, 0.8726570820763088, 0.6650049068733719, 0.653036889162294, 0.9374457599034023, 0.8068954043422238, 0.08663305551417921]
b_ys = [0.90104618488335, 0.9961534911686879, 0.7135675484774522, 0.6274281145154397, 0.6881091124705155, 0.7559564304821617, 0.26997753895074517]'''


path_a = xy_to_vectors(a_xs, a_ys)
path_b = xy_to_vectors(b_xs, b_ys)

input1 = CellMatrix(path_a, path_b, traverse=0)
print(input1)

sample1 = input1.sample_l(-1, n_p)
#sample1 = input1.sample([0.6322879582403348], n_p)

sample_to_matplotlib(sample1, plot_3d=False, plot_traversals=False, plot_critical_traversals=True, plot_l_lines=False, plot_borders=False)
