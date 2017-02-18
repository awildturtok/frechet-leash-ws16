from src.Ellipsen.v2.Algorithm import CellMatrix
from src.Ellipsen.v2.Geometry import *
from src.Ellipsen.v2.Graphics import sample_to_matplotlib

import sys
import csv
import getopt

options, __ = getopt.getopt(sys.argv[1:], "n:l:s:")

n_heights = 7
specific_l = []
n_samples = 100

for name, value in options:
    if name == "-n":
        n_samples = int(value)
    elif name == "-l":
        specific_l = [float for val in value.split(",")]
    elif name == "-s":
        n_samples = int(value)




reader = csv.reader(iter(sys.stdin.readline, ''))

# firstln = next(reader)

# heights = int(firstln[0])
# samples = int(firstln[1])

paths = []

for row in reader:
    path = []
    for cell in row:
        coords = cell.split(";")
        path.append(Vector(float(coords[0]), float(coords[1])))

    paths.append(path)

    if len(paths) == 2:
        break

patha = paths[0]
pathb = paths[1]

input1 = CellMatrix(patha, pathb)

if len(specific_l) > 0:
    sample1 = input1.sample(specific_l, n_samples)
else:
    sample1 = input1.sample_l(n_heights, n_samples)


#print(sample1)
sample_heatmap1 = input1.sample_heatmap_a(n_samples)
sample1["heatmap"] = sample_heatmap1
#print(sample_heatmap1)
sample_to_matplotlib(sample1, plot_3d=True, show_legend=True)