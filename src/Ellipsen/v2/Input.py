from src.Ellipsen.v2.Algorithm import CellMatrix
from src.Ellipsen.v2.Geometry import *
from src.Ellipsen.v2.Graphics import sample_to_matplotlib

import sys
import csv

reader = csv.reader(iter(sys.stdin.readline, ''))

firstln = next(reader)

heights = firstln[0]
samples = firstln[1]

paths = []

for row in reader:
    path = []
    for cell in row:
        coords = cell.split(";")
        path.append(Vector(float(coords[0]), float(coords[1])))

    paths.append(path)

    if len(paths) == 2:
        break

print(paths)

patha = paths[0]
pathb = paths[1]

input1 = CellMatrix(patha, pathb)
print(input1)
sample1 = input1.sample(heights, samples)
#print(sample1)
sample_heatmap1 = input1.sample_heatmap_a(samples)
sample1["heatmap"] = sample_heatmap1
#print(sample_heatmap1)
sample_to_matplotlib(sample1, plot_3d=True, show_legend=True)