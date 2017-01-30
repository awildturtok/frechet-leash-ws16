from math import floor, ceil

from src.Ellipsen.v2.Algorithm import TwoLineSegments
from src.Types import *
from src.Ellipsen.v2.Geometry import Vector, LineSegment
from src.Ellipsen.v2.Algorithm import intersection

P = list(map(Vector.from_tuple, [(0,0), (1,2), (2, 3), (0,7)]))
Q = list(map(Vector.from_tuple, [(0,0), (1,1), (1, 2), (1,1), (0,7)]))


class Curve:

    def __init__(self, points: List[Vector]):
        self.points = points

    def at(self, x: float) -> Vector:
        base = floor(x)
        return self.points[base] + (self.points[base + 1] - self.points[base]) * (x - base)




def descent(a: Curve, b: Curve, ax: float, bx: float):
    left = LineSegment(a.at(ax), a.at(ceil(ax)))
    right = LineSegment(b.at(bx), b.at(ceil(ax)))

    pair = TwoLineSegments(left, right)
    cell = pair.cell()

    cross = cell.l_hor.d.cross_product(Vector(ax, bx))

    if(cross == 0): # is on line
        print()
    else: #over or under line
        print()


