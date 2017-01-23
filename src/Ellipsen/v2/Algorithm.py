########################################################################################################################
#                                                                                                                      #
#                                    Software Projekt: Frechet Distanz                                                 #
#                                    Teilgebiet: Ellipsen-Alg. einer Zelle                                             #
#                                    Erstellt: WS 16/17 FU Berlin                                                      #
#                                                                                                                      #
#                                    Team: Josephine Mertens, Jana Kirschner,                                          #
#                                    Alexander Korzech, Fabian Kovacs, Alexander                                       #
#                                    Timme, Kilian Kraatz, Anton Begehr                                                #
#                                                                                                                      #
########################################################################################################################

# -*- coding: utf-8 -*-

from Geometry import *
import plotly.plotly as py
import plotly.graph_objs as go


class OneLineSegment(LineSegment):  # one of the two input line segments
    def __init__(self, ls: LineSegment, s: Vector):
        super().__init__(ls.p1, ls.p2)

        self.intersects = self.contains_point(s)  # intersects S?
        self.dir = (self.r_point(s) <= 1)  # points in direction  of S?
        self.rs = self.r_point(s)  # S in parametric terms of line segment

    def __str__(self):
        return str(self.p1) + "->" + str(self.p2) + ": m=" + str(self.m) + " n=" + str(self.n) + \
               "\n   Intersects:" + str(self.intersects) + " Direction:" + str(self.dir) + \
               " r(S): " + str(self.rs)  # + " Offset:" + str(self.offset)


class Cell:
    def __init__(self, c: Vector, d: Vector, offset: Vector, bounds_xy: (float, float), bounds_l: (float, float)):
        self.norm_ellipsis = Ellipse(offset, c, d)  # Ellipsis for l = 1
        self.bounds_xy = bounds_xy  # cell bounds
        self.bounds_l = bounds_l  # line length bounds

    def __str__(self):
        return "Cell:\n   Norm-" + str(self.norm_ellipsis) + '\n' + \
               "   Bounds l: " + str(self.bounds_l) + '\n' + \
               "   Bounds XY: " + str(self.bounds_xy)

    def sample(self, n1: int, n2: int, rel_bounds: ((float, float), (float, float)) = ((0, 1), (0, 1)))\
            -> [(float, [Vector])]:
        # n1: # of ellipsis, n2: # of points per ellipsis, bounds: relative xy bounds

        bounds = ((rel_bounds[0][0] * self.bounds_xy[0], rel_bounds[0][1] * self.bounds_xy[0]),
                  (rel_bounds[1][0] * self.bounds_xy[1], rel_bounds[1][1] * self.bounds_xy[1]))
        ellipses_sample = []  # holds ellipses in form: (l, [Points])

        for i1 in range(n1):
            l =  1 # self.bounds_l[0] + (float(i1) / (n1 - 1)) * (self.bounds_l[1] - self.bounds_l[0])
            ellipsis = self.norm_ellipsis * l
            ellipsis_sample = []
            print(ellipsis)
            for i2 in range(n2):
                t = (i2 / (n2 - 1)) * (2 * math.pi)
                p = ellipsis.p(t)
                print(str(t) + "->" + str(p))
                if p.in_bounds(bounds):
                    ellipsis_sample.append(p)
            ellipses_sample.append((l, ellipsis_sample))

        return ellipses_sample


def intersection(a: LineSegment, b: LineSegment):  # calculates the intersection point of two line segments
    if math.isinf(a.m) and not math.isinf(b.m):
        return Vector(a.n, b.fx(a.n))
    elif math.isinf(b.m) and not math.isinf(a.m):
        return Vector(b.n, a.fx(b.n))
    elif (a.m - b.m) != 0:
        x = (b.n - a.n) / (a.m - b.m)
    else:
        x = float("nan")
    return Vector(x, a.fx(x))


class TwoLineSegments:  # Input: two line segments
    def __init__(self, a: LineSegment, b: LineSegment):
        # calculate shortest and longest possible line length - FALSCH: ex. 0,0->0,2 und 1,1->2,2
        self.bounds_l = (min(a.d_ls_point(b.p1), a.d_ls_point(b.p2), b.d_ls_point(a.p1), b.d_ls_point(a.p2)),
                         max(a.p1.d(b.p1), a.p1.d(b.p2), a.p2.d(b.p1), a.p2.d(b.p2)))

        self.parallel = a.m == b.m
        if not self.parallel:  # case 1: lines are not parallel
            self.s = intersection(a, b)  # intersection point S

            self.a = OneLineSegment(a, self.s)  # line segment a
            self.b = OneLineSegment(b, self.s)  # line sqgment b

            self.intersect = self.a.intersects and self.b.intersects  # does the intersection point lie on A and B
            if self.intersect:  # if line segments intersect, set min length to 0
                self.bounds_l = (0.0, self.bounds_l[1])

        else:  # case 2: lines are parallel
            """self.dirB = (a.rP(a.p1+b.d) >= 0)  # do A and B point in the same direction
            self.stWinkel = math.atan(a.m) # Steigungswinkel der Geraden
            # Distanz der 2 Geraden
            if self.a.m == 0 or math.isinf(self.a.m): self.dist = abs(self.a.n - self.b.n)
            else: self.dist = math.sin(self.stWinkel) * abs(a.fy(b.n))
            self.ob = 0 # Offset Strecke B im Vergleich zu A"""
            # Hier weiter Spezialfall: Parallel implementieren

    def __str__(self):
        return "two Line Segments:\nLine Segment A: " + str(self.a) + "\nLine Segment B: " + str(
            self.b) + "\nIntersect (" + str(self.intersect) + "): " + str(self.s)

    def cell(self) -> Cell:
        angle_a = self.a.d.angle()

        # Vectors a and b rotated so a lies on x-axis
        rot_a = self.a.d.rotate(-angle_a)
        rot_b = self.b.d.rotate(-angle_a)

        # Vectors a and b normalised
        norm_a = rot_a.norm()
        norm_b = rot_b.norm()

        # Calculate ellipsis vectors c and d for l = 1
        c = (norm_b + norm_a) * (0.5 / (norm_b - norm_a).l)
        d = (norm_b - norm_a) * (0.5 / (norm_b + norm_a).l)

        # Calculate ellipsis offset:
        offset_a = rot_a * (-abs(self.a.rs))
        offset_b = rot_b * (-abs(self.b.rs))
        offset = offset_a + offset_b

        # set cell bounds: length of a and b
        bounds_xy = (self.a.d.l, self.b.d.l)

        return Cell(c, d, offset, bounds_xy, self.bounds_l)


def sample_to_plotly(sample):  # send sample to plotly
    data = []
    for s in sample:
        l = s[0]
        x = []
        y = []
        for p in s[1]:
            x.append(p.x)
            y.append(p.y)
        trace = go.Scatter(x=x, y=y, name=str(l))
        data.append(trace)
    fig = dict(data=data)
    py.plot(fig, filename='Ellipsen-Test')


st1 = LineSegment(Vector(0, 0), Vector(1, 0))
st2 = LineSegment(Vector(0, 0), Vector(1/math.sqrt(2), 1/math.sqrt(2)))
eingabe1 = TwoLineSegments(st1, st2)
print(eingabe1)
cell1 = eingabe1.cell()
print(cell1)
sample1 = cell1.sample(8, 100, ((-2, 3), (-2, 3)))
#print(sample1)
sample_to_plotly(sample1)
