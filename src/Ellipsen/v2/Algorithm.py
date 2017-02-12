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

from Geometry import *

import plotly.plotly as py
import plotly.graph_objs as go

import matplotlib.pyplot as plt


class OneLineSegment(LineSegment):  # one of the two input line segments
    def __init__(self, ls: LineSegment, s: Vector):
        super().__init__(ls.p1, ls.p2)

        self.intersects = self.contains_point(s)  # intersects S?
        self.dir = (self.r_point(s) <= 1)  # points in direction  of S?
        self.rs = self.r_point(s)  # S in parametric terms of line segment

    def __str__(self):
        return str(self.p1) + "->" + str(self.p2) + ": l=" + str(self.l) + " m=" + str(self.m) + " n=" + str(self.n) + \
               "\n      Intersects:" + str(self.intersects) + " Direction:" + str(self.dir) + " r(S): " + str(self.rs)


class Cell:
    def __init__(self, a: OneLineSegment, b: OneLineSegment, c: Vector, d: Vector, m: Vector, bounds_xy: (float, float),
                 bounds_l: (float, float), offset: Vector = Vector(0, 0), do_traverse: bool = False,
                 start: Vector = Vector(0, 0)):
        self.a = a
        self.b = b
        self.norm_ellipsis = Ellipse(m, c, d)  # Ellipsis for l = 1
        self.bounds_xy = bounds_xy  # cell bounds
        self.bounds_l = bounds_l  # line length bounds
        self.offset = offset
        self.end = Vector(bounds_xy[0], bounds_xy[1])
        self.end_l = self.lp(self.end)

        # steepest decent lines l=l_ver and l'=l_hor
        t_l_ver = math.atan(d.x / c.x)
        t_l_hor = math.atan(d.y / c.y)
        self.l_ver = LineSegment(m, self.norm_ellipsis.p(t_l_ver))
        self.l_hor = LineSegment(m, self.norm_ellipsis.p(t_l_hor))
        self.l_ver_cut = Vector(bounds_xy[0], self.l_ver.fx(bounds_xy[0]))
        self.l_hor_cut = Vector(self.l_hor.fy(bounds_xy[1]), bounds_xy[1])

        self.do_traverse = do_traverse
        if do_traverse:
            self.traversals = self.traverse((self.lp(start-offset), [start]))

    def __str__(self):
        return "    Norm-" + str(self.norm_ellipsis) + '\n' + \
               "    Offset: " + str(self.offset) + '\n' + \
               "    End: " + str(self.end) + '\n' + \
               "    Steepest Decent Lines:\n" + \
               "      l: " + str(self.l_ver.d) + '\n' + \
               "      l': " + str(self.l_hor.d) + '\n' + \
               "    Bounds l: " + str(self.bounds_l) + '\n' + \
               "    Bounds XY: " + str(self.bounds_xy)

    def sample_l(self, nl: int, np: int, rel_bounds: ((float, float), (float, float)) = ((0, 1), (0, 1))) \
            -> [(float, [Vector])]:
        ls = []
        for i in range(nl):
            l = self.bounds_l[0] + (float(i) / (nl - 1)) * (self.bounds_l[1] - self.bounds_l[0])
            ls.append(l)

        return self.sample(ls, np, rel_bounds)

    def sample(self, ls: [float], n: int, rel_bounds: ((float, float), (float, float)) = ((0, 1), (0, 1))) \
            -> [(float, [Vector])]:
        # n: # of points per ellipsis, bounds: relative xy bounds

        bounds = ((rel_bounds[0][0] * self.bounds_xy[0], rel_bounds[0][1] * self.bounds_xy[0]),
                  (rel_bounds[1][0] * self.bounds_xy[1], rel_bounds[1][1] * self.bounds_xy[1]))

        # Sample Ellipses Ellipses old
        '''ellipses_sample = []  # holds ellipses in form: (l, [Points])
        for l in ls:
            if l < self.bounds_l[0] or l > self.bounds_l[1]:
                continue
            ellipsis = self.norm_ellipsis * l
            ellipsis_sample = []
            tps = {}
            for i2 in range(n):
                t = (i2 / (n - 1)) * (2 * math.pi)
                p = ellipsis.p(t)
                if p.in_bounds(bounds) and p not in tps.values():
                    tps[t] = p
            for t in ellipsis.cuts_bounds_t(bounds):
                tps[t] = ellipsis.p(t)
            for t in sorted(tps.keys()):
                ellipsis_sample.append(tps[t] + self.offset)
            ellipses_sample.append((0, ellipsis_sample))  # (l, ellipsis_sample))'''

        # Sample Ellipses
        ellipses_sample = []  # holds sample ellipses & lines in form: (name, [Vector])

        for l in ls:

            if l < self.bounds_l[0] or l > self.bounds_l[1]:
                continue

            if l == 0:
                p = self.norm_ellipsis.m
                if p.in_bounds(bounds):
                    ellipses_sample.append((0, [p + self.offset]))
                continue

            ellipsis = self.norm_ellipsis * l
            ts = ellipsis.cuts_bounds_t(bounds)
            l_ts = len(ts)

            if l_ts == 0:
                ts.append(0)

            for i in range(len(ts)):

                t1 = ts[i-1]
                t2 = ts[i]
                if t1 >= t2:
                    t2 += 2*math.pi

                d_t = t2 - t1
                mid_t = t1 + 0.5*d_t

                if ellipsis.p(mid_t).in_bounds(bounds):
                    ellipsis_sample = [ellipsis.p(t1) + self.offset]

                    np1 = math.ceil((t1 / (2 * math.pi)) * n)
                    np2 = math.ceil((t2 / (2 * math.pi)) * n)

                    for ip in range(np1, np2):
                        t = (ip / n) * 2*math.pi
                        ellipsis_sample.append(ellipsis.p(t) + self.offset)

                    ellipsis_sample.append(ellipsis.p(t2) + self.offset)

                    ellipses_sample.append((0, ellipsis_sample))

        # Sample Ellipsis center-point
        # ellipses_sample.append(("S", [self.norm_ellipsis.m]))

        # Sample Ellipsis axis
        '''ellipses_sample.append(("c", LineSegment(self.norm_ellipsis.m,
                                                 self.norm_ellipsis.m + self.norm_ellipsis.a).cuts_bounds(bounds)))
        ellipses_sample.append(("d", LineSegment(self.norm_ellipsis.m,
                                                 self.norm_ellipsis.m + self.norm_ellipsis.b).cuts_bounds(bounds)))'''

        # Sample steepest decent lines
        '''ellipses_sample.append(("l", [p + self.offset for p in self.l_ver.cuts_bounds(bounds)]))
        ellipses_sample.append(("l'", [p + self.offset for p in self.l_hor.cuts_bounds(bounds)]))'''

        return ellipses_sample

    def lp(self, p: Vector) -> float:  # l for given point
        return self.a.frl(p.x).d(self.b.frl(p.y))

    def traverse_right(self, traversal: (float, [Vector])) -> (float, [Vector]):
        l1x = self.l_ver_cut

        start = traversal[1][-1] - self.offset

        max_l = -1
        end = -self.offset

        if not math.isclose(l1x.y, self.bounds_xy[1], rel_tol=1e-13):
            if l1x.y <= start.y:
                end = Vector(self.bounds_xy[0], start.y)
                max_l = max(traversal[0], self.lp(end))
            elif start.y < l1x.y < self.bounds_xy[1]:
                end = l1x
                max_l = max(traversal[0], self.lp(end))

        return max_l, traversal[1] + [end + self.offset]

    def traverse_top(self, traversal: (float, [Vector])) -> (float, [Vector]):
        l2x = self.l_hor_cut

        start = traversal[1][-1] - self.offset

        max_l = -1
        end = -self.offset

        if not math.isclose(l2x.x, self.bounds_xy[0], rel_tol=1e-13):
            if l2x.x <= start.x:
                end = Vector(start.x, self.bounds_xy[1])
                max_l = max(traversal[0], self.lp(end))
            elif start.x < l2x.x < self.bounds_xy[0]:
                end = l2x
                max_l = max(traversal[0], self.lp(end))

        return max_l, traversal[1] + [end + self.offset]

    def traverse_top_right(self, traversal: (float, [Vector])) -> (float, [Vector]):
        l1x = self.l_ver_cut
        l2x = self.l_hor_cut

        max_l = traversal[0]
        start = traversal[1][-1] - self.offset

        if (l1x.y > self.bounds_xy[1] and l2x.x > self.bounds_xy[0]) or \
                math.isclose(l1x.y, self.bounds_xy[1], rel_tol=1e-13) or \
                math.isclose(l2x.x, self.bounds_xy[0], rel_tol=1e-13):
            end = Vector(self.bounds_xy[0], self.bounds_xy[1])
            max_l = max(max_l, self.lp(end))
        else:
            max_l = -1
            end = -self.offset

        return max_l, traversal[1] + [end + self.offset]

    def traverse(self, traversal: (float, [Vector])) -> [(float, [Vector])]:
        traversal_right = self.traverse_right(traversal)
        traversal_top = self.traverse_top(traversal)
        traversal_top_right = self.traverse_top_right(traversal)

        traversals = [traversal_right, traversal_top_right, traversal_top]

        return traversals


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
        # calculate shortest and longest possible line length
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
        return "    Line Segment A: " + str(self.a) + "\n    Line Segment B: " + str(
            self.b) + "\n     Intersect (" + str(self.intersect) + "): " + str(self.s)

    def cell(self, offset: Vector = Vector(0, 0), do_traverse: bool = False, start:Vector = Vector(0, 0)) -> Cell:
        # Vectors a and b normalised
        norm_a = self.a.d.norm()
        norm_b = self.b.d.norm()

        # x- and y-coordinates for ellipsis vectors c and d for l = 1
        c_xy = 1 / (norm_a - norm_b).l
        d_xy = 1 / (norm_a + norm_b).l

        # Ellipsis vectors c and d for l = 1
        c = Vector(c_xy, c_xy)
        d = Vector(-d_xy, d_xy)

        # Calculate ellipsis offset:
        offset_a = Vector(self.a.l, 0) * self.a.rs
        offset_b = Vector(0, self.b.l) * self.b.rs
        m = offset_a + offset_b

        # set cell bounds: length of a and b
        bounds_xy = (self.a.d.l, self.b.d.l)

        return Cell(self.a, self.b, c, d, m, bounds_xy, self.bounds_l, offset=offset, do_traverse=do_traverse,
                    start=start)


class Input:  # Input: two paths
    def __init__(self, points_a: [Vector], points_b: [Vector]):
        self.path_a = [LineSegment(points_a[i], points_a[i+1]) for i in range(len(points_a)-1)]
        self.path_b = [LineSegment(points_b[i], points_b[i+1]) for i in range(len(points_b)-1)]
        self.count_a = len(self.path_a)
        self.count_b = len(self.path_b)
        self.length_a = 0
        self.length_b = 0
        self.bounds_l = (math.inf, 0)

        self.lengths_a = [0] * self.count_a
        self.lengths_b = [0] * self.count_b
        self.offsets_a = [0] * (self.count_a + 1)
        self.offsets_b = [0] * (self.count_b + 1)

        for i in range(self.count_a):
            a = self.path_a[i]
            self.length_a += a.l
            self.offsets_a[i + 1] = self.offsets_a[i] + a.l
            self.lengths_a[i] = a.l

        for i in range(self.count_b):
            b = self.path_b[i]
            self.length_b += b.l
            self.offsets_b[i + 1] = self.offsets_b[i] + b.l
            self.lengths_b[i] = b.l

        self.twoLSs = []
        self.cells = []

        for i_a in range(self.count_a):
            self.twoLSs.append([])
            self.cells.append([])
            for i_b in range(self.count_b):
                two_line_segments = TwoLineSegments(self.path_a[i_a], self.path_b[i_b])
                self.bounds_l = (min(self.bounds_l[0], two_line_segments.bounds_l[0]),
                                 max(self.bounds_l[1], two_line_segments.bounds_l[1]))
                self.twoLSs[i_a].append(two_line_segments)
                self.cells[i_a].append(two_line_segments.cell(offset=Vector(self.offsets_a[i_a], self.offsets_b[i_b])))

        # traverse
        self.traversals = self.traverse(0, 0, (0, [Vector(0, 0)]))
        print(self.traversals)

    def traverse(self, i_a: int, i_b: int, traversal: (float, [Vector])) -> [(float, [Vector])]:
        if i_a >= self.count_a and i_b >= self.count_b:
            return [traversal]

        traversals = []

        if i_a >= self.count_a or i_b >= self.count_b:
            if i_a < self.count_a:
                cell = self.cells[i_a][self.count_b - 1]
                traversals += self.traverse(i_a + 1, i_b, (max(traversal[0], cell.end_l),
                                                           traversal[1] + [cell.end + cell.offset]))
            if i_b < self.count_b:
                cell = self.cells[self.count_a - 1][i_b]
                traversals += self.traverse(i_a, i_b + 1, (max(traversal[0], cell.end_l),
                                                           traversal[1] + [cell.end + cell.offset]))
        else:
            next_traversal = self.cells[i_a][i_b].traverse(traversal)

            if next_traversal[0][0] != -1:
                traversals += self.traverse(i_a + 1, i_b, next_traversal[0])
            if next_traversal[1][0] != -1:
                traversals += self.traverse(i_a + 1, i_b + 1, next_traversal[1])
            if next_traversal[2][0] != -1:
                traversals += self.traverse(i_a, i_b + 1, next_traversal[2])

        return traversals

    def __str__(self):
        desc = "Input (" + str(self.count_a) + "x" + str(self.count_a) + ") (" + \
                 str(self.length_a) + "x" + str(self.length_b) + "):\n"
        desc += " Path A:\n"
        for i in range(self.count_a):
            desc += "  " + str(i) + ": " + str(self.path_a[i]) + '\n'
        desc += " Path B:\n"
        for i in range(self.count_b):
            desc += "  " + str(i) + ": " + str(self.path_b[i]) + '\n'
        desc += "==>\n"
        desc += " Bounds_l: " + str(self.bounds_l) + '\n'
        desc += " Two Line Segments & Cells:\n"
        for i_a in range(self.count_a):
            for i_b in range(self.count_b):
                desc += " Cell " + str(i_a) + "x" + str(i_b) + '\n'
                desc += str(self.twoLSs[i_a][i_b]) + '\n'
                desc += str(self.cells[i_a][i_b]) + '\n'

        return desc

    def sample(self, nl: int, np: int) -> [(str, [Vector])]:
        ls = []
        samples = []

        for i in range(nl + 1):
            l = self.bounds_l[0] + (float(i) / (nl - 1)) * (self.bounds_l[1] - self.bounds_l[0])
            ls.append(l)

        # sample cells
        for columns in self.cells:
            for cell in columns:
                samples += cell.sample(ls, np)

        # sample cell borders
        '''for i in range(1, self.count_a):  # vertical
            samples.append(("border", [Vector(self.offsets_a[i], 0), Vector(self.offsets_a[i], self.length_b)]))
        for i in range(1, self.count_b):  # horizontal
            samples.append(("border", [Vector(0, self.offsets_b[i]), Vector(self.length_a, self.offsets_b[i])]))'''

        # sample traversals
        for traversal in self.traversals:
            if traversal[0] != -1:
                samples.append(("t: " + str(traversal[0]), traversal[1]))

        return samples


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


def sample_to_matplotlib(sample):  # plot sample with matplotlib
    for s in sample:
        l = s[0]
        x = []
        y = []
        for p in s[1]:
            if isinstance(p, Vector):
                x.append(p.x)
                y.append(p.y)
            else:
                print("Wrong Type: " + str(l) + ": " + str(p))
        len_x = len(x)
        if len_x > 1 and l != 0:
            plt.plot(x, y, label=l)
        elif len_x > 1:
            plt.plot(x, y)
        elif len_x == 1:
            plt.plot(x, y, 'x')
        else:
            plt.plot(x, y, 'o', label=str(l))
    #plt.legend()
    plt.show()


ap1 = Vector(0, 0)
ap2 = Vector(2, 1)
ap3 = Vector(1.5, 2.5)
ap4 = Vector(1, -2)
ap5 = Vector(-3, 1)

bp1 = Vector(2, 0)
bp2 = Vector(0, 1)
bp3 = Vector(-3, 1)
bp4 = Vector(-1, 0)
bp5 = Vector(-3, 0)

patha = [ap1, ap2, ap3, ap4, ap5]
pathb = [bp1, bp2, bp3, bp4, bp5]

input1 = Input(patha, pathb)
print(input1)
sample1 = input1.sample(23, 500)
print(sample1)
sample_to_matplotlib(sample1)

'''eingabe1 = TwoLineSegments(sta1, stb1)
print(eingabe1)
cell1 = eingabe1.cell(offset=Vector(10, 10), start=Vector(10, 10))
print(cell1)
# sample1 = cell1.sample_l(7, 100)  # , ((-2, 3), (-4, 6)))
sample1 = cell1.sample_l(7,
                         100)  # , rel_bounds=((0.1, 0.6), (0, 0.3)))  # ([0,2,4,6,8,10,12,14,math.sqrt(200), 15], 20)
# print(sample1)
sample_to_matplotlib(sample1)'''

'''cell2 = Cell(Vector(1.307, 1.307), Vector(-0.5412, 0.5412), Vector(0, 0), (1, 1), (0, 0.77))
print(cell2)
sample2 = cell2.sample(1, 100, ((-2, 3), (-2, 3)))
print(sample2)
sample_to_plotly(sample2)'''
