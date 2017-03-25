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

import matplotlib.pyplot as plt  # DEBUG
from Geometry import *

CellCoord = (int, int)
CM_Point = (Vector, CellCoord)
TraversalType = (float, Vector, [float], [Vector], CellCoord)  # old


# DEBUG:
# Todo: create Sample in Graphics
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


class Cell:
    def __init__(self, parallel: bool, p: LineSegment, q: LineSegment, norm_ellipsis: Ellipse,
                 bounds_xy: Bounds1D, bounds_l: (float, float), offset: Vector = Vector(0, 0)):
        self.parallel = parallel
        # line segments
        self.p = p
        self.q = q

        self.norm_ellipsis = norm_ellipsis  # normed ellipsis (l=1)
        self.bounds_xy = bounds_xy  # local cell bounds
        self.bounds_hor = Bounds1D(offset.x, offset.x + bounds_xy[0])  # global bounds horizontal
        self.bounds_ver = Bounds1D(offset.y, offset.y + bounds_xy[1])  # global bounds vertical
        self.bounds_xy_global = (offset.to_tuple(), (self.bounds_hor[1], self.bounds_ver[1]))  # global cell bounds
        self.bounds_l = bounds_l  # l length bounds
        self.offset = offset  # offset in cell-matrix

        # variables for traversal
        self.top_right = Vector(bounds_xy[0], bounds_xy[1])  # point at top right of cell local
        self.top_right_global = self.top_right + self.offset  # point at top right of cell in cell-matrix
        self.top_right_l = self.lp(self.top_right)  # l for top right of cell
        self.bottom_left = Vector(0, 0)  # point at bottom left of cell local
        self.bottom_left_global = self.offset  # point at bottom left of cell in cell-matrix
        self.bottom_left_l = self.lp(self.bottom_left_global)  # l for bottom left of cell

        # steepest decent lines l=l_ver and l'=l_hor
        if not self.parallel:  # case 1: lines are not parallel
            c = norm_ellipsis.a
            d = norm_ellipsis.b
            m = norm_ellipsis.m
            t_l_ver = math.atan(d.x / c.x)
            t_l_hor = math.atan(d.y / c.y)
            self.l_ver = LineSegment(m, self.norm_ellipsis.p(t_l_ver))
            self.l_hor = LineSegment(m, self.norm_ellipsis.p(t_l_hor))
        else:  # case 2: lines are parallel
            c = norm_ellipsis.a
            m = norm_ellipsis.m
            self.l_ver = LineSegment(m, m + c)
            self.l_hor = LineSegment(m, m + c)
        # intersection points of l and l' with cell borders local
        self.l_ver_cut_right = Vector(bounds_xy[0], self.l_ver.fx(bounds_xy[0]))
        self.l_hor_cut_top = Vector(self.l_hor.fy(bounds_xy[1]), bounds_xy[1])
        self.l_ver_cut_left = Vector(0, self.l_ver.fx(0))
        self.l_hor_cut_bottom = Vector(self.l_hor.fy(0), 0)
        # ls for these intersection points
        self.l_ver_cut_right_l = self.lp(self.l_ver_cut_right)
        self.l_hor_cut_top_l = self.lp(self.l_hor_cut_top)
        self.l_ver_cut_left_l = self.lp(self.l_ver_cut_left)
        self.l_hor_cut_bottom_l = self.lp(self.l_hor_cut_bottom)
        # intersection points in cell-matrix
        self.l_ver_cut_right_global = self.l_ver_cut_right + self.offset
        self.l_hor_cut_top_global = self.l_hor_cut_top + self.offset
        self.l_ver_cut_left_global = self.l_ver_cut_left + self.offset
        self.l_hor_cut_bottom_global = self.l_hor_cut_bottom + self.offset

        # which steepest decent traversals are possible
        self.traverses_right = not about_equal(self.l_ver_cut_right.y, self.bounds_xy[1])
        self.traverses_top = not about_equal(self.l_hor_cut_top.x, self.bounds_xy[0])
        self.traverses_top_right = self.l_ver_cut_right.y > self.bounds_xy[1] or \
                                   self.l_hor_cut_top.x > self.bounds_xy[0] or \
                                   not self.traverses_right or not self.traverses_top
        self.traverses_left = not about_equal(self.l_ver_cut_left.y, 0)
        self.traverses_bottom = not about_equal(self.l_hor_cut_bottom.x, 0)
        self.traverses_bottom_left = self.l_ver_cut_left.y < 0 \
                                     or self.l_hor_cut_bottom.x < 0 \
                                     or not self.traverses_left or not self.traverses_bottom

    def __str__(self):
        return "    Norm-" + str(self.norm_ellipsis) + '\n' + \
               "    Offset: " + str(self.offset) + '\n' + \
               "    End: " + str(self.top_right) + " l: " + str(self.top_right_l) + '\n' + \
               "    Steepest Decent Lines:\n" + \
               "      l: " + str(self.l_ver.d) + '\n' + \
               "      l': " + str(self.l_hor.d) + '\n' + \
               "    Bounds l: " + str(self.bounds_l) + '\n' + \
               "    Bounds XY: " + str(self.bounds_xy) + '\n' + \
               "    Traverses:\n" + \
               "      right: " + str(self.traverses_right) + '\n' + \
               "      top_right: " + str(self.traverses_top_right) + '\n' + \
               "      top: " + str(self.traverses_top) + '\n'

    def sample_l(self, n_l: int, n_p: int, rel_bounds: Bounds_2D = ((0, 1), (0, 1))) -> {}:
        # sample cell for n_l: # of ls, n_p: # of points per ellipsis, rel_bounds: relative xy bounds
        ls = []
        for i in range(n_l):
            l = self.bounds_l[0] + (float(i) / (n_l - 1)) * (self.bounds_l[1] - self.bounds_l[0])
            ls.append(l)

        return self.sample(ls, n_p, rel_bounds)

    def sample(self, ls: [float], n: int, rel_bounds: Bounds_2D = ((0, 1), (0, 1))) \
            -> {}:
        # sample cell with ls: array of ls to sample, n: # of points per ellipsis, rel_bounds: relative xy bounds

        bounds = ((rel_bounds[0][0] * self.bounds_xy[0], rel_bounds[0][1] * self.bounds_xy[0]),
                  (rel_bounds[1][0] * self.bounds_xy[1], rel_bounds[1][1] * self.bounds_xy[1]))

        # holds ellipses, steepest decent lines and axis in form: (name, [Vector])
        sample = {"ellipses": [], "l-lines": [], "axis": []}

        # Sample Ellipses
        for l in ls:
            if l < self.bounds_l[0] or l > self.bounds_l[1]:
                continue  # skip if l is not in bounds

            if not self.parallel:  # case 1: lines are not parallel
                if l == 0:  # plot points as points
                    p = self.norm_ellipsis.m
                    if p.in_bounds(bounds):
                        sample["ellipses"].append((l, [p + self.offset]))
                    continue

                ellipsis = self.norm_ellipsis * l

                # sample ellipsis only in bounds of cell
                ts = ellipsis.cuts_bounds_t(bounds)

                if len(ts) == 0:
                    ts.append(0)

                for i in range(len(ts)):
                    t1 = ts[i - 1]
                    t2 = ts[i]
                    if t1 >= t2:
                        t2 += 2 * math.pi

                    d_t = t2 - t1
                    mid_t = t1 + 0.5 * d_t

                    if ellipsis.p(mid_t).in_bounds(bounds):
                        ellipsis_sample = [ellipsis.p(t1) + self.offset]

                        np1 = math.ceil((t1 / (2 * math.pi)) * n)
                        np2 = math.ceil((t2 / (2 * math.pi)) * n)

                        for ip in range(np1, np2):
                            t = (ip / n) * 2 * math.pi
                            ellipsis_sample.append(ellipsis.p(t) + self.offset)

                        ellipsis_sample.append(ellipsis.p(t2) + self.offset)

                        sample["ellipses"].append((l, ellipsis_sample))

            else:  # case 2: lines are parallel
                ellipsis = self.norm_ellipsis * l
                for ps in ellipsis.cuts_bounds_p(bounds):
                    ellipsis_sample = [p + self.offset for p in ps]
                    sample["ellipses"].append((l, ellipsis_sample))

        # Sample Ellipsis axis
        if not self.parallel:
            sample["axis"].append(("c", [p + self.offset for p in
                                         LineSegment(self.norm_ellipsis.m,
                                                     self.norm_ellipsis.m + self.norm_ellipsis.a).cuts_bounds(bounds)]))
            sample["axis"].append(("d", [p + self.offset for p in
                                         LineSegment(self.norm_ellipsis.m,
                                                     self.norm_ellipsis.m + self.norm_ellipsis.b).cuts_bounds(bounds)]))

        # Sample steepest decent lines
        sample["l-lines"].append(("l", [p + self.offset for p in self.l_ver.cuts_bounds(bounds)]))
        sample["l-lines"].append(("l'", [p + self.offset for p in self.l_hor.cuts_bounds(bounds)]))

        return sample

    def lp(self, p: Vector) -> float:  # epsilon for given point
        return self.p.frl(p.x).d(self.q.frl(p.y))

    def hyperbola_horizontal(self, y: float) -> Hyperbola:  # hyperbola for set height y
        bounds = self.bounds_ver
        if y not in bounds:
            return Hyperbola.nan()
        elif y == bounds[0]:
            return self.hyperbola_bottom
        elif y == bounds[1]:
            return self.hyperbola_top
        else:
            return self.q.hyperbola_with_point(self.p.fr(y - self.offset.y))

    def hyperbola_vertical(self, x: float) -> Hyperbola:  # hyperbola for set spot x
        bounds = self.bounds_hor
        if x not in bounds:
            return Hyperbola.nan()
        elif x == bounds[0]:
            return self.hyperbola_left
        elif x == bounds[1]:
            return self.hyperbola_right
        else:
            return self.p.hyperbola_with_point(self.q.fr(x - self.offset.x))

    def free_bounds_horizontal(self, y: float, epsilon: float) -> Bounds1D:  # free interval for epsilon on height y
        ret_bounds = Bounds1D.nan()
        hyperbola = self.hyperbola_horizontal(y)
        if not hyperbola.is_nan():
            ys = hyperbola.fy(epsilon)
            if len(ys) == 1:
                ret_bounds = (ys[0], ys[0])
            elif len(ys) == 2:
                ret_bounds = (min(ys), max(ys))
        return ret_bounds

    def free_bounds_vertical(self, x: float, epsilon: float) -> Bounds1D:  # free interval for epsilon on spot x
        ret_bounds = Bounds1D.nan()
        hyperbola = self.hyperbola_vertical(x)
        if not hyperbola.is_nan():
            ys = hyperbola.fy(epsilon)
            if len(ys) == 1:
                ret_bounds = (ys[0], ys[0])
            elif len(ys) == 2:
                ret_bounds = (min(ys), max(ys))
        return ret_bounds

    def traverse_right(self, traversal: TraversalType) -> \
            TraversalType:  # steepest decent traversal to right
        max_l = math.inf
        end = self.l_ver_cut_right_global
        end_l = self.l_ver_cut_right_l

        if self.traverses_right:
            start = traversal[1]
            if self.l_ver_cut_right_global.y <= start.y:
                end = Vector(self.top_right_global.x, start.y)
                end_l = self.lp(end - self.offset)
                max_l = max(traversal[0], end_l)
            elif start.y < self.l_ver_cut_right_global.y < self.top_right_global.y:

                max_l = max(traversal[0], end_l)

        return max_l, end, traversal[2] + [end_l], traversal[3] + [end], (traversal[4][0] + 1, traversal[4][1])

    def traverse_top(self, traversal: TraversalType) -> \
            TraversalType:  # steepest decent traversal to top
        max_l = math.inf
        end = self.l_hor_cut_top_global
        end_l = self.l_hor_cut_top_l

        if self.traverses_top:
            start = traversal[1]
            if self.l_hor_cut_top_global.x <= start.x:
                end = Vector(start.x, self.top_right_global.y)
                end_l = self.lp(end - self.offset)
                max_l = max(traversal[0], end_l)
            elif start.x < self.l_hor_cut_top_global.x < self.top_right_global.x:
                max_l = max(traversal[0], end_l)

        return max_l, end, traversal[2] + [end_l], traversal[3] + [end], (traversal[4][0], traversal[4][1] + 1)

    def traverse_top_right(self, traversal: TraversalType) -> \
            TraversalType:  # steepest decent traversal to top-right
        max_l = math.inf
        end = self.top_right_global
        end_l = self.top_right_l

        if self.traverses_top_right:
            max_l = max(traversal[0], end_l)

        return max_l, end, traversal[2] + [end_l], traversal[3] + [end], (traversal[4][0] + 1, traversal[4][1] + 1)

    def traverse_top_right_force(self, traversal: TraversalType) -> \
            TraversalType:  # force traversal to top-right
        end = self.top_right_global
        end_l = self.top_right_l
        max_l = max(traversal[0], end_l)

        return max_l, end, traversal[2] + [end_l], traversal[3] + [end], traversal[4]

    def traverse(self, traversal: TraversalType) -> \
            [TraversalType]:
        # steepest decent traversal to all possible destination
        traversal_right = self.traverse_right(traversal)
        traversal_top = self.traverse_top(traversal)
        traversal_top_right = self.traverse_top_right(traversal)

        traversals = [traversal_right, traversal_top_right, traversal_top]

        return traversals


class OneLineSegment(LineSegment):  # one line segment and intersection point
    def __init__(self, ls: LineSegment, s: Vector):
        super().__init__(ls.p1, ls.p2)

        self.intersects = self.contains_point(s)  # intersects S?
        self.dir = (self.r_point(s) <= 1)  # points in direction  of S?
        self.rs = self.r_point(s)  # S in parametric terms of line segment

    def __str__(self):
        return str(self.p1) + "->" + str(self.p2) + ": l=" + str(self.l) + " m=" + str(self.m) + " n=" + str(self.n) + \
               "\n      Intersects:" + str(self.intersects) + " Direction:" + str(self.dir) + " r(S): " + str(self.rs)


class TwoLineSegments:  # calculates and saves parameters of two line segments
    def __init__(self, a, b):
        # calculate shortest and longest possible line length
        self.bounds_l = (min(a.d_ls_point(b.p1), a.d_ls_point(b.p2), b.d_ls_point(a.p1), b.d_ls_point(a.p2)),
                         max(a.p1.d(b.p1), a.p1.d(b.p2), a.p2.d(b.p1), a.p2.d(b.p2)))

        self.parallel = about_equal(a.m, b.m)
        if not self.parallel:  # case 1: lines are not parallel
            self.s = a.intersection_p(b)  # intersection point S

            self.a = OneLineSegment(a, self.s)  # line segment a
            self.b = OneLineSegment(b, self.s)  # line segment b

            self.intersect = self.a.intersects and self.b.intersects  # does the intersection point lie on A and B
            if self.intersect:  # if line segments intersect, set min length to 0
                self.bounds_l = (0.0, self.bounds_l[1])

        else:  # case 2: lines are parallel
            # parallel lines don't intersect
            self.s = Vector(math.nan, math.nan)
            self.intersect = False
            self.a = a  # line segment a
            self.b = b  # line segment b

            self.dirB = (a.r_point(a.p1 + b.d) >= 0)  # do A and B point in the same direction
            b1_p = a.project_p(b.p1)  # projection of b1 on a
            self.dist = b.p1.d(b1_p)  # distance of the two lines
            self.anchor = Vector(a.r_point(b1_p) * a.l, 0)  # point with b(0) and minimal l

    def __str__(self):
        return "    Line Segment A: " + str(self.a) + "\n    Line Segment B: " + str(self.b) + '\n' + \
               "       Parallel: " + str(self.parallel) + '\n' + \
               "       Intersect (" + str(self.intersect) + "): " + str(self.s)

    def cell(self, offset: Vector = Vector(0, 0)) -> Cell:
        # Vectors a and b normalised
        norm_a = self.a.d.norm()
        norm_b = self.b.d.norm()

        if not self.parallel:  # case 1: lines are not parallel
            # x- and y-coordinates for ellipsis vectors c and d for l = 1
            c_xy = 1 / (norm_a - norm_b).l
            d_xy = 1 / (norm_a + norm_b).l

            # Ellipsis vectors c and d for l = 1
            c = Vector(c_xy, c_xy)
            d = Vector(-d_xy, d_xy)

            # Calculate ellipsis midpoint:
            m_a = Vector(self.a.l, 0) * self.a.rs
            m_b = Vector(0, self.b.l) * self.b.rs
            m = m_a + m_b

            norm_ellipsis = Ellipse(m, c, d)

        else:  # case 2: lines are parallel
            if self.dirB:
                a = Vector(1, 1)
            else:
                a = Vector(-1, 1)

            norm_ellipsis = EllipseInfinite(self.anchor, a, self.dist)

        # set cell bounds: length of a and b
        bounds_xy = Bounds1D(self.a.d.l, self.b.d.l)

        return Cell(self.parallel, self.a, self.b, norm_ellipsis, bounds_xy, self.bounds_l, offset=offset)


# functions for handling traversals

def traverse_do(t1: TraversalType,
                t2: TraversalType) -> \
        TraversalType:  # traverse traversal t1 further by t2
    return max(t1[0], t2[0]), t2[3][-1], t1[2] + t2[2], t1[3] + t2[3], (t1[4][0] + t2[4][0], t1[4][1] + t2[4][1])


def traverse_a(t: TraversalType) -> \
        TraversalType:  # traversal: move cell indicator one in direction of A
    return t[0], t[1], t[2], t[3], (t[4][0] + 1, t[4][1])


def traverse_b(t: TraversalType) -> \
        TraversalType:  # traversal: move cell indicator one in direction of B
    return t[0], t[1], t[2], t[3], (t[4][0], t[4][1] + 1)


class Traversal:
    def __init__(self, cell_matrix: 'CellMatrix', a_cm: CM_Point, b_cm: CM_Point, path: [Vector], epsilon: float,
                 epsilons: [float]):
        self.cell_matrix = cell_matrix
        
        self.a_cm = a_cm
        self.a = a_cm[0]
        self.cell_a = a_cm[1]
        self.b_cm = b_cm
        self.b = b_cm[0]
        self.cell_b = b_cm[1]
        
        self.path = path
        self.epsilon = epsilon
        self.epsilons = epsilons

    def __str__(self):
        return "    " + str(self.epsilon) + "-Traversal:" + '\n' + \
               "      A: " + str(self.a) + " -> B: " + str(self.b) + '\n' + \
               "      Cell_A: " + str(self.cell_a) + " -> Cell_B: " + str(self.cell_b) + '\n' + \
               "      Path: " + str([str(point) for point in self.path]) + '\n' + \
               "      Epsilon: " + str(self.epsilon) + '\n' + \
               "      Epsilons: " + str(self.epsilons)


class CriticalEvents:
    def __init__(self):
        self.dict = {}

    def __str__(self):
        desc = ""
        for traversal in self.list():
            desc += '\n' + str(traversal)
        return desc

    def __getitem__(self, item) -> [Traversal]:
        if item not in self.dict:
            return []
        return self.dict[item]

    def append(self, traversal: Traversal):
        epsilon = traversal.epsilon
        if epsilon not in self.dict:
            self.dict[epsilon] = []
        self.dict[epsilon].append(traversal)

    def list(self) -> [Traversal]:
        sorted_events = []
        for epsilon in self.epsilons():
            sorted_events += self[epsilon]
        return sorted_events

    def epsilons(self) -> [float]:
        return sorted(self.dict.keys())


class CrossSection:
    def __init__(self, path: Path, point: Vector):
        self.path = path
        self.point = point

        self.hyperbolas = []
        for i_segment in range(len(path.segments)):
            segment = path.segments[i_segment]
            hyperbola = segment.hyperbola_with_point(point).move_x(path.offsets[i_segment])
            self.hyperbolas.append(hyperbola)

        self.hyperbolas_overload = [self.hyperbolas[0].reflect_x(self.path.offsets[0])] + self.hyperbolas + \
                                   [self.hyperbolas[-1].reflect_x(self.path.offsets[-1])]
        self._minima = []
        self._minima_no_borders = []
        self._minima_borders = []
        self._maxima = []

    def __str__(self):
        desc = ""
        for i in range(self.path.count):
            bounds = (self.path.offsets[i], self.path.offsets[i + 1])
            desc += "     " + str(i) + ": " + str(bounds) + ": " + str(self.hyperbolas[i]) + '\n'
        return desc

    def __len__(self):
        return self.path.count

    def __getitem__(self, item) -> Hyperbola:
        assert 0 <= int(item) <= self.path.count, "Error: No Hyperbola with index: " + str(item) + '\n' + \
                                                  "Hyperbolas: \n" + str(self)
        return self.hyperbolas[item]

    def minima(self) -> [float]:  # returns all local minima
        if len(self._minima) == 0:
            self._minima = self.minima_no_borders() + self.minima_borders()
        return self._minima

    def minima_no_borders(self) -> [float]:  # returns local minima that don't lie on borders
        if len(self._minima_no_borders) == 0:
            for i in range(self.path.count):
                hyperbola = self.hyperbolas[i]
                bounds = (self.path.offsets[i], self.path.offsets[i + 1])
                x = hyperbola.s.x
                if bounds[0] < x < bounds[1]:
                    self._minima_no_borders.append(x)
        return self._minima_no_borders

    def minima_borders(self) -> [float]:  # returns local minima that lie on borders
        if len(self._minima_borders) == 0:
            for i in range(self.path.count + 1):
                x = self.path.offsets[i]
                left_hyperbola = self.hyperbolas_overload[i]
                right_hyperbola = self.hyperbolas_overload[i + 1]
                if left_hyperbola.orientation(x) <= 0 <= right_hyperbola.orientation(x):
                    self._minima_borders.append(x)
        return self._minima_borders

    def is_minima(self, x: float) -> bool:  # checks if local minima at set x
        if len(self._minima) == 0:
            self._minima = self.minima()
        return x in self._minima

    def maxima(self) -> [float]:  # returns all local maxima (always on borders)
        if len(self._maxima) == 0:
            for i in range(self.path.count + 1):
                x = self.path.offsets[i]
                left_hyperbola = self.hyperbolas_overload[i]
                right_hyperbola = self.hyperbolas_overload[i + 1]
                if left_hyperbola.orientation(x) > 0 > right_hyperbola.orientation(x):
                    self._maxima.append(x)
        return self._maxima

    def is_maxima(self, x: float) -> bool:  # checks if a point is a local maxima
        if len(self._maxima) == 0:
            self._maxima = self.maxima()
        return x in self._maxima


class CellMatrix:
    def __init__(self, points_p: [Vector], points_q: [Vector], traverse: int = 0, fast: bool = False):
        # paths
        self.p = Path(points_p)
        self.q = Path(points_q)

        # real bounds of length l over the whole cell matrix
        self.bounds_l = (math.inf, 0)

        # border hyperbolas
        self.cross_sections_hor = self.calculate_cross_sections(self.p, self.q.points)
        self.cross_sections_ver = self.calculate_cross_sections(self.q, self.p.points)

        # generate and save TwoLineSegments & Cells
        self.twoLSs: [[TwoLineSegments]] = []
        self.cells: [[Cell]] = []
        for i_p in range(self.p.count):
            self.twoLSs.append([])
            self.cells.append([])
            for i_q in range(self.q.count):
                two_line_segments = TwoLineSegments(self.p.segments[i_p], self.q.segments[i_q])
                self.bounds_l = (min(self.bounds_l[0], two_line_segments.bounds_l[0]),
                                 max(self.bounds_l[1], two_line_segments.bounds_l[1]))
                self.twoLSs[i_p].append(two_line_segments)
                cell = two_line_segments.cell(offset=Vector(self.p.offsets[i_p], self.q.offsets[i_q]))
                cell.hyperbola_left = self.cross_sections_ver[i_p][i_q]
                cell.hyperbola_bottom = self.cross_sections_hor[i_q][i_p]
                cell.hyperbola_right = self.cross_sections_ver[i_p + 1][i_q]
                cell.hyperbola_top = self.cross_sections_hor[i_q + 1][i_p]
                self.cells[i_p].append(cell)

        # critical events
        self.critical_events_hor = self.calculate_critical_events(horizontal=True)
        self.critical_events_ver = self.calculate_critical_events(horizontal=False)

        # DEBUG: sample hyperbolas:
        fig_both = plt.figure(figsize=plt.figaspect(0.5))
        ax_hor = fig_both.add_subplot(2, 1, 1, aspect=1, ylim=self.bounds_l, xlabel="p", ylabel="ε")
        ax_ver = fig_both.add_subplot(2, 1, 2, aspect=1, ylim=self.bounds_l, xlabel="p", ylabel="ε")
        # horizontal
        fig_hor = plt.figure(figsize=plt.figaspect(0.5))
        for i_q in range(self.q.count + 1):
            ax = fig_hor.add_subplot(self.q.count + 1, 1, self.q.count - i_q + 1, aspect=1, ylim=self.bounds_l,
                                     xlabel="p", ylabel="ε")
            all_points = []
            for i_p in range(self.p.count):
                points = []
                bounds = (self.p.offsets[i_p], self.p.offsets[i_p + 1])
                n_points = math.ceil(100 * (self.p.lengths[i_p] / self.p.length))
                sample = self.cross_sections_hor[i_q][i_p].sample(bounds, n_points)
                points += sample
                all_points += sample
                x, y = vectors_to_xy(points)
                ax.plot(x, y)
            ax_hor.plot(*vectors_to_xy(all_points), label="q = " + str(i_q))
        # vertical
        fig_ver = plt.figure(figsize=plt.figaspect(0.5))
        for i_p in range(self.p.count + 1):
            ax = fig_ver.add_subplot(1, self.p.count + 1, i_p + 1, aspect=1, xlim=self.bounds_l, xlabel="ε", ylabel="p")
            all_points = []
            for i_q in range(self.q.count):
                points = []
                bounds = (self.q.offsets[i_q], self.q.offsets[i_q + 1])
                n_points = math.ceil(100 * (self.q.lengths[i_q] / self.q.length))
                sample = self.cross_sections_ver[i_p][i_q].sample(bounds, n_points)
                points += sample
                all_points += sample
                x, y = vectors_to_xy(points)
                ax.plot(y, x)
            ax_ver.plot(*vectors_to_xy(all_points), label="p = " + str(i_p))
        ax_hor.legend()
        ax_ver.legend()

        # set fast=True to quit when global minimum is reached, local minimum is disregarded
        self.fast = fast
        self.global_minimum_reached = False
        # lower limit for smallest globally reachable l
        self.min_global_l = max(self.p.points[0].d(self.q.points[0]), self.p.points[-1].d(self.q.points[-1]))

        # traverse
        self.traverse = traverse

        # Old Traversal Call:
        # traverse
        self.traverse = traverse
        if self.traverse > 0:
            self.lowest_l = math.inf
            self.traversals = self.traverse_best(self.traverse, delete_duplicates=True)

    def __str__(self):
        desc = "Input (" + str(self.p.count) + "x" + str(self.q.count) + ") (" + \
               str(self.p.length) + "x" + str(self.q.length) + "):\n"
        desc += " Path P " + str(self.p)
        desc += " Path Q " + str(self.q)
        desc += '\n'
        desc += "==>\n"
        desc += " Two Line Segments & Cells:\n"
        for i_p in range(self.p.count):
            for i_q in range(self.q.count):
                desc += " Cell " + str(i_p) + "x" + str(i_q) + '\n'
                desc += str(self.twoLSs[i_p][i_q]) + '\n'
                desc += str(self.cells[i_p][i_q]) + '\n'
        desc += '\n'
        desc += " Bounds_l: " + str(self.bounds_l) + '\n'
        desc += '\n'
        desc += " Border Hyperbolas: " + '\n'
        desc += "  Horizontal: " + '\n'
        for i_q in range(self.q.count + 1):
            desc += "   " + str(i_q) + ".\n" + str(self.cross_sections_hor[i_q]) + '\n'
        desc += "  Vertical: " + '\n'
        for i_p in range(self.p.count + 1):
            desc += "   " + str(i_p) + ".\n" + str(self.cross_sections_ver[i_p]) + '\n'
        desc += '\n'
        desc += " Critical Events: " + '\n'
        desc += "  Horizontal: " + str(self.critical_events_hor) + '\n'
        desc += "  Vertical: " + str(self.critical_events_ver) + '\n'
        desc += '\n'
        if self.traverse > 0:
            desc += " Lowest_l: " + str(self.lowest_l) + '\n'
            desc += " Traversals:\n"
            for traversal in self.traversals:
                desc += "   " + str(traversal) + '\n'

        return desc

    @staticmethod
    def calculate_cross_sections(path, points) -> [CrossSection]:
        cross_sections = []

        for point in points:
            cross_section = CrossSection(path, point)
            cross_sections.append(cross_section)

        return cross_sections

    @staticmethod
    def calculate_critical_points(cross_sections: [CrossSection], other_cross_sections: [CrossSection]) -> [[Vector]]:
        critical_points = []
        n_cross_sections = len(cross_sections)
        path = cross_sections[0].path
        other_path = other_cross_sections[0].path

        # events of type a
        for i_cross_section in range(1, len(cross_sections) - 1):
            cross_section = cross_sections[i_cross_section]

            # no borders
            minima_no_borders = cross_section.minima_no_borders()
            for x in minima_no_borders:
                other_cross_section = CrossSection(other_path, path.p_rl(x))
                y = other_path.offsets[i_cross_section]
                if other_cross_section.is_maxima(y):
                    critical_points.append([Vector(x, y)])

            # borders
            minima_borders = cross_section.minima_borders()
            for x in minima_borders:
                other_cross_section = other_cross_sections[path.i_rl(x)]
                y = other_path.offsets[i_cross_section]
                if other_cross_section.is_maxima(y):
                    critical_points.append([Vector(x, y)])

        # events of type b
        for i_section in range(path.count):
            bounds = (path.offsets[i_section], path.offsets[i_section + 1])
            for i_start in range(n_cross_sections):
                cross_section_start = cross_sections[i_start]
                hyperbola_start = cross_section_start[i_section]
                for i_end in range(n_cross_sections - 1, i_start, -1):
                    cross_section_end = cross_sections[i_end]
                    hyperbola_end = cross_section_end[i_section]

                    intersection = hyperbola_start.intersects_hyperbola_in_bounds_critical(hyperbola_end, bounds)
                    x = intersection.x
                    if math.isnan(x):
                        continue

                    other_cross_section = CrossSection(other_path, path.p_rl(x))
                    if not other_cross_section.is_maxima(other_path.offsets[i_start]) or \
                            not other_cross_section.is_maxima(other_path.offsets[i_end]):
                        continue

                    is_critical = True
                    i_between = i_start + 1
                    points = [Vector(x, other_path.offsets[i_start])]
                    while is_critical and i_between < i_end:
                        cross_section_between = cross_sections[i_between]
                        hyperbola_between = cross_section_between[i_section]
                        points.append((Vector(x, other_path.offsets[i_between])))
                        is_critical = is_critical and intersection.y >= hyperbola_between.fx(intersection.x)
                        i_between += 1
                    points.append(Vector(x, other_path.offsets[i_end]))

                    if is_critical:
                        critical_points.append(points)

        return critical_points

    def calculate_critical_events(self, horizontal: bool) -> CriticalEvents:
        critical_events = CriticalEvents()
        critical_points = []

        if horizontal:
            critical_points += self.calculate_critical_points(self.cross_sections_ver, self.cross_sections_hor)
        else:
            critical_points += self.calculate_critical_points(self.cross_sections_hor, self.cross_sections_ver)

        for points in critical_points:
            if horizontal:
                points = [point.x_to_y() for point in points]

            traversal = self.traversal_from_points(points)
            critical_events.append(traversal)

        return critical_events

    def decide_critical_traversal(self, a1_cm: CM_Point, traversal: Traversal, b2_cm: CM_Point) -> bool:
        epsilon = traversal.epsilon
        b1_cm = traversal.a_cm
        a2_cm = traversal.b_cm
        decision1 = self.decide_traversal(a1_cm, b1_cm, epsilon)
        decision2 = self.decide_traversal(a2_cm, b2_cm, epsilon)
        return decision1 and decision2
    
    def decide_traversal(self, a_cm: CM_Point, b_cm: CM_Point, epsilon: float) -> bool:
        a = a_cm[0]
        cell_a = a_cm[1]
        b = b_cm[0]
        cell_b = b_cm[1]
        
        if not a < b:
            print("Error: Traversal is impossible, because not a < b: a=" + str(a) + " and b=" + str(b))
            return False
        
        range_p = range(cell_a[0], cell_b[0] + 1)
        range_q = range(cell_a[1], cell_b[1] + 1)

        start_cell = self.cells[cell_a[0]][cell_a[1]]
        start_bounds_p = Bounds1D(a.x, start_cell.bounds_x[1])
        start_bounds_q = Bounds1D(a.y, start_cell.bounds_y[1])
        reachable_bounds_p = start_cell.free_bounds_horizontal(a.y, epsilon)
        reachable_bounds_q = start_cell.free_bounds_vertical(a.x, epsilon)
        reachable_p = start_bounds_p.cut(reachable_bounds_p)
        reachable_q = start_bounds_q.cut(reachable_bounds_q)
        
        for i_p in range_p:
            reachable_p.append([])
            reachable_q.append([])
            for i_q in range_q:
                print("fghjk")  #TODO: next!!!

        return True

    def traversal_from_points(self, points) -> Traversal:
        start = points[0]
        end = points[-1]

        cell_a = (self.p.i_rl(start.x), self.q.i_rl(start.y))
        cell_b = (self.p.i_rl(end.x), self.q.i_rl(end.y))

        epsilons = []
        for point in points:
            tmp_epsilon = self.p.p_rl(point.x).d(self.q.p_rl(point.y))
            epsilons.append(tmp_epsilon)
        epsilon = max(epsilons)

        return Traversal(self, (start, cell_a), (end, cell_b), points, epsilon, epsilons)

    # Old Traversal (v2):
    def traverse_best(self, criteria: int = 2, delete_duplicates: bool = True, start: Vector = Vector(0, 0)) \
            -> [TraversalType]:  # traverses cell-matrix and chooses best traversal depending on criteria
        i_p = 0
        i_q = 0
        while start.x > self.p.offsets[i_p + 1]:
            i_p += 1
        while start.y > self.q.offsets[i_q + 1]:
            i_q += 1

        print("Lower Limit for minimum global l: " + str(self.min_global_l))

        l_start = self.cells[i_p][i_q].lp(start)
        traversals0 = self.traverse_rec((l_start, start, [l_start], [start], (i_p, i_q)))
        traversals = traversals0

        # select best traversal(s):
        if criteria > 1:
            # 1. lowest l
            lowest_l = math.inf
            for traversal in traversals0:
                lowest_l = min(lowest_l, traversal[0])
            traversals1 = []
            for traversal in traversals0:
                if traversal[0] <= lowest_l or about_equal(traversal[0], lowest_l):
                    traversals1.append(traversal)
            traversals = traversals1

            if criteria > 2:
                # 2. lowest average of ls
                lowest_avg_ls = math.inf
                for traversal in traversals1:
                    avg_ls = sum(traversal[2]) / len(traversal[2])
                    lowest_avg_ls = min(lowest_avg_ls, avg_ls)
                traversals2 = []
                for traversal in traversals1:
                    avg_ls = sum(traversal[2]) / len(traversal[2])
                    if avg_ls <= lowest_avg_ls or about_equal(avg_ls, lowest_avg_ls):
                        traversals2.append(traversal)
                traversals = traversals2

        traversals_set = []
        if len(traversals) > 1 and delete_duplicates:
            for i1 in range(len(traversals)):
                unique = True
                traversal1 = traversals[i1]
                for i2 in range(i1 + 1, len(traversals)):
                    traversal2 = traversals[i2]
                    if traversal1[3] == traversal2[3]:
                        unique = False
                        break
                if unique:
                    traversals_set.append(traversal1)
            return traversals_set

        return traversals

    # Old Traversal rec
    def traverse_rec(self, traversal: TraversalType) -> \
            [TraversalType]:  # traverse cells recursive
        # indicate which cell to traverse next
        i_p = traversal[4][0]
        i_q = traversal[4][1]

        if self.global_minimum_reached and self.fast:
            return []

        # complete recursive call on reach of cell-matrix top right
        if i_p >= self.p.count and i_q >= self.q.count:
            if traversal[0] < self.lowest_l:
                self.lowest_l = traversal[0]
                print("momentary lowest_l = " + str(self.lowest_l))
                if self.lowest_l <= self.min_global_l and not self.global_minimum_reached:
                    print("global minimum reached: " + str(self.lowest_l))
                    self.global_minimum_reached = True
            return [traversal]

        # stores possible traversals of this cell and beyond
        traversals = []

        # if traversal hits top- or right-side of cell-matrix move to top-right
        if i_p >= self.p.count or i_q >= self.q.count:
            if i_p < self.p.count:
                cell = self.cells[i_p][self.q.count - 1]
                if cell.top_right_l <= self.lowest_l + tol:
                    next_traversal = traverse_a(cell.traverse_top_right_force(traversal))
                    traversals += self.traverse_rec(next_traversal)
            if i_q < self.q.count:
                cell = self.cells[self.p.count - 1][i_q]
                if cell.top_right_l <= self.lowest_l + tol:
                    next_traversal = traverse_b(cell.traverse_top_right_force(traversal))
                    traversals += self.traverse_rec(next_traversal)

        else:  # if traversal is not at the top- or right-side yet

            # traverse cell by steepest decent
            sd_traversal = self.cells[i_p][i_q].traverse(traversal)
            # traverse to right side?
            if self.lowest_l + tol >= sd_traversal[0][0] and not math.isinf(sd_traversal[0][0]):
                traversals += self.traverse_rec(sd_traversal[0])
            # traverse to top-right side?
            if self.lowest_l + tol >= sd_traversal[1][0] and not math.isinf(sd_traversal[1][0]):
                traversals += self.traverse_rec(sd_traversal[1])
            # traverse to top side?
            if self.lowest_l + tol >= sd_traversal[2][0] and not math.isinf(sd_traversal[2][0]):
                traversals += self.traverse_rec(sd_traversal[2])

            # traverse by critical traversal paths
            p = traversal[1]  # last point of traversal
            # horizontal
            if i_p < self.p.count - 1:
                for critical_traversal_horizontal in self.critical_events_hor[i_p + 1][i_q]:
                    max_l = critical_traversal_horizontal[0]
                    start = critical_traversal_horizontal[1]
                    if p < start and self.lowest_l + tol >= max_l and \
                            (not self.fast or max_l >= self.min_global_l - tol):
                        next_traversal = traverse_a(traverse_do(traversal, critical_traversal_horizontal))
                        traversals += self.traverse_rec(next_traversal)
            # vertical
            if i_q < self.q.count - 1:
                for critical_traversal_vertical in self.critical_events_ver[i_p][i_q + 1]:
                    max_l = critical_traversal_vertical[0]
                    start = critical_traversal_vertical[1]
                    if p < start and self.lowest_l + tol >= max_l and \
                            (not self.fast or max_l >= self.min_global_l - tol):
                        next_traversal = traverse_b(traverse_do(traversal, critical_traversal_vertical))
                        traversals += self.traverse_rec(next_traversal)

        return traversals

    def sample_l(self, n_l: int, n_p: int, heatmap: int = 500, traversals_n: int = 10) -> {}:  # DEBUG: heatmap = 100
        # sample cell-matrix with n_l: number of ls and n_p: points per ellipses
        ls = []

        if n_l > 0:
            for i in range(n_l + 1):
                l = self.bounds_l[0] + (float(i) / n_l) * (self.bounds_l[1] - self.bounds_l[0])
                ls.append(l)

        return self.sample(ls, n_p, heatmap=heatmap, traversals_n=traversals_n)

    def sample(self, ls: [float], n_p: int, heatmap: int = 100, traversals_n: int = 10) -> {}:
        # sample cell-matrix for given ls and n_p: points per ellipses
        samples = {"bounds-l": [], "borders-v": [], "borders-h": [], "cells": [], "traversals": []}

        # are all ls in bounds
        for l in ls:
            if l < self.bounds_l[0] or l > self.bounds_l[1]:
                print("l: " + str(l) + " is not in bounds_l: " + str(self.bounds_l))
                ls.remove(l)

        # include bounds_l in sample
        samples["bounds-l"] = [self.bounds_l[0], self.bounds_l[1]]

        # sample input
        samples["input"] = [self.p.points, self.q.points]

        # sample cells
        for i_p in range(self.p.count):
            for i_q in range(self.q.count):
                cell = self.cells[i_p][i_q]
                samples["cells"].append((str(i_p) + "x" + str(i_q), cell.sample(ls, n_p)))

        # sample cell borders
        for i in range(1, self.p.count):  # vertical
            samples["borders-v"].append(("border-v: " + str(i), [Vector(self.p.offsets[i], 0),
                                                                 Vector(self.p.offsets[i], self.q.length)]))
        for i in range(1, self.q.count):  # horizontal
            samples["borders-h"].append(("border-h: " + str(i), [Vector(0, self.q.offsets[i]),
                                                                 Vector(self.p.length, self.q.offsets[i])]))

        # include size in sample
        samples["size"] = (self.p.length, self.q.length)

        # sample critical traversals
        samples["critical-traversals"] = self.critical_events_hor.list() + self.critical_events_ver.list()

        # Old Traversal Sample:
        # sample traversals
        if self.traverse > 0:
            for traversal in self.traversals:
                if traversal[0] != -1:
                    samples["traversals"].append(traversal)

        # sample a traversal
        if self.traverse > 0:
            traversal = self.traversals[0]
            samples["traversal"] = self.sample_traversal(traversal, traversals_n * max(self.p.count, self.q.count))

        # sample heatmap
        if heatmap > 0:
            samples["heatmap"] = self.sample_heatmap_a(heatmap)

        return samples

    def sample_heatmap_a(self, n_a: int) -> []:  # sample heat map with squares scaled by n_a divisions on a axis
        n_b = math.floor(n_a * (self.q.length / self.p.length))
        return self.sample_heatmap(n_a, int(n_b))

    def sample_heatmap(self, n_a: int, n_b: int) -> []:  # sample heatmap by n_a rectangles on a-axis and n_b on b-axis
        # coordinates
        xs = [[]]
        ys = [[]]
        zs = [[]]
        # step size on a and b
        s_a = self.p.length / n_a
        s_b = self.q.length / n_b
        # x- & y-coordinate to iterate through
        x, y = 0, 0
        # counters of x- & and y-coordinates
        i_x, i_y = 0, 0
        # counters for active cell
        c_a, c_b = 0, 0

        # iterate through cells & a-/b-axis
        while c_b < self.q.count:
            while y <= self.q.offsets[c_b + 1] or (c_b >= self.q.count - 1 and i_y <= n_b):
                while x <= self.p.offsets[c_a + 1] or (c_a >= self.p.count - 1 and i_x <= n_a):
                    xs[-1].append(x)
                    ys[-1].append(y)
                    z = self.cells[c_a][c_b].lp(Vector(x - self.p.offsets[c_a], y - self.q.offsets[c_b]))
                    zs[-1].append(z)

                    i_x += 1
                    x += s_a
                c_a += 1

                if c_a >= self.p.count and c_b < self.q.count:
                    xs.append([])
                    ys.append([])
                    zs.append([])
                    x = 0
                    c_a = 0
                    i_x = 0
                    y += s_b
                    i_y += 1
            c_b += 1

        del xs[-1]
        del ys[-1]
        del zs[-1]

        return [xs, ys, zs]

    # Old Traversal Sampling
    def points_for_traversal_point(self, traversal_p: Vector, c_a: int = 0, c_b: int = 0) -> [Vector, Vector]:

        r_a = traversal_p.x
        r_b = traversal_p.y

        if r_a < self.p.offsets[c_a]:
            c_a = 0
        if r_b < self.q.offsets[c_b]:
            c_b = 0

        while r_a > self.p.offsets[c_a + 1] and c_a < self.p.count - 1:
            c_a += 1
        while r_b > self.q.offsets[c_b + 1] and c_b < self.q.count - 1:
            c_b += 1

        r_a -= self.p.offsets[c_a]
        r_b -= self.q.offsets[c_b]
        a = self.p.segments[c_a]
        b = self.q.segments[c_b]
        pa = a.frl(r_a)
        pb = b.frl(r_b)

        return [pa, pb]

    # Old Traversal Sampling
    def sample_traversal(self, traversal: TraversalType, n: int) -> {}:
        # sample a specific traversal for lines between paths and 3d-plot with n points

        sample = {"in-traversal-l": [], "in-traversal": [], "traversal-3d": [], "traversal-3d-l": []}

        max_l = traversal[0]
        ls = traversal[2]
        traversal_ls = []
        traversal_length = 0

        for i in range(1, len(traversal[3])):
            p1 = traversal[3][i - 1]
            p2 = traversal[3][i]
            if p1 != p2:
                linesegment = LineSegment(p1, p2)
                traversal_ls.append(linesegment)
                traversal_length += linesegment.l

        x, y, z = [], [], []  # arrays for traversal 3d-plot
        x_l, y_l, z_l = [], [], []  # arrays for traversal 3d-plot max_l

        c_a, c_b = 0, 0  # counters for active cell

        for c_t in range(len(traversal_ls)):

            t_ls = traversal_ls[c_t]
            n_t = math.ceil((t_ls.l / traversal_length) * n)
            i_t = 0

            p1 = t_ls.p1
            while p1.x > self.p.offsets[c_a + 1] and c_a < self.p.count - 1:
                c_a += 1
            while p1.y > self.q.offsets[c_b + 1] and c_b < self.q.count - 1:
                c_b += 1

            if ls[c_t] >= max_l:
                line = self.points_for_traversal_point(p1, c_a=c_a, c_b=c_b)
                sample["in-traversal-l"].append(line)

                x_l.append(p1.x)
                y_l.append(p1.y)
                z_l.append(ls[c_t])

                x.append(p1.x)
                y.append(p1.y)
                z.append(ls[c_t])

                i_t += 1

            while i_t < n_t:
                t = i_t / n_t
                t_p = t_ls.fr(t)

                line = self.points_for_traversal_point(t_p, c_a=c_a, c_b=c_b)
                sample["in-traversal"].append(line)

                x.append(t_p.x)
                y.append(t_p.y)
                z.append(line[0].d(line[1]))

                i_t += 1

        if ls[-1] >= max_l:
            sample["in-traversal-l"].append([self.p.points[-1], self.q.points[-1]])
            x_l.append(self.p.length)
            y_l.append(self.q.length)
            z_l.append(ls[-1])
        else:
            sample["in-traversal"].append([self.p.points[-1], self.q.points[-1]])
        x.append(self.p.length)
        y.append(self.q.length)
        z.append(ls[-1])

        sample["traversal-3d"] = [x, y, z]
        sample["traversal-3d-l"] = [x_l, y_l, z_l]

        return sample
