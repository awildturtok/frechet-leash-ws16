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
        self.bounds_hor = Bounds1D(offset.x, offset.x + bounds_xy[0])  # global bounds horizontal
        self.bounds_ver = Bounds1D(offset.y, offset.y + bounds_xy[1])  # global bounds vertical
        self.bounds_xy = (self.bounds_hor, self.bounds_ver)  # global cell bounds
        self.bounds_l = bounds_l  # l length bounds
        self.offset = offset  # offset in cell-matrix

        # border hyperbolas:
        self.hyperbola_bottom = p.hyperbola_with_point(q.p1).move_x(offset.x)
        self.hyperbola_top = p.hyperbola_with_point(q.p2).move_x(offset.x)
        self.hyperbola_left = q.hyperbola_with_point(p.p1).move_x(offset.y)
        self.hyperbola_right = q.hyperbola_with_point(p.p2).move_x(offset.y)

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
        self.acute = self.l_hor.d.orientation() and self.l_ver.d.orientation()

    def __str__(self):
        return "    Offset: " + str(self.offset) + '\n' + \
               "    Norm-" + str(self.norm_ellipsis) + '\n' + \
               "    Steepest Decent Lines:\n" + \
               "      l: " + str(self.l_ver) + '\n' + \
               "      l': " + str(self.l_hor) + '\n' + \
               "    Bounds l: " + str(self.bounds_l) + '\n' + \
               "    Bounds XY: " + str(self.bounds_xy) + '\n'

    def lp(self, p: Vector) -> float:  # epsilon for given point
        return self.p.frl(p.x).d(self.q.frl(p.y))

    def steepest_decent(self, a: Vector, direction: int) -> (Vector, float, (Hyperbola, Hyperbola)):
        # calculates data needed for steepest decent to the top-right
        assert a.x in self.bounds_hor and a.y in self.bounds_ver, "Error: Cannot do Steepest Decent.\n" +\
                                                                  "A:" + str(a) + " not in Cell:\n" + str(self)
        assert direction == 1 or direction == -1, "Error: Cannot do Steepest Decent.\n" +\
                                                  "Invalid direction (" + str(direction) + ") given."

        # define default return values
        a_epsilon = self.lp(a - self.offset)
        hyperbola_hor = Hyperbola.nan()
        hyperbola_ver = Hyperbola.nan()

        # steepest decent lines
        lv = self.l_ver
        lh = self.l_hor

        # which case for steepest decent
        if direction == 1:
            case_hor = lv.point_above(a)
            case_ver = lh.point_right(a)
            lh_left = lh.point_left(a)
            lv_below = lv.point_below(a)
        else:
            case_hor = lv.point_below(a)
            case_ver = lh.point_left(a)
            lh_left = lh.point_right(a)
            lv_below = lv.point_above(a)
        case_eq = lh_left and lv_below
        on_lv = lv.point_on(a)
        on_lh = lh.point_on(a)
        on_l = on_lv or on_lh

        # determine direction of traversal
        if not self.parallel and self.norm_ellipsis.m == a:
            print("Error: Cannot do Steepest Decent. Lowest Point reached.\nA:" + str(a) + " Cell:\n" + str(self))
            return a, a_epsilon, (hyperbola_hor, hyperbola_ver)
        elif on_lv:
            if self.acute:
                r = lv.d.norm_dir() * direction
            else:
                r = Vector(direction, 0)
        elif on_lh:
            if self.acute:
                r = lh.d.norm_dir() * direction
            else:
                r = Vector(0, direction)
        elif case_hor:
            assert lh_left, "Error: Cannot do Steepest Decent.\n" +\
                            "A:" + str(a) + " above l & right of l'. Cell:\n" + str(self)
            r = Vector(direction, 0)
        elif case_ver:
            assert lv_below, "Error: Cannot do Steepest Decent.\n" +\
                             "A:" + str(a) + " above l & right of l'. Cell:\n" + str(self)
            r = Vector(0, direction)
        elif case_eq:
            r = Vector(math.sqrt(0.5), math.sqrt(0.5)) * direction
        else:
            print("Error: no direction for steepest decent found!")
            r = Vector(direction, direction)

        # determine endpoint a2 of traversal
        a_r = LineSegment(a, a + r)

        if not on_l:
            lv_rl = a_r.intersection_rl(lv)
            lh_rl = a_r.intersection_rl(lh)
            cut_l_rl = np.nanmin([lv_rl, lh_rl])
        elif not self.parallel and (direction == 1 and a < self.norm_ellipsis.m) or\
                (direction == -1 and a > self.norm_ellipsis.m):
            cut_l_rl = a.d(self.norm_ellipsis.m)
        else:
            cut_l_rl = math.inf
        '''elif not self.acute:
            bisector = LineSegment(self.norm_ellipsis.m, self.norm_ellipsis.m + Vector(-1, 1))
            cut_l_rl = a_r.intersection_rl(bisector)
            if about_equal(cut_l_rl, 0):
                cut_l_rl = math.inf
        else:
            cut_l_rl = math.inf'''  # DEBUG

        cut_hor_rl = a_r.rly(self.bounds_ver[direction])
        cut_ver_rl = a_r.rlx(self.bounds_hor[direction])
        cut_border_rl = np.nanmin([cut_hor_rl, cut_ver_rl])

        if cut_border_rl <= cut_l_rl:
            a2_rl = cut_border_rl
        else:
            a2_rl = cut_l_rl

        a2 = a_r.frl(a2_rl)
        #a2 = a2.to_bounds(self.bounds_xy)

        if a == a2:
            print("Error: Steepest decent not possible.")
            return a, a_epsilon, (hyperbola_hor, hyperbola_ver)

        a2_epsilon = self.lp(a2 - self.offset)
        a2_epsilon = self.bounds_l.to_bounds(a2_epsilon)

        # determine horizontal and vertical hyperbolas
        d_x = abs(a2.x - a.x)
        d_y = abs(a2.y - a.y)
        a_a2 = LineSegment(a, a2)
        if about_equal(d_x, 0):
            hyperbola_ver = self.hyperbola_vertical(a2.x)
        elif about_equal(d_y, 0):
            hyperbola_hor = self.hyperbola_horizontal(a2.y)
        else:
            p_cut = LineSegment(self.p.frl(min(a.x, a2.x) - self.offset.x), self.p.frl(max(a.x, a2.x) - self.offset.x))
            q_cut = LineSegment(self.q.frl(min(a.y, a2.y) - self.offset.y), self.q.frl(max(a.y, a2.y) - self.offset.y))
            hyperbola = p_cut.hyperbola_with_line(q_cut)
            hyperbola_hor = hyperbola.scaled(d_x / a_a2.l).move_x(min(a.x, a2.x))
            hyperbola_ver = hyperbola.scaled(d_y / a_a2.l).move_x(min(a.y, a2.y))

        return a2, a2_epsilon, (hyperbola_hor, hyperbola_ver)

    def hyperbola_horizontal(self, y: float) -> Hyperbola:  # hyperbola for set height y
        bounds = self.bounds_ver
        if y not in bounds:
            return Hyperbola.nan()
        elif about_equal(y, bounds[0]):
            return self.hyperbola_bottom
        elif about_equal(y, bounds[1]):
            return self.hyperbola_top
        else:
            return self.q.hyperbola_with_point(self.p.frl(y - self.offset.y)).move_x(self.offset.x)

    def hyperbola_vertical(self, x: float) -> Hyperbola:  # hyperbola for set spot x
        bounds = self.bounds_hor
        if x not in bounds:
            return Hyperbola.nan()
        elif about_equal(x, bounds[0]):
            return self.hyperbola_left
        elif about_equal(x, bounds[1]):
            return self.hyperbola_right
        else:
            return self.p.hyperbola_with_point(self.q.frl(x - self.offset.x)).move_x(self.offset.y)

    def free_bounds_horizontal(self, y: float, epsilon: float) -> Bounds1D:  # free interval for epsilon on height y
        ret_bounds = Bounds1D.nan()
        hyperbola = self.hyperbola_horizontal(y)
        if not hyperbola.is_nan():
            ys = hyperbola.fy(epsilon)
            if len(ys) == 1:
                ret_bounds = Bounds1D(ys[0], ys[0])
            elif len(ys) == 2:
                ret_bounds = Bounds1D(min(ys), max(ys))
        return self.bounds_hor.cut(ret_bounds)

    def free_bounds_vertical(self, x: float, epsilon: float) -> Bounds1D:  # free interval for epsilon on spot x
        ret_bounds = Bounds1D.nan()
        hyperbola = self.hyperbola_vertical(x)
        if not hyperbola.is_nan():
            ys = hyperbola.fy(epsilon)
            if len(ys) == 1:
                ret_bounds = Bounds1D(ys[0], ys[0])
            elif len(ys) == 2:
                ret_bounds = Bounds1D(min(ys), max(ys))
        return self.bounds_ver.cut(ret_bounds)

    def sample_l(self, n_l: int, n_p: int, rel_bounds: Bounds_2D = ((0, 1), (0, 1))) -> {}:
        # sample cell for n_l: # of ls, n_p: # of points per ellipsis, rel_bounds: relative xy bounds
        ls = []
        for i in range(n_l):
            l = self.bounds_l[0] + (float(i) / (n_l - 1)) * (self.bounds_l[1] - self.bounds_l[0])
            ls.append(l)

        return self.sample(ls, n_p, rel_bounds)

    def sample(self, ls: [float], n: int) -> {}:
        # sample cell with ls: array of ls to sample, n: # of points per ellipsis, rel_bounds: relative xy bounds

        bounds = self.bounds_xy

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
                        sample["ellipses"].append((l, [p]))
                    continue

                ellipsis_to_sample = self.norm_ellipsis * l

                # sample ellipsis only in bounds of cell
                ts = ellipsis_to_sample.cuts_bounds_t(bounds)

                if len(ts) == 0:
                    ts.append(0)

                for i in range(len(ts)):
                    t1 = ts[i - 1]
                    t2 = ts[i]
                    if t1 >= t2:
                        t2 += 2 * math.pi

                    d_t = t2 - t1
                    mid_t = t1 + 0.5 * d_t

                    if ellipsis_to_sample.p(mid_t).in_bounds(bounds):
                        ellipsis_sample = [ellipsis_to_sample.p(t1)]

                        np1 = math.ceil((t1 / (2 * math.pi)) * n)
                        np2 = math.ceil((t2 / (2 * math.pi)) * n)

                        for ip in range(np1, np2):
                            t = (ip / n) * 2 * math.pi
                            ellipsis_sample.append(ellipsis_to_sample.p(t))

                        ellipsis_sample.append(ellipsis_to_sample.p(t2))

                        sample["ellipses"].append((l, ellipsis_sample))

            else:  # case 2: lines are parallel
                ellipsis_to_sample = self.norm_ellipsis * l
                for ellipsis_sample in ellipsis_to_sample.cuts_bounds_p(bounds):
                    sample["ellipses"].append((l, ellipsis_sample))

        # Sample Ellipsis axis
        if not self.parallel:
            sample["axis"].append(("c", LineSegment(self.norm_ellipsis.m, self.norm_ellipsis.m +
                                                    self.norm_ellipsis.a).cuts_bounds(bounds)))
            sample["axis"].append(("d", LineSegment(self.norm_ellipsis.m, self.norm_ellipsis.m +
                                                    self.norm_ellipsis.b).cuts_bounds(bounds)))

        # Sample steepest decent lines
        sample["l-lines"].append(("l", self.l_ver.cuts_bounds(bounds)))
        sample["l-lines"].append(("l'", self.l_hor.cuts_bounds(bounds)))

        return sample


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
        self.bounds_l = Bounds1D(min(a.d_ls_point(b.p1), a.d_ls_point(b.p2), b.d_ls_point(a.p1), b.d_ls_point(a.p2)),
                         max(a.p1.d(b.p1), a.p1.d(b.p2), a.p2.d(b.p1), a.p2.d(b.p2)))

        self.parallel = about_equal(a.m, b.m)
        if not self.parallel:  # case 1: lines are not parallel
            self.s = a.intersection_p(b)  # intersection point S

            self.a = OneLineSegment(a, self.s)  # line segment a
            self.b = OneLineSegment(b, self.s)  # line segment b

            self.intersect = self.a.intersects and self.b.intersects  # does the intersection point lie on A and B
            if self.intersect:  # if line segments intersect, set min length to 0
                self.bounds_l = Bounds1D(0.0, self.bounds_l[1])

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

            norm_ellipsis = Ellipse(m + offset, c, d)

        else:  # case 2: lines are parallel
            if self.dirB:
                a = Vector(1, 1)
            else:
                a = Vector(-1, 1)

            norm_ellipsis = EllipseInfinite(self.anchor + offset, a, self.dist)

        # set cell bounds: length of a and b
        bounds_xy = Bounds1D(self.a.d.l, self.b.d.l)

        return Cell(self.parallel, self.a, self.b, norm_ellipsis, bounds_xy, self.bounds_l, offset=offset)


class Traversal:
    def __init__(self, cell_matrix: 'CellMatrix', a_cm: CM_Point, b_cm: CM_Point, points: [Vector], epsilon: float,
                 epsilons: [float]):
        self.cell_matrix = cell_matrix

        self.a_cm = a_cm
        self.a = a_cm[0]
        self.cell_a = a_cm[1]
        self.b_cm = b_cm
        self.b = b_cm[0]
        self.cell_b = b_cm[1]

        self.points = points
        self.epsilon = epsilon
        self.epsilons = epsilons

    def __str__(self):
        return "    " + str(self.epsilon) + "-Traversal:" + '\n' + \
               "      A: " + str(self.a) + " -> B: " + str(self.b) + '\n' + \
               "      Cell_A: " + str(self.cell_a) + " -> Cell_B: " + str(self.cell_b) + '\n' + \
               "      Points: " + str([str(point) for point in self.points]) + '\n' + \
               "      Epsilon: " + str(self.epsilon) + '\n' + \
               "      Epsilons: " + str(self.epsilons) + '\n' + \
               "      Decision: " + str(self.cell_matrix.decide_critical_traversal(self.cell_matrix.a_cm, self,
                                                                                   self.cell_matrix.b_cm))  # DEBUG

    @staticmethod
    def nan() -> 'Traversal':
        return Traversal(None, (Vector.nan(), (math.nan, math.nan)), (Vector.nan(), (math.nan, math.nan)),
                         [Vector.nan()], math.inf, [math.inf])

    def is_nan(self) -> bool:
        return math.isinf(self.epsilon)

    def __add__(self, other: 'Traversal') -> 'Traversal':

        if self.is_nan():
            return other
        if other.is_nan():
            return self

        assert self.b_cm == other.a_cm, "Cannot combine Traversals (B1 != A2): \n" + str(self) + '\n' + str(other)

        return Traversal(self.cell_matrix, self.a_cm, other.b_cm, self.points[:-1] + other.points,
                         max(self.epsilon, other.epsilon), self.epsilons[:-1] + other.epsilons)


class CriticalEvents:
    def __init__(self, dictionary: {} = None):
        if dictionary is None:
            dictionary = {}
        self.dictionary = dictionary

    def __str__(self):
        desc = ""
        for traversal in self.list():
            desc += '\n' + str(traversal)
        return desc

    def __getitem__(self, item) -> [Traversal]:
        if item not in self.dictionary:
            return []
        return self.dictionary[item]

    def __len__(self):
        return len(self.dictionary)

    def append(self, traversal: Traversal):
        epsilon = traversal.epsilon
        if epsilon not in self.dictionary:
            self.dictionary[epsilon] = []
        self.dictionary[epsilon].append(traversal)

    def list(self) -> [Traversal]:
        sorted_events = []
        for epsilon in self.epsilons():
            sorted_events += self[epsilon]
        return sorted_events

    def epsilons(self) -> [float]:
        return sorted(self.dictionary.keys())

    def in_bounds(self, a1: Vector, b2: Vector) -> 'CriticalEvents':
        critical_events_cut = CriticalEvents()

        for epsilon in self.epsilons():
            for traversal in self[epsilon]:
                b1 = traversal.a
                a2 = traversal.b
                if a1 < b1 and a2 < b2 and not (a1 == b1 or a2 == b2):
                    critical_events_cut.append(traversal)

        return critical_events_cut

    def in_and_on_bounds(self, a1: Vector, b2: Vector) -> 'CriticalEvents':
        critical_events_cut = CriticalEvents()

        for epsilon in self.epsilons():
            for traversal in self[epsilon]:
                b1 = traversal.a
                a2 = traversal.b
                if a1 < b1 and a2 < b2:
                    critical_events_cut.append(traversal)

        return critical_events_cut

    def critical(self, cell_matrix: 'CellMatrix', a_cm: CM_Point, b_cm: CM_Point) -> (float, [Traversal]):
        traversals = []

        if len(self) == 0:
            return -1, traversals

        critical_epsilon = self.critical_helper(cell_matrix, a_cm, b_cm, 0, len(self) - 1)
        for traversal in self[critical_epsilon]:
            if cell_matrix.decide_critical_traversal(a_cm, traversal, b_cm):
                traversals.append(traversal)

        if len(traversals) == 0:
            return -1, traversals

        return critical_epsilon, traversals

    def critical_helper(self, cell_matrix: 'CellMatrix', a_cm: CM_Point, b_cm: CM_Point, i_start_epsilon: int,
                        i_end_epsilon: int) -> [Traversal]:

        epsilons = self.epsilons()
        if i_start_epsilon + 1 >= i_end_epsilon:
            start_epsilon = epsilons[i_start_epsilon]
            if cell_matrix.decide_traversal(a_cm, b_cm, start_epsilon):
                return start_epsilon
            return epsilons[i_end_epsilon]

        i_mid_epsilon = i_start_epsilon + math.floor(0.5 * (i_end_epsilon - i_start_epsilon))
        mid_epsilon = epsilons[i_mid_epsilon]

        decision = cell_matrix.decide_traversal(a_cm, b_cm, mid_epsilon)
        if decision:
            return self.critical_helper(cell_matrix, a_cm, b_cm, i_start_epsilon, i_mid_epsilon)
        else:
            return self.critical_helper(cell_matrix, a_cm, b_cm, i_mid_epsilon, i_end_epsilon)


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
        desc += "     => Minima: " + str(self.minima())
        desc += "     => Maxima: " + str(self.maxima())
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
                orientation_left = left_hyperbola.orientation(x)
                orientation_right = right_hyperbola.orientation(x)
                if orientation_left >= 0 > orientation_right or orientation_left > 0 >= orientation_right:
                    self._maxima.append(x)
        return self._maxima

    def is_maxima(self, x: float) -> bool:  # checks if local maxima at set x
        if len(self._maxima) == 0:
            self._maxima = self.maxima()
        return x in self._maxima


class CellMatrix:
    def __init__(self, points_p: [Vector], points_q: [Vector], traverse: int = 0):
        # paths
        self.p = Path(points_p)
        self.q = Path(points_q)

        # Cell Matrix Points
        self.a_cm = (Vector(0, 0), (0, 0))
        self.b_cm = (Vector(self.p.length, self.q.length), (self.p.count - 1, self.q.count - 1))

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
                self.cells[i_p].append(cell)

        # critical events
        self.critical_events = self.calculate_critical_events()

        # DEBUG: sample hyperbolas:
        fig_both = plt.figure(figsize=plt.figaspect(0.5))
        ax_hor = fig_both.add_subplot(2, 1, 1, aspect=1, ylim=self.bounds_l, xlabel="p", ylabel="ε")
        ax_ver = fig_both.add_subplot(2, 1, 2, aspect=1, ylim=self.bounds_l, xlabel="q", ylabel="ε")
        # horizontal
        # fig_hor = plt.figure(figsize=plt.figaspect(0.5))
        for i_q in range(self.q.count + 1):
            # ax = fig_hor.add_subplot(self.q.count + 1, 1, self.q.count - i_q + 1, aspect=1, ylim=self.bounds_l, xlabel="p", ylabel="ε")
            all_points = []
            for i_p in range(self.p.count):
                points = []
                bounds = (self.p.offsets[i_p], self.p.offsets[i_p + 1])
                n_points = math.ceil(100 * (self.p.lengths[i_p] / self.p.length))
                sample = self.cross_sections_hor[i_q][i_p].sample(bounds, n_points)
                points += sample
                all_points += sample
                x, y = vectors_to_xy(points)
                # ax.plot(x, y)
            ax_hor.plot(*vectors_to_xy(all_points), label="q = " + str(i_q))
        # vertical
        # fig_ver = plt.figure(figsize=plt.figaspect(0.5))
        for i_p in range(self.p.count + 1):
            # ax = fig_ver.add_subplot(1, self.p.count + 1, i_p + 1, aspect=1, xlim=self.bounds_l, xlabel="ε", ylabel="q")
            all_points = []
            for i_q in range(self.q.count):
                points = []
                bounds = (self.q.offsets[i_q], self.q.offsets[i_q + 1])
                n_points = math.ceil(100 * (self.q.lengths[i_q] / self.q.length))
                sample = self.cross_sections_ver[i_p][i_q].sample(bounds, n_points)
                points += sample
                all_points += sample
                x, y = vectors_to_xy(points)
                # ax.plot(y, x)
            ax_ver.plot(*vectors_to_xy(all_points), label="p = " + str(i_p))
        ax_hor.legend()
        ax_ver.legend()

        # traverse
        self.traverse = traverse
        if traverse > 0:
            self.max_epsilon, self.traversals = self.do_traverse()

    def __str__(self):
        desc = "Input (" + str(self.p.count) + "x" + str(self.q.count) + ") (" + \
               str(self.p.length) + "x" + str(self.q.length) + "):\n"
        desc += " Points P:\n"
        desc += "  p_xs = " + str([point.x for point in self.p.points]) + '\n'
        desc += "  p_ys = " + str([point.y for point in self.p.points]) + '\n'
        desc += " Points Q:\n"
        desc += "  q_xs = " + str([point.x for point in self.q.points]) + '\n'
        desc += "  q_ys = " + str([point.y for point in self.q.points]) + '\n'
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
        if len(self.critical_events) > 0:
            desc += " Critical Events: " + str(self.critical_events) + '\n'
        else:
            desc += " Critical Events: None\n"
        if self.traverse > 0:
            desc += '\n'
            desc += " Global Epsilon: " + str(self.max_epsilon) + '\n'
            desc += " Traversals:\n"
            for traversal in self.traversals:
                desc += str(traversal) + '\n'

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

        # events of type b
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
                other_cross_section = other_cross_sections[path.i_rl_point(x)]
                y = other_path.offsets[i_cross_section]
                if other_cross_section.is_maxima(y):
                    critical_points.append([Vector(x, y)])

        # events of type c
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

    def calculate_critical_events(self) -> CriticalEvents:
        critical_events = CriticalEvents()

        # horizontal
        critical_points_hor = self.calculate_critical_points(self.cross_sections_ver, self.cross_sections_hor)
        for points in critical_points_hor:
            points = [point.x_to_y() for point in points]
            traversal = self.traversal_from_points(points)
            critical_events.append(traversal)

        # vertical
        critical_points_ver = self.calculate_critical_points(self.cross_sections_hor, self.cross_sections_ver)
        for points in critical_points_ver:
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

        if a == b:
            return True

        if not a < b:
            print("Error: Traversal is impossible, because not a < b: a=" + str(a) + " and b=" + str(b))
            return False

        cells = self.cells.copy()
        cells.append(cells[-1])
        cells = [cells_col + [cells_col[-1]] for cells_col in cells]

        start_i_p = cell_a[0]
        end_i_p = cell_b[0]
        start_i_q = cell_a[1]
        end_i_q = cell_b[1]

        while end_i_p > 0 and b.x <= self.p.offsets[end_i_p]:
            end_i_p -= 1
        while start_i_p < end_i_p and a.x >= self.p.offsets[start_i_p + 1]:
            start_i_p += 1
        while end_i_q > 0 and b.y <= self.q.offsets[end_i_q]:
            end_i_q -= 1
        while start_i_q < end_i_q and a.y >= self.q.offsets[start_i_q + 1]:
            start_i_q += 1

        d_p = end_i_p - start_i_p + 1
        d_q = end_i_q - start_i_q + 1
        if about_equal(a.x, b.x):
            d_p = 0
        if about_equal(a.y, b.y):
            d_q = 0

        offsets_hor = self.p.offsets[start_i_p: end_i_p + 2]
        offsets_hor[0] = a.x
        offsets_hor[-1] = b.x
        offsets_ver = self.q.offsets[start_i_q: end_i_q + 2]
        offsets_ver[0] = a.y
        offsets_ver[-1] = b.y

        bounds_hor = [Bounds1D(offsets_hor[i_p], offsets_hor[i_p + 1]) for i_p in range(d_p)]
        bounds_ver = [Bounds1D(offsets_ver[i_q], offsets_ver[i_q + 1]) for i_q in range(d_q)]

        reachable_hor = [[Bounds1D.nan() for i_q in range(d_q + 1)] for i_p in range(d_p)]
        reachable_ver = [[Bounds1D.nan() for i_q in range(d_q)] for i_p in range(d_p + 1)]

        start_cell = cells[start_i_p][start_i_q]

        if d_p > 0:
            reachable_bottom = bounds_hor[0].cut(start_cell.free_bounds_horizontal(offsets_ver[0], epsilon))
            if offsets_hor[0] in reachable_bottom:
                reachable_hor[0][0] = reachable_bottom
            for i_p in range(1, d_p):
                cell = cells[start_i_p + i_p][start_i_q]
                reachable_bottom = bounds_hor[i_p].cut(cell.free_bounds_horizontal(offsets_ver[0], epsilon))
                if reachable_bottom.start in reachable_hor[i_p - 1][0]:
                    reachable_hor[i_p][0] = reachable_bottom
            if d_q == 0:
                return b.x in reachable_hor[d_p - 1][0]

        if d_q > 0:
            reachable_left = bounds_ver[0].cut(start_cell.free_bounds_vertical(offsets_hor[0], epsilon))
            if offsets_ver[0] in reachable_left:
                reachable_ver[0][0] = reachable_left
            for i_q in range(1, d_q):
                cell = cells[start_i_p][start_i_q + i_q]
                reachable_left = bounds_ver[i_q].cut(cell.free_bounds_vertical(offsets_hor[0], epsilon))
                if reachable_left.start in reachable_ver[0][i_q - 1]:
                    reachable_ver[0][i_q] = reachable_left
            if d_p == 0:
                return b.y in reachable_ver[0][d_q - 1]

        for i_p in range(d_p):
            for i_q in range(d_q):
                cell = cells[start_i_p + i_p][start_i_q + i_q]

                reachable_left = reachable_ver[i_p][i_q]
                reachable_bottom = reachable_hor[i_p][i_q]

                free_top = bounds_hor[i_p].cut(cell.free_bounds_horizontal(offsets_ver[i_q + 1], epsilon))
                reachable_top = Bounds1D.nan()
                if not reachable_left.is_nan():
                    reachable_top = free_top
                elif not reachable_bottom.cut(free_top).is_nan():
                    reachable_top = bounds_hor[i_p].cut(Bounds1D(reachable_bottom.start, free_top.end))
                reachable_hor[i_p][i_q + 1] = reachable_top

                free_right = bounds_ver[i_q].cut(cell.free_bounds_vertical(offsets_hor[i_p + 1], epsilon))
                reachable_right = Bounds1D.nan()
                if not reachable_bottom.is_nan():
                    reachable_right = free_right
                elif not reachable_left.cut(free_right).is_nan():
                    reachable_right = bounds_ver[i_q].cut(Bounds1D(reachable_left.start, free_right.end))
                reachable_ver[i_p + 1][i_q] = reachable_right

        '''if epsilon == 1.4142135623730951:  # DEBUG
            print("A: " + str([str(p) for p in a_cm]))
            print("B: " + str([str(p) for p in b_cm]))
            print("start_i_p: " + str(start_i_p) + " bis end_i_p: " + str(end_i_p))
            print("start_i_q: " + str(start_i_q) + " bis end_i_q: " + str(end_i_q))
            print("d_p: " + str(d_p))
            print("d_q: " + str(d_q))
            print("reachable_hor: " + str([[str(b) for b in bs] for bs in reachable_hor]))
            print("reachable_ver: " + str([[str(b) for b in bs] for bs in reachable_ver]))
            print("offsets_hor: " + str([str(offset) for offset in offsets_hor]))
            print("offsets_ver: " + str([str(offset) for offset in offsets_ver]))
            print("bounds_hor: " + str([str(bound) for bound in bounds_hor]))
            print("bounds_ver: " + str([str(bound) for bound in bounds_ver]))
            print("==============")'''

        return (b.x in reachable_hor[d_p - 1][d_q]) or (b.y in reachable_ver[d_p][d_q - 1])

    def cm_point_a(self, a: Vector) -> CM_Point:
        cell_a = (self.p.i_rl_path(a.x), self.q.i_rl_path(a.y))
        return a, cell_a

    def cm_point_b(self, b: Vector) -> CM_Point:
        cell_b = (self.p.i_rl_path(b.x), self.q.i_rl_path(b.y))
        if about_equal(self.p.offsets[cell_b[0]], b.x) and 0 < cell_b[0]:
            cell_b = (cell_b[0] - 1, cell_b[1])
        if about_equal(self.q.offsets[cell_b[1]], b.y) and 0 < cell_b[1]:
            cell_b = (cell_b[0], cell_b[1] - 1)
        return b, cell_b

    def traversal_from_points(self, points: [Vector]) -> Traversal:
        start = points[0]
        end = points[-1]

        start_cm = self.cm_point_b(start)
        end_cm = self.cm_point_a(end)

        epsilons = []
        for point in points:
            tmp_epsilon = self.epsilon_from_point(point)
            epsilons.append(tmp_epsilon)
        epsilon = max(epsilons)

        return Traversal(self, start_cm, end_cm, points, epsilon, epsilons)

    def epsilon_from_point(self, point: Vector) -> float:
        epsilon = self.p.p_rl(point.x).d(self.q.p_rl(point.y))
        if about_equal(epsilon, 0):
            return 0.0
        return epsilon

    def epsilon_from_cm_point(self, point_cm: CM_Point) -> float:
        cell = point_cm[1]
        point = point_cm[0] - Vector(self.p.offsets[cell[0]], self.q.offsets[cell[1]])
        epsilon = self.p[cell[0]].frl(point.x).d(self.q[cell[1]].frl(point.y))
        if about_equal(epsilon, 0):
            return 0.0
        return epsilon

    def do_traverse(self) -> (float, [Traversal]):
        traversals = self.traverse_recursive(self.a_cm, self.critical_events, self.b_cm)
        assert len(traversals) > 0, "Error: No Traversal was found !!? Traversals: " + str(traversals)
        return traversals[0].epsilon, traversals

    def traverse_recursive(self, a_cm: CM_Point, critical_events: CriticalEvents, b_cm: CM_Point) -> [Traversal]:

        a = a_cm[0]
        b = b_cm[0]

        cc_a = a_cm[1]
        cc_b = b_cm[1]

        a_epsilon = self.epsilon_from_cm_point(a_cm)
        b_epsilon = self.epsilon_from_cm_point(b_cm)
        max_ab_epsilon = max(a_epsilon, b_epsilon)

        if a == b:
            return [Traversal(self, a_cm, b_cm, [a], max_ab_epsilon, [a_epsilon])]

        critical_event = critical_events.critical(self, a_cm, b_cm)
        critical_epsilon = critical_event[0]

        '''print("========")  # DEBUG
        print("Critical Epsilon: " + str(critical_epsilon))
        print("A: " + str(a1) + " Cell: " + str(a1_cm[1]) + " a_epsilon: " + str(a_epsilon))
        print("B: " + str(b2) + " Cell: " + str(b2_cm[1]) + " b_epsilon: " + str(b_epsilon))
        print("=")'''

        if not about_equal(max_ab_epsilon, critical_epsilon) and \
                (max_ab_epsilon > critical_epsilon or (max_ab_epsilon < critical_epsilon
                                                       and self.decide_traversal(a_cm, b_cm, max_ab_epsilon))):
            if about_equal(a.x, b.x) or about_equal(a.y, b.y) or a.x > b.x or a.y > b.y:
                print("********** Quad!!?")  # DEBUG
                return [Traversal(self, a_cm, b_cm, [a, b], max_ab_epsilon,
                                  [a_epsilon, b_epsilon])]  # Todo:berichtigen!

            a_cell = self.cells[cc_a[0]][cc_a[1]]
            b_cell = self.cells[cc_b[0]][cc_b[1]]

            if about_equal(a_epsilon, b_epsilon):
                # do steepest decent from A and B
                a2, a2_epsilon, (a_hyperbola_hor, a_hyperbola_ver) = a_cell.steepest_decent(a, 1)
                b2, b2_epsilon, (b_hyperbola_hor, b_hyperbola_ver) = b_cell.steepest_decent(b, -1)
                a2_cm = self.cm_point_a(a2)
                b2_cm = self.cm_point_b(b2)
            elif a_epsilon > b_epsilon:
                # do steepest decent from A
                a2, a2_epsilon, (a_hyperbola_hor, a_hyperbola_ver) = a_cell.steepest_decent(a, 1)
                b2, b2_epsilon, (b_hyperbola_hor, b_hyperbola_ver) = b, b_epsilon, (Hyperbola.nan(), Hyperbola.nan())
                a2_cm = self.cm_point_a(a2)
                b2_cm = b_cm
            else:
                # do steepest decent from B
                a2, a2_epsilon, (a_hyperbola_hor, a_hyperbola_ver) = a, a_epsilon, (Hyperbola.nan(), Hyperbola.nan())
                b2, b2_epsilon, (b_hyperbola_hor, b_hyperbola_ver) = b_cell.steepest_decent(b, -1)
                a2_cm = a_cm
                b2_cm = self.cm_point_b(b2)

            # set bounds
            a_bounds_hor = Bounds1D(min(a.x, a2.x), max(a.x, a2.x))
            a_bounds_ver = Bounds1D(min(a.y, a2.y), max(a.y, a2.y))
            b_bounds_hor = Bounds1D(min(b.x, b2.x), max(b.x, b2.x))
            b_bounds_ver = Bounds1D(min(b.y, b2.y), max(b.y, b2.y))

            # if possible decent to equal height
            if not about_equal(a2_epsilon, b2_epsilon) and b2_epsilon > a2_epsilon:
                a2_x = a.x
                a2_y = a.y
                if not a_hyperbola_hor.is_nan():
                    xs = a_bounds_hor.in_bounds(a_hyperbola_hor.fy(b2_epsilon))
                    if len(xs) > 0:
                        a2_x = xs[0]
                if not a_hyperbola_ver.is_nan():
                    ys = a_bounds_ver.in_bounds(a_hyperbola_ver.fy(b2_epsilon))
                    if len(ys) > 0:
                        a2_y = ys[0]
                a2 = Vector(a2_x, a2_y)
                a2_cm = self.cm_point_a(a2)
                a2_epsilon = b2_epsilon
            elif not about_equal(a2_epsilon, b2_epsilon) and a2_epsilon > b2_epsilon:
                b2_x = b.x
                b2_y = b.y
                if not b_hyperbola_hor.is_nan():
                    xs = b_bounds_hor.in_bounds(b_hyperbola_hor.fy(a2_epsilon))
                    if len(xs) > 0:
                        b2_x = xs[0]
                if not b_hyperbola_ver.is_nan():
                    ys = b_bounds_ver.in_bounds(b_hyperbola_ver.fy(a2_epsilon))
                    if len(ys) > 0:
                        b2_y = ys[0]
                b2 = Vector(b2_x, b2_y)
                b2_cm = self.cm_point_b(b2)
                b2_epsilon = a2_epsilon


            # create traversals
            traversal_a = Traversal.nan()
            traversal_b = Traversal.nan()
            if a != a2:
                traversal_a = Traversal(self, a_cm, a2_cm, [a, a2], max(a_epsilon, a2_epsilon), [a_epsilon, a2_epsilon])
            if b != b2:
                traversal_b = Traversal(self, b2_cm, b_cm, [b2, b], max(b_epsilon, b2_epsilon), [b2_epsilon, b_epsilon])

            if a == a2 and b == b2:  # DEBUG: just connect if no steepest decent possible
                return [Traversal(self, a_cm, b_cm, [a, b], max(a_epsilon, b_epsilon), [a_epsilon, b_epsilon])]

            new_critical_events = critical_events.in_and_on_bounds(a2, b2)
            rec_traversals = self.traverse_recursive(a2_cm, new_critical_events, b2_cm)
            traversals = []
            for rec_traversal in rec_traversals:
                traversals.append(traversal_a + rec_traversal + traversal_b)
            return traversals

        critical_traversals = critical_event[1]

        '''print("Critical Traversals:")  # DEBUG
        print("critical_events_count: " + str(len(critical_events)))  # DEBUG
        print("critical_traversals_count: " + str(len(critical_traversals)))  # DEBUG
        print("critical_traversals: " + str([str(ct) for ct in critical_traversals]))  # DEBUG
        print("critical_traversals: " + str([str(ct) for ct in critical_traversals]))  # DEBUG
        print("===")'''

        traversals = []
        for critical_traversal in critical_traversals:

            b1_cm = critical_traversal.a_cm
            a1_cm = critical_traversal.b_cm
            b1 = b1_cm[0]
            a1 = a1_cm[0]

            critical_events_1 = critical_events.in_bounds(a, b1)
            critical_events_2 = critical_events.in_bounds(a1, b)

            traversals_1 = self.traverse_recursive(a_cm, critical_events_1, b1_cm)
            traversals_2 = self.traverse_recursive(a1_cm, critical_events_2, b_cm)

            for traversal_1 in traversals_1:
                for traversal_2 in traversals_2:
                    traversal = traversal_1 + critical_traversal + traversal_2
                    traversals.append(traversal)

        return traversals

    def sample_l(self, n_l: int, n_p: int, heatmap: int = 100, traversals_n: int = 10) -> {}:
        # sample cell-matrix with n_l: number of ls and n_p: points per ellipses
        ls = []

        if n_l > 0:
            for i in range(n_l + 1):
                l = self.bounds_l[0] + (float(i) / n_l) * (self.bounds_l[1] - self.bounds_l[0])
                ls.append(l)
        elif n_l == -1:
            ls = [self.max_epsilon]

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
        samples["critical-traversals"] = self.critical_events.list()

        # sample traversals
        if self.traverse > 0:
            samples["traversals"] = self.traversals

        # sample a traversal over the input paths
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

    def points_for_traversal_point(self, traversal_p: Vector, c_a: int = 0, c_b: int = 0) -> [Vector, Vector]:
        # sample lines for sampling over the input paths
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

    def sample_traversal(self, traversal: Traversal, n: int) -> {}:
        # sample a specific traversal for lines between paths and 3d-plot with n points

        sample = {"in-traversal-l": [], "in-traversal": [], "traversal-3d": [], "traversal-3d-l": []}

        epsilon = traversal.epsilon
        epsilons = traversal.epsilons
        traversal_segments = []
        traversal_length = 0

        for i in range(1, len(traversal.points)):
            p1 = traversal.points[i - 1]
            p2 = traversal.points[i]
            if p1 != p2:
                linesegment = LineSegment(p1, p2)
                traversal_segments.append(linesegment)
                traversal_length += linesegment.l

        x, y, z = [], [], []  # arrays for traversal 3d-plot
        x_l, y_l, z_l = [], [], []  # arrays for traversal 3d-plot epsilon

        c_a, c_b = 0, 0  # counters for active cell

        for c_t in range(len(traversal_segments)):

            t_ls = traversal_segments[c_t]
            n_t = math.ceil((t_ls.l / traversal_length) * n)
            i_t = 0

            p1 = t_ls.p1
            while p1.x > self.p.offsets[c_a + 1] and c_a < self.p.count - 1:
                c_a += 1
            while p1.y > self.q.offsets[c_b + 1] and c_b < self.q.count - 1:
                c_b += 1

            if epsilons[c_t] >= epsilon:
                line = self.points_for_traversal_point(p1, c_a=c_a, c_b=c_b)
                sample["in-traversal-l"].append(line)

                x_l.append(p1.x)
                y_l.append(p1.y)
                z_l.append(epsilons[c_t])

                x.append(p1.x)
                y.append(p1.y)
                z.append(epsilons[c_t])

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

        if epsilons[-1] >= epsilon:
            sample["in-traversal-l"].append([self.p.points[-1], self.q.points[-1]])
            x_l.append(self.p.length)
            y_l.append(self.q.length)
            z_l.append(epsilons[-1])
        else:
            sample["in-traversal"].append([self.p.points[-1], self.q.points[-1]])
        x.append(self.p.length)
        y.append(self.q.length)
        z.append(epsilons[-1])

        sample["traversal-3d"] = [x, y, z]
        sample["traversal-3d-l"] = [x_l, y_l, z_l]

        return sample
