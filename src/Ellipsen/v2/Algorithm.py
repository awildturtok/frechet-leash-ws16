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


TraversalType = (float, Vector, [float], [Vector], (float, float))

class Cell:
    def __init__(self, parallel: bool, a: LineSegment, b: LineSegment, norm_ellipsis: Ellipse,
                 bounds_xy: (float, float), bounds_l: (float, float), offset: Vector = Vector(0, 0)):
        self.parallel = parallel
        # line segments
        self.a = a
        self.b = b

        self.norm_ellipsis = norm_ellipsis  # normed ellipsis (l=1)
        self.bounds_xy = bounds_xy  # cell bounds
        self.bounds_l = bounds_l  # l length bounds
        self.offset = offset  # offset in cell-matrix

        # variables for traversal
        self.end = Vector(bounds_xy[0], bounds_xy[1])  # point at top right of cell local
        self.end_global = self.end + self.offset  # point at top right of cell in cell-matrix
        self.end_l = self.lp(self.end)  # l for top right of cell

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
        # intersection points of l and l' with top-right cell sides local
        self.l_ver_cut = Vector(bounds_xy[0], self.l_ver.fx(bounds_xy[0]))
        self.l_hor_cut = Vector(self.l_hor.fy(bounds_xy[1]), bounds_xy[1])
        # ls for these intersection points
        self.l_ver_cut_l = self.lp(self.l_ver_cut)
        self.l_hor_cut_l = self.lp(self.l_hor_cut)
        # intersection points in cell-matrix
        self.l_ver_cut_global = self.l_ver_cut + self.offset
        self.l_hor_cut_global = self.l_hor_cut + self.offset

        # which steepest decent traversals are possible
        self.traverses_right = not math.isclose(self.l_ver_cut.y, self.bounds_xy[1], abs_tol=1e-13)
        self.traverses_top = not math.isclose(self.l_hor_cut.x, self.bounds_xy[0], abs_tol=1e-13)
        self.traverses_top_right = self.l_ver_cut.y > self.bounds_xy[1] or self.l_hor_cut.x > self.bounds_xy[0] \
                                   or not self.traverses_right or not self.traverses_top

    def __str__(self):
        return "    Norm-" + str(self.norm_ellipsis) + '\n' + \
               "    Offset: " + str(self.offset) + '\n' + \
               "    End: " + str(self.end) + " l: " + str(self.end_l) + '\n' + \
               "    Steepest Decent Lines:\n" + \
               "      l: " + str(self.l_ver.d) + '\n' + \
               "      l': " + str(self.l_hor.d) + '\n' + \
               "    Bounds l: " + str(self.bounds_l) + '\n' + \
               "    Bounds XY: " + str(self.bounds_xy) + '\n' + \
               "    Traverses:\n" + \
               "      right: " + str(self.traverses_right) + '\n' + \
               "      top_right: " + str(self.traverses_top_right) + '\n' + \
               "      top: " + str(self.traverses_top) + '\n'

    def sample_l(self, nl: int, np: int, rel_bounds: ((float, float), (float, float)) = ((0, 1), (0, 1))) -> {}:
        # sample cell for nl: # of ls, np: # of points per ellipsis, rel_bounds: relative xy bounds
        ls = []
        for i in range(nl):
            l = self.bounds_l[0] + (float(i) / (nl - 1)) * (self.bounds_l[1] - self.bounds_l[0])
            ls.append(l)

        return self.sample(ls, np, rel_bounds)

    def sample(self, ls: [float], n: int, rel_bounds: ((float, float), (float, float)) = ((0, 1), (0, 1))) \
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

    def lp(self, p: Vector) -> float:  # length l for given point
        return self.a.frl(p.x).d(self.b.frl(p.y))

    def traverse_right(self, traversal: TraversalType) -> \
            TraversalType:  # steepest decent traversal to right
        max_l = math.inf
        end = self.l_ver_cut_global
        end_l = self.l_ver_cut_l

        if self.traverses_right:
            start = traversal[1]
            if self.l_ver_cut_global.y <= start.y:
                end = Vector(self.end_global.x, start.y)
                end_l = self.lp(end - self.offset)
                max_l = max(traversal[0], end_l)
            elif start.y < self.l_ver_cut_global.y < self.end_global.y:

                max_l = max(traversal[0], end_l)

        return max_l, end, traversal[2] + [end_l], traversal[3] + [end], (traversal[4][0] + 1, traversal[4][1])

    def traverse_top(self, traversal: TraversalType) -> \
            TraversalType:  # steepest decent traversal to top
        max_l = math.inf
        end = self.l_hor_cut_global
        end_l = self.l_hor_cut_l

        if self.traverses_top:
            start = traversal[1]
            if self.l_hor_cut_global.x <= start.x:
                end = Vector(start.x, self.end_global.y)
                end_l = self.lp(end - self.offset)
                max_l = max(traversal[0], end_l)
            elif start.x < self.l_hor_cut_global.x < self.end_global.x:
                max_l = max(traversal[0], end_l)

        return max_l, end, traversal[2] + [end_l], traversal[3] + [end], (traversal[4][0], traversal[4][1] + 1)

    def traverse_top_right(self, traversal: TraversalType) -> \
            TraversalType:  # steepest decent traversal to top-right
        max_l = math.inf
        end = self.end_global
        end_l = self.end_l

        if self.traverses_top_right:
            max_l = max(traversal[0], end_l)

        return max_l, end, traversal[2] + [end_l], traversal[3] + [end], (traversal[4][0] + 1, traversal[4][1] + 1)

    def traverse_top_right_force(self, traversal: TraversalType) -> \
            TraversalType:  # force traversal to top-right
        end = self.end_global
        end_l = self.end_l
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

        self.parallel = math.isclose(a.m, b.m, abs_tol=1e-13)
        if not self.parallel:  # case 1: lines are not parallel
            self.s = a.intersection(b)  # intersection point S

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
        bounds_xy = (self.a.d.l, self.b.d.l)

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


class CellMatrix:
    def __init__(self, points_a: [Vector], points_b: [Vector], traverse: int = 2, fast: bool = False):
        # points of path
        self.points_a = points_a
        self.points_b = points_b
        # line segments of path
        self.path_a = [LineSegment(points_a[i], points_a[i + 1]) for i in range(len(points_a) - 1)]
        self.path_b = [LineSegment(points_b[i], points_b[i + 1]) for i in range(len(points_b) - 1)]
        # count of line segments
        self.count_a = len(self.path_a)
        self.count_b = len(self.path_b)

        # real bounds of length l over the whole cell matrix
        self.bounds_l = (math.inf, 0)
        # full length of path
        self.length_a = 0
        self.length_b = 0
        # length of line segment i
        self.lengths_a = [0] * self.count_a
        self.lengths_b = [0] * self.count_b
        # summed up length of line segments up to i
        self.offsets_a = [0] * (self.count_a + 1)
        self.offsets_b = [0] * (self.count_b + 1)

        # calculate offsets & lengths
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

        # generate and save TwoLineSegments/Cells
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

        # set fast=True to quit when global minimum is reached, local minimum is disregarded
        self.fast = fast
        self.global_minimum_reached = False
        # lower limit for smallest globally reachable l
        self.min_global_l = max(self.points_a[0].d(self.points_b[0]), self.points_a[-1].d(self.points_b[-1]))

        # critical traversals
        self.critical_traversals_horizontal, self.critical_traversals_vertical = self.calculate_critical_traversals()

        # traverse
        self.traverse = traverse
        if self.traverse > 0:
            self.lowest_l = math.inf
            self.traversals = self.traverse_best(self.traverse, delete_duplicates=True)

    def __str__(self):
        desc = "Input (" + str(self.count_a) + "x" + str(self.count_a) + ") (" + \
               str(self.length_a) + "x" + str(self.length_b) + "):\n"
        desc += " Path A (l=" + str(self.length_a) + "):\n"
        for i in range(self.count_a):
            desc += "  " + str(i) + ": " + str(self.path_a[i]) + '\n'
        desc += " Path B (l=" + str(self.length_b) + "):\n"
        for i in range(self.count_b):
            desc += "  " + str(i) + ": " + str(self.path_b[i]) + '\n'
        desc += "==>\n"
        desc += " Bounds_l: " + str(self.bounds_l) + '\n'
        if self.traverse > 0:
            desc += " Lowest_l: " + str(self.lowest_l) + '\n'
            desc += " Traversals:\n"
            for traversal in self.traversals:
                desc += "   " + str(traversal) + '\n'
        desc += " Two Line Segments & Cells:\n"
        for i_a in range(self.count_a):
            for i_b in range(self.count_b):
                desc += " Cell " + str(i_a) + "x" + str(i_b) + '\n'
                desc += str(self.twoLSs[i_a][i_b]) + '\n'
                desc += str(self.cells[i_a][i_b]) + '\n'

        return desc

    def calculate_critical_traversals(self) -> ([], []):  # calculates critical traversal paths

        # horizontal
        critical_traversals_horizontal = [[[] for b in range(self.count_b)] for a in range(self.count_a)]
        for i_b in range(self.count_b):
            b = self.path_b[i_b]
            for i1_a in range(self.count_a):
                for i2_a in range(self.count_a, i1_a, - 1):
                    p1_a = self.points_a[i1_a]
                    p2_a = self.points_a[i2_a]
                    if b.d.acute(p2_a - p1_a):
                        continue
                    b_p = b.p_for_equal_dist_to_points(p1_a, p2_a)
                    if math.isnan(b_p.x):
                        continue
                    b_d = p1_a.d(b_p)
                    b_rl = b.rl_point(b_p)
                    bb_rl = self.offsets_b[i_b] + b_rl
                    if 0 - 1e-13 <= b_rl <= b.l + 1e-13:
                        max_l = b_d
                        start_p = Vector(self.offsets_a[i1_a], bb_rl)
                        traversal_ls = [b_d]
                        traversal_ps = [start_p]
                        for i3_a in range(i1_a + 1, i2_a + 1):
                            p_a = self.points_a[i3_a]
                            p_l = p_a.d(b_p)
                            max_l = max(max_l, p_l)
                            next_p = Vector(self.offsets_a[i3_a], bb_rl)
                            traversal_ls.append(p_l)
                            traversal_ps.append(next_p)
                        if b_d + 1e-13 >= max_l and (not self.fast or max_l >= self.min_global_l):
                            critical_traversals_horizontal[i1_a][i_b].append(
                                (max_l, start_p, traversal_ls, traversal_ps, (i2_a - i1_a, 0)))

        # vertical
        critical_traversals_vertical = [[[] for b in range(self.count_b)] for a in range(self.count_a)]
        for i_a in range(self.count_a):
            a = self.path_a[i_a]
            for i1_b in range(self.count_b):
                for i2_b in range(self.count_b, i1_b, - 1):
                    p1_b = self.points_b[i1_b]
                    p2_b = self.points_b[i2_b]
                    if a.d.acute(p2_b - p1_b):
                        continue
                    a_p = a.p_for_equal_dist_to_points(p1_b, p2_b)
                    if math.isnan(a_p.x):
                        continue
                    a_d = p1_b.d(a_p)
                    a_rl = a.rl_point(a_p)
                    aa_rl = self.offsets_a[i_a] + a_rl
                    if 0 - 1e-13 <= a_rl <= a.l + 1e-13:
                        max_l = a_d
                        start_p = Vector(aa_rl, self.offsets_b[i1_b])
                        traversal_ls = [a_d]
                        traversal_ps = [start_p]
                        for i3_b in range(i1_b + 1, i2_b + 1):
                            p_b = self.points_b[i3_b]
                            p_l = p_b.d(a_p)
                            max_l = max(max_l, p_l)
                            next_p = Vector(aa_rl, self.offsets_b[i3_b])
                            traversal_ls.append(p_l)
                            traversal_ps.append(next_p)
                        if a_d + 1e-13 >= max_l and (not self.fast or max_l >= self.min_global_l):
                            critical_traversals_vertical[i_a][i1_b].append(
                                (max_l, start_p, traversal_ls, traversal_ps, (0, i2_b - i1_b)))

        return critical_traversals_horizontal, critical_traversals_vertical

    def traverse_best(self, criteria: int = 2, delete_duplicates: bool = True, start: Vector = Vector(0, 0)) \
            -> [(float, [float], [Vector])]:  # traverses cell-matrix and chooses best traversal depending on criteria
        i_a = 0
        i_b = 0
        while start.x > self.offsets_a[i_a + 1]:
            i_a += 1
        while start.y > self.offsets_b[i_b + 1]:
            i_b += 1

        print("Lower Limit for minimum global l: " + str(self.min_global_l))

        l_start = self.cells[i_a][i_b].lp(start)
        traversals0 = self.traverse_rec((l_start, start, [l_start], [start], (i_a, i_b)))
        traversals = traversals0

        # select best traversal(s):
        if criteria > 0:
            # 1. lowest l
            lowest_l = math.inf
            for traversal in traversals0:
                lowest_l = min(lowest_l, traversal[0])
            traversals1 = []
            for traversal in traversals0:
                if traversal[0] <= lowest_l or math.isclose(traversal[0], lowest_l, abs_tol=1e-13):
                    traversals1.append(traversal)
            traversals = traversals1

            if criteria > 1:
                # 2. lowest average of ls
                lowest_avg_ls = math.inf
                for traversal in traversals1:
                    avg_ls = sum(traversal[2]) / len(traversal[2])
                    lowest_avg_ls = min(lowest_avg_ls, avg_ls)
                traversals2 = []
                for traversal in traversals1:
                    avg_ls = sum(traversal[2]) / len(traversal[2])
                    if avg_ls <= lowest_avg_ls or math.isclose(avg_ls, lowest_avg_ls, abs_tol=1e-10):
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

    def traverse_rec(self, traversal: TraversalType) -> \
            [TraversalType]:  # traverse cells recursive
        # indicate which cell to traverse next
        i_a = traversal[4][0]
        i_b = traversal[4][1]

        if self.global_minimum_reached and self.fast:
            return []

        # complete recursive call on reach of cell-matrix top right
        if i_a >= self.count_a and i_b >= self.count_b:
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
        if i_a >= self.count_a or i_b >= self.count_b:
            if i_a < self.count_a:
                cell = self.cells[i_a][self.count_b - 1]
                if cell.end_l <= self.lowest_l + 1e-13:
                    next_traversal = traverse_a(cell.traverse_top_right_force(traversal))
                    traversals += self.traverse_rec(next_traversal)
            if i_b < self.count_b:
                cell = self.cells[self.count_a - 1][i_b]
                if cell.end_l <= self.lowest_l + 1e-13:
                    next_traversal = traverse_b(cell.traverse_top_right_force(traversal))
                    traversals += self.traverse_rec(next_traversal)

        else:  # if traversal is not at the top- or right-side yet

            # traverse cell by steepest decent
            sd_traversal = self.cells[i_a][i_b].traverse(traversal)
            # traverse to right side?
            if self.lowest_l + 1e-13 >= sd_traversal[0][0] and not math.isinf(sd_traversal[0][0]):
                traversals += self.traverse_rec(sd_traversal[0])
            # traverse to top-right side?
            if self.lowest_l + 1e-13 >= sd_traversal[1][0] and not math.isinf(sd_traversal[1][0]):
                traversals += self.traverse_rec(sd_traversal[1])
            # traverse to top side?
            if self.lowest_l + 1e-13 >= sd_traversal[2][0] and not math.isinf(sd_traversal[2][0]):
                traversals += self.traverse_rec(sd_traversal[2])

            # traverse by critical traversal paths
            p = traversal[1]  # last point of traversal
            # horizontal
            if i_a < self.count_a - 1:
                for critical_traversal_horizontal in self.critical_traversals_horizontal[i_a + 1][i_b]:
                    max_l = critical_traversal_horizontal[0]
                    start = critical_traversal_horizontal[1]
                    if p < start and self.lowest_l + 1e-13 >= max_l and\
                            (not self.fast or max_l >= self.min_global_l - 1e-13):
                        next_traversal = traverse_a(traverse_do(traversal, critical_traversal_horizontal))
                        traversals += self.traverse_rec(next_traversal)
            # vertical
            if i_b < self.count_b - 1:
                for critical_traversal_vertical in self.critical_traversals_vertical[i_a][i_b + 1]:
                    max_l = critical_traversal_vertical[0]
                    start = critical_traversal_vertical[1]
                    if p < start and self.lowest_l + 1e-13 >= max_l and\
                            (not self.fast or max_l >= self.min_global_l - 1e-13):
                        next_traversal = traverse_b(traverse_do(traversal, critical_traversal_vertical))
                        traversals += self.traverse_rec(next_traversal)

        return traversals

    def sample_l(self, nl: int, np: int, heatmap: int = 100, traversals_n: int = 10) -> {}:
        # sample cell-matrix with nl: number of ls and np: points per ellipses
        ls = []

        for i in range(nl + 1):
            l = self.bounds_l[0] + (float(i) / nl) * (self.bounds_l[1] - self.bounds_l[0])
            ls.append(l)

        return self.sample(ls, np, heatmap=heatmap, traversals_n=traversals_n)

    def sample(self, ls: [float], np: int, heatmap: int = 100, traversals_n: int = 10) -> {}:
        # sample cell-matrix for given ls and np: points per ellipses
        samples = {"bounds-l": [], "borders-v": [], "borders-h": [], "cells": [], "traversals": []}

        # are all ls in bounds
        for l in ls:
            if l < self.bounds_l[0] or l > self.bounds_l[1]:
                print("l: " + str(l) + " is not in bounds_l: " + str(self.bounds_l))
                ls.remove(l)

        # include bounds_l in sample
        samples["bounds-l"] = [self.bounds_l[0], self.bounds_l[1]]

        # sample input
        samples["input"] = [self.points_a, self.points_b]

        # sample cells
        for i_a in range(self.count_a):
            for i_b in range(self.count_b):
                cell = self.cells[i_a][i_b]
                samples["cells"].append((str(i_a) + "x" + str(i_b), cell.sample(ls, np)))

        # sample cell borders
        for i in range(1, self.count_a):  # vertical
            samples["borders-v"].append(("border-v: " + str(i), [Vector(self.offsets_a[i], 0),
                                                                 Vector(self.offsets_a[i], self.length_b)]))
        for i in range(1, self.count_b):  # horizontal
            samples["borders-h"].append(("border-h: " + str(i), [Vector(0, self.offsets_b[i]),
                                                                 Vector(self.length_a, self.offsets_b[i])]))

        # include size in sample
        samples["size"] = (self.length_a, self.length_b)

        # sample critical traversals
        samples["critical-traversals"] = [item for sublist in [item for sublist in self.critical_traversals_horizontal
                                                               for item in sublist] for item in sublist] + \
                                         [item for sublist in [item for sublist in self.critical_traversals_vertical
                                                               for item in sublist] for item in sublist]

        # sample traversals
        if self.traverse > 0:
            for traversal in self.traversals:
                if traversal[0] != -1:
                    samples["traversals"].append(traversal)

        # sample a traversal
        if self.traverse > 0:
            traversal = self.traversals[0]
            samples["traversal"] = self.sample_traversal(traversal, traversals_n * max(self.count_a, self.count_b))

        # sample heatmap
        if heatmap > 0:
            samples["heatmap"] = self.sample_heatmap_a(heatmap)

        return samples

    def sample_heatmap_a(self, n_a: int) -> []:  # sample heat map with squares scaled by n_a divisions on a axis
        n_b = math.floor(n_a * (self.length_b / self.length_a))
        return self.sample_heatmap(n_a, int(n_b))

    def sample_heatmap(self, n_a: int, n_b: int) -> []:  # sample heatmap by n_a rectangles on a-axis and n_b on b-axis
        # coordinates
        xs = [[]]
        ys = [[]]
        zs = [[]]
        # step size on a and b
        s_a = self.length_a / n_a
        s_b = self.length_b / n_b
        # x- & y-coordinate to iterate through
        x, y = 0, 0
        # counters of x- & and y-coordinates
        i_x, i_y = 0, 0
        # counters for active cell
        c_a, c_b = 0, 0

        # iterate through cells & a-/b-axis
        while c_b < self.count_b:
            while y <= self.offsets_b[c_b + 1] or (c_b >= self.count_b - 1 and i_y <= n_b):
                while x <= self.offsets_a[c_a + 1] or (c_a >= self.count_a - 1 and i_x <= n_a):
                    xs[-1].append(x)
                    ys[-1].append(y)
                    z = self.cells[c_a][c_b].lp(Vector(x - self.offsets_a[c_a], y - self.offsets_b[c_b]))
                    zs[-1].append(z)

                    i_x += 1
                    x += s_a
                c_a += 1

                if c_a >= self.count_a and c_b < self.count_b:
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

        r_a = traversal_p.x
        r_b = traversal_p.y

        if r_a < self.offsets_a[c_a]:
            c_a = 0
        if r_b < self.offsets_b[c_b]:
            c_b = 0

        while r_a > self.offsets_a[c_a + 1] and c_a < self.count_a - 1:
            c_a += 1
        while r_b > self.offsets_b[c_b + 1] and c_b < self.count_b - 1:
            c_b += 1

        r_a -= self.offsets_a[c_a]
        r_b -= self.offsets_b[c_b]
        a = self.path_a[c_a]
        b = self.path_b[c_b]
        pa = a.frl(r_a)
        pb = b.frl(r_b)

        return [pa, pb]

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
            while p1.x > self.offsets_a[c_a + 1] and c_a < self.count_a - 1:
                c_a += 1
            while p1.y > self.offsets_b[c_b + 1] and c_b < self.count_b - 1:
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
            sample["in-traversal-l"].append([self.points_a[-1], self.points_b[-1]])
            x_l.append(self.length_a)
            y_l.append(self.length_b)
            z_l.append(ls[-1])
        else:
            sample["in-traversal"].append([self.points_a[-1], self.points_b[-1]])
        x.append(self.length_a)
        y.append(self.length_b)
        z.append(ls[-1])

        sample["traversal-3d"] = [x, y, z]
        sample["traversal-3d-l"] = [x_l, y_l, z_l]

        return sample
