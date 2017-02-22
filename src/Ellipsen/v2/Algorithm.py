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
    def __init__(self, parallel: bool, a: LineSegment, b: LineSegment, norm_ellipsis: Ellipse,
                 bounds_xy: (float, float), bounds_l: (float, float), offset: Vector = Vector(0, 0),
                 do_traverse: bool = False, start: Vector = Vector(0, 0)):
        self.parallel = parallel
        self.a = a
        self.b = b
        self.norm_ellipsis = norm_ellipsis
        self.bounds_xy = bounds_xy  # cell bounds
        self.bounds_l = bounds_l  # line length bounds
        self.offset = offset
        self.end = Vector(bounds_xy[0], bounds_xy[1])  # point at top right of cell
        self.end_l = self.lp(self.end)  # l for top right of cell
        self.cps_a = []  # critical points on a
        self.cps_b = []  # critical points on b

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
            -> {}:
        ls = []
        for i in range(nl):
            l = self.bounds_l[0] + (float(i) / (nl - 1)) * (self.bounds_l[1] - self.bounds_l[0])
            ls.append(l)

        return self.sample(ls, np, rel_bounds)

    def sample(self, ls: [float], n: int, rel_bounds: ((float, float), (float, float)) = ((0, 1), (0, 1))) \
            -> {}:
        # n: # of points per ellipsis, bounds: relative xy bounds

        bounds = ((rel_bounds[0][0] * self.bounds_xy[0], rel_bounds[0][1] * self.bounds_xy[0]),
                  (rel_bounds[1][0] * self.bounds_xy[1], rel_bounds[1][1] * self.bounds_xy[1]))

        # holds ellipses, steepest decent lines and axis in form: (name, [Vector])
        sample = {"ellipses": [], "l-lines": [], "axis": []}

        # Sample Ellipses
        for l in ls:

            if l < self.bounds_l[0] or l > self.bounds_l[1]:
                continue

            if not self.parallel:  # case 1: lines are not parallel
                if l == 0:
                    p = self.norm_ellipsis.m
                    if p.in_bounds(bounds):
                        sample["ellipses"].append((l, [p + self.offset]))
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

    def lp(self, p: Vector) -> float:  # l for given point
        return self.a.frl(p.x).d(self.b.frl(p.y))

    def traverse_right(self, traversal: (float, [float], [Vector])) -> (float, [float], [Vector]):
        l1x = self.l_ver_cut

        start = traversal[2][-1] - self.offset

        max_l = -1
        end_l = -1
        end = -self.offset

        if not math.isclose(l1x.y, self.bounds_xy[1], abs_tol=1e-13):
            if l1x.y <= start.y:
                end = Vector(self.bounds_xy[0], start.y)
                end_l = self.lp(end)
                max_l = max(traversal[0], end_l)
            elif start.y < l1x.y < self.bounds_xy[1]:
                end = l1x
                end_l = self.lp(end)  # precalculate !
                max_l = max(traversal[0], end_l)

        return max_l, traversal[1] + [end_l], traversal[2] + [end + self.offset]

    def traverse_top(self, traversal: (float, [float], [Vector])) -> (float, [float], [Vector]):
        l2x = self.l_hor_cut

        start = traversal[2][-1] - self.offset

        max_l = -1
        end_l = -1
        end = -self.offset

        if not math.isclose(l2x.x, self.bounds_xy[0], abs_tol=1e-13):
            if l2x.x <= start.x:
                end = Vector(start.x, self.bounds_xy[1])
                end_l = self.lp(end)
                max_l = max(traversal[0], end_l)
            elif start.x < l2x.x < self.bounds_xy[0]:
                end = l2x
                end_l = self.lp(end)
                max_l = max(traversal[0], end_l)

        return max_l, traversal[1] + [end_l], traversal[2] + [end + self.offset]

    def traverse_top_right(self, traversal: (float, [float], [Vector])) -> (float, [float], [Vector]):
        l1x = self.l_ver_cut
        l2x = self.l_hor_cut

        max_l = -1
        end_l = -1
        end = -self.offset

        if (l1x.y > self.bounds_xy[1] or l2x.x > self.bounds_xy[0]) or \
                math.isclose(l1x.y, self.bounds_xy[1], abs_tol=1e-13) or \
                math.isclose(l2x.x, self.bounds_xy[0], abs_tol=1e-13):
            end = Vector(self.bounds_xy[0], self.bounds_xy[1])
            end_l = self.lp(end)
            max_l = max(traversal[0], end_l)

        return max_l, traversal[1] + [end_l], traversal[2] + [end + self.offset]

    def traverse(self, traversal: (float, [float], [Vector])) -> [(float, [float], [Vector])]:
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
    def __init__(self, a, b):
        # calculate shortest and longest possible line length
        self.bounds_l = (min(a.d_ls_point(b.p1), a.d_ls_point(b.p2), b.d_ls_point(a.p1), b.d_ls_point(a.p2)),
                         max(a.p1.d(b.p1), a.p1.d(b.p2), a.p2.d(b.p1), a.p2.d(b.p2)))

        self.parallel = math.isclose(a.m, b.m, abs_tol=1e-13)
        if not self.parallel:  # case 1: lines are not parallel
            self.s = intersection(a, b)  # intersection point S

            self.a = OneLineSegment(a, self.s)  # line segment a
            self.b = OneLineSegment(b, self.s)  # line segment b

            self.intersect = self.a.intersects and self.b.intersects  # does the intersection point lie on A and B
            if self.intersect:  # if line segments intersect, set min length to 0
                self.bounds_l = (0.0, self.bounds_l[1])

        else:  # case 2: lines are parallel
            self.s = Vector(math.nan, math.nan)
            self.a = a  # line segment a
            self.b = b  # line segment b
            self.intersect = False

            self.dirB = (a.r_point(a.p1+b.d) >= 0)  # do A and B point in the same direction
            b1_p = a.project_p(b.p1)  # projection of b1 on a
            self.dist = b.p1.d(b1_p)  # distance of the two lines
            self.anchor = Vector(a.r_point(b1_p) * a.l, 0)  # point with b(0) and minimal l

    def __str__(self):
        return "    Line Segment A: " + str(self.a) + "\n    Line Segment B: " + str(self.b) + '\n' + \
               "       Parallel: " + str(self.parallel) + '\n' + \
               "       Intersect (" + str(self.intersect) + "): " + str(self.s)

    def cell(self, offset: Vector = Vector(0, 0), do_traverse: bool = False, start: Vector = Vector(0, 0)) -> Cell:
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

            # Calculate ellipsis offset:
            offset_a = Vector(self.a.l, 0) * self.a.rs
            offset_b = Vector(0, self.b.l) * self.b.rs
            m = offset_a + offset_b

            norm_ellipsis = Ellipse(m, c, d)

        else:  # case 2: lines are parallel
            if self.dirB:
                a = Vector(1, 1)
            else:
                a = Vector(-1, 1)

            norm_ellipsis = EllipseInfinite(self.anchor, a, self.dist)

        # set cell bounds: length of a and b
        bounds_xy = (self.a.d.l, self.b.d.l)

        return Cell(self.parallel, self.a, self.b, norm_ellipsis, bounds_xy, self.bounds_l, offset=offset,
                    do_traverse=do_traverse, start=start)


class CellMatrix:  # : Matrix of Cells
    def __init__(self, points_a: [Vector], points_b: [Vector]):
        self.points_a = points_a
        self.points_b = points_b
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

        # critical points
        self.cps_a, self.cps_b = self.calculate_critical_points()
        self.set_critical_points()

        # traverse
        self.lowest_l = math.inf
        self.traversals = self.traverse_best(2)

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
        desc += " Lowest_l: " + str(self.lowest_l) + '\n'
        desc += " Critical Points \n"
        desc += "   A: " + str(self.cps_a) + '\n'
        desc += "   B: " + str(self.cps_b) + '\n'
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

    def calculate_critical_points(self) -> ([float], [float]):

        cps_a = []
        for i in range(self.count_a):
            cps_a.append([self.offsets_a[i]])
        cps_b = []
        for i in range(self.count_b):
            cps_b.append([self.offsets_b[i]])

        for i_a in range(self.count_a):
            for i_b in range(self.count_b):
                a = self.path_a[i_a]
                b = self.path_b[i_b]

                m_a = a.fr(0.5)
                m_b = b.fr(0.5)
                n_a = LineSegment(m_a, m_a + a.d.rotate_90_l())
                n_b = LineSegment(m_b, m_b + b.d.rotate_90_l())

                p_a = intersection(a, n_b)
                p_b = intersection(b, n_a)

                d_a = 0.5 * (p_a.d(b.p1) + p_a.d(b.p2))
                if a.contains_point(p_a) and d_a < a.p1.d(b.p1) and d_a < a.p2.d(b.p2):
                    cps_a[i_a].append(self.offsets_a[i_a] + a.rl_point(p_a))
                d_b = 0.5 * (p_b.d(a.p1) + p_b.d(a.p2))
                if b.contains_point(p_b) and d_b < a.p1.d(b.p1) and d_b < a.p2.d(b.p2):
                    cps_b[i_b].append(self.offsets_b[i_b] + b.rl_point(p_b))

        for i in range(len(cps_a)):
            cps_a[i] = sorted(cps_a[i])
        for i in range(len(cps_b)):
            cps_b[i] = sorted(cps_b[i])

        print("cps_a = " + str(cps_a))
        print("cps_b = " + str(cps_b))

        return cps_a, cps_b

    def set_critical_points(self):
        for i_a in range(self.count_a):
            for i_b in range(self.count_b):
                cell = self.cells[i_a][i_b]
                for cp_a in self.cps_a[i_a]:
                    cp = Vector(cp_a, self.offsets_b[i_b + 1])
                    l_cp = cell.lp(cp - cell.offset)
                    cell.cps_a.append((l_cp, cp))
                for cp_b in self.cps_b[i_b]:
                    cp = Vector(self.offsets_a[i_a + 1], cp_b)
                    l_cp = cell.lp(cp - cell.offset)
                    cell.cps_b.append((l_cp, cp))

    def traverse_best(self, criteria: int = 2, delete_duplicates: bool = True, start: Vector = Vector(0, 0))\
            -> [(float, [float], [Vector])]:
        i_a = 0
        i_b = 0
        while start.x > self.offsets_a[i_a + 1]:
            i_a += 1
        while start.y > self.offsets_b[i_b + 1]:
            i_b += 1

        l_start = self.cells[i_a][i_b].lp(start)
        traversals0 = self.traverse(i_a, i_b, (l_start, [l_start], [start]))
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
                # 2. lowest sum of ls
                lowest_avg_ls = math.inf
                for traversal in traversals1:
                    avg_ls = sum(traversal[1])/len(traversal[1])
                    lowest_avg_ls = min(lowest_avg_ls, avg_ls)
                traversals2 = []
                for traversal in traversals1:
                    avg_ls = sum(traversal[1])/len(traversal[1])
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
                    if traversal1[2] == traversal2[2]:
                        unique = False
                        break
                if unique:
                    traversals_set.append(traversal1)
            return traversals_set

        print("lowest_l = " + str(self.lowest_l))
        print("lowest_l = " + str(self.lowest_l))
        print("traversals0 = " + str(traversals0))
        print("traversals1 = " + str(traversals0))
        print("traversals2 = " + str(traversals0))

        return traversals

    def traverse(self, i_a: int, i_b: int, traversal: (float, [float], [Vector])) -> [(float, [float], [Vector])]:
        if i_a >= self.count_a and i_b >= self.count_b:
            if traversal[0] < self.lowest_l:
                self.lowest_l = traversal[0]
                print("lowest_l = " + str(self.lowest_l))
            return [traversal]

        traversals = []

        if i_a >= self.count_a or i_b >= self.count_b:
            if i_a < self.count_a:
                cell = self.cells[i_a][self.count_b - 1]
                if cell.end_l <= self.lowest_l + 1e-13:
                    traversals += self.traverse(i_a + 1, i_b, (max(traversal[0], cell.end_l), traversal[1] +
                                                               [cell.end_l], traversal[2] + [cell.end + cell.offset]))
            if i_b < self.count_b:
                cell = self.cells[self.count_a - 1][i_b]
                if cell.end_l <= self.lowest_l + 1e-13:
                    traversals += self.traverse(i_a, i_b + 1, (max(traversal[0], cell.end_l), traversal[1] +
                                                               [cell.end_l], traversal[2] + [cell.end + cell.offset]))
        else:
            next_traversal = self.cells[i_a][i_b].traverse(traversal)

            if self.lowest_l + 1e-13 >= next_traversal[0][0] != -1:
                traversals += self.traverse(i_a + 1, i_b, next_traversal[0])
            if self.lowest_l + 1e-13 >= next_traversal[1][0] != -1:
                traversals += self.traverse(i_a + 1, i_b + 1, next_traversal[1])
            if self.lowest_l + 1e-13 >= next_traversal[2][0] != -1:
                traversals += self.traverse(i_a, i_b + 1, next_traversal[2])

            p = traversal[2][-1]
            cell = self.cells[i_a][i_b]

            for cp in cell.cps_a:
                if cp[1].x >= p.x and cp[0] <= self.lowest_l + 1e-13:
                    traversals += self.traverse(i_a, i_b + 1, (max(traversal[0], cp[0]), traversal[1] + [cp[0]],
                                                               traversal[2] + [cp[1]]))
            for cp in cell.cps_b:
                if cp[1].y >= p.y and cp[0] <= self.lowest_l + 1e-13:
                    traversals += self.traverse(i_a + 1, i_b, (max(traversal[0], cp[0]), traversal[1] + [cp[0]],
                                                               traversal[2] + [cp[1]]))

        return traversals

    def sample_l(self, nl: int, np: int) -> {}:  # sample cell with nl: number of ls and np: points per ellipses
        ls = []

        for i in range(nl + 1):
            l = self.bounds_l[0] + (float(i) / nl) * (self.bounds_l[1] - self.bounds_l[0])
            ls.append(l)

        return self.sample(ls, np)

    def sample(self, ls: [float], np: int) -> {}:  # sample cell for given ls and np: points per ellipses
        samples = {"borders-v": [], "borders-h": [], "cells": [], "traversals": []}

        # are all ls in bounds
        for l in ls:
            if l < self.bounds_l[0] or l > self.bounds_l[1]:
                print("l: " + str(l) + " is not in bounds_l: " + str(self.bounds_l))
                ls.remove(l)

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

        # sample traversals
        for traversal in self.traversals:
            if traversal[0] != -1:
                samples["traversals"].append(traversal)

        # include size in sample
        samples["size"] = (self.length_a, self.length_b)

        return samples

    def sample_heatmap_a(self, n_a: int) -> []:
        n_b = math.floor(n_a * (self.length_b/self.length_a))
        return self.sample_heatmap(n_a, int(n_b))

    def sample_heatmap(self, n_a: int, n_b: int) -> []:
        xs = [[]]
        ys = [[]]
        zs = [[]]

        s_a = self.length_a / n_a
        s_b = self.length_b / n_b

        x = 0
        y = 0

        i_x = 0
        i_y = 0

        c_a = 0
        c_b = 0

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

    def sample_traversal(self, traversal: (float, [float], [Vector]), n: int) -> []:

        traversal_ls = []
        traversal_length = 0
        traversal_offsets = [0]
        for i in range(1, len(traversal[2])):
            p1 = traversal[2][i-1]
            p2 = traversal[2][i]
            if p1 != p2:
                ls = LineSegment(p1, p2)
                traversal_ls.append(ls)
                traversal_length += ls.l
                traversal_offsets.append(traversal_offsets[-1] + ls.l)

        step = traversal_length / n
        i_t = 0
        c_t = 0
        c_a = 0
        c_b = 0

        lines = []

        while i_t <= traversal_length + 1e-13:

            while i_t <= traversal_offsets[c_t + 1] + 1e-13:

                ls = traversal_ls[c_t]
                vec = ls.fr((i_t - traversal_offsets[c_t]) / ls.l)
                r_a = vec.x
                r_b = vec.y

                while r_a > self.offsets_a[c_a + 1] and c_a < self.count_a - 1:
                    c_a += 1
                while r_b > self.offsets_b[c_b + 1] and c_b < self.count_b - 1:
                    c_b += 1

                cell = self.cells[c_a][c_b]
                r_a -= self.offsets_a[c_a]
                r_b -= self.offsets_b[c_b]
                pa = cell.a.fr(r_a/cell.a.l)
                pb = cell.b.fr(r_b/cell.b.l)
                lines.append([pa, pb])

                i_t += step

            c_t += 1

        return lines
