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

import math
from typing import Tuple
import numpy as np
from bisect import bisect


Bounds_1D = (float, float)
Bounds_2D = (Bounds_1D, Bounds_1D)
tol = 1e-13  # global absolute tolerance


def about_equal(x1: float, x2: float) -> bool:
    return math.isclose(x1, x2, abs_tol=tol)


def in_bounds(x: float, bounds: Bounds_1D) -> bool:
    return bounds[0] <= x <= bounds[1]


class Vector:
    def __init__(self, x: float, y: float):
        self.x = float(x)
        self.y = float(y)
        self.l = math.sqrt(x ** 2 + y ** 2)

    @staticmethod
    def from_tuple(t: Tuple[float, float]) -> 'Vector':
        return Vector(t[0], t[1])

    def to_tuple(self) -> Tuple[float, float]:
        return self.x, self.y

    def __hash__(self):
        return hash((self.x, self.y))

    def __str__(self):
        return "(%s,%s)" % (self.x, self.y)

    @staticmethod
    def nan():
        return Vector(float("nan"), float("nan"))

    # Vector Arithmetic

    def __lt__(self, other: 'Vector') -> int:
        return self.x <= other.x and self.y <= other.y

    def __add__(self, other: 'Vector') -> 'Vector':
        x = self.x + other.x
        y = self.y + other.y
        return Vector(x, y)

    def __neg__(self) -> 'Vector':
        return Vector(-self.x, -self.y)

    def __sub__(self, other: 'Vector') -> 'Vector':
        x = self.x - other.x
        y = self.y - other.y
        return Vector(x, y)

    def __mul__(self, s: float) -> 'Vector':
        x = self.x * s
        y = self.y * s
        return Vector(x, y)

    def __abs__(self) -> float:
        return self.l

    def __eq__(self, other) -> bool:
        return about_equal(self.x, other.x) and about_equal(self.y, other.y)

    def __ne__(self, other) -> bool:
        return not self.__eq__(other)

    def d(self, other) -> float:  # distance to other point
        d = other - self
        return d.l

    def x_to_y(self) -> 'Vector':  # switches x and y coordinates
        return Vector(self.y, self.x)

    def norm(self) -> 'Vector':  # normalized vector (lenght = 1)
        return self * (1 / self.l)

    def scalar(self, other: 'Vector') -> bool:
        return about_equal(self.x / other.x, self.y / other.y)

    def dot_product(self, other: 'Vector') -> float:
        return self.x * other.x + self.y * other.y

    def cross_product(self, other: 'Vector') -> float:
        return self.x * other.y - self.y * other.x

    def acute(self, other) -> bool:  # returns true if angle between vectors is < 90°
        return self.dot_product(other) > 0

    def angle(self) -> float:
        return math.atan2(self.y, self.x)

    def rotate(self, a: float):  # rotate vector by angle a counterclockwise
        x = self.x * math.cos(a) - self.y * math.sin(a)
        y = self.x * math.sin(a) + self.y * math.cos(a)
        return Vector(x, y)

    def rotate_90_l(self) -> 'Vector':  # rotate 90° to the left
        return Vector(self.y, -self.x)

    def rotate_90_r(self) -> 'Vector':  # rotate 90° to the right
        return Vector(-self.y, self.x)

    def in_bounds(self, bounds: ((float, float), (float, float))) -> bool:
        return self.in_bounds_x(bounds[0]) and self.in_bounds_y(bounds[1])

    def in_bounds_x(self, bounds: (float, float)) -> bool:
        return bounds[0] - tol <= self.x <= bounds[1] + tol

    def in_bounds_y(self, bounds: (float, float)) -> bool:
        return bounds[0] - tol <= self.y <= bounds[1] + tol


class Circle:
    def __init__(self, m: Vector, r: float):
        self.m = m
        self.r = r

    def __str__(self):
        return "Circle: M" + str(self.m) + " r:" + str(self.r)

    def p_no_offset(self, t: float) -> Vector:  # point on ellipse for set parameter t without regard for offset
        return Vector(math.cos(t) * self.r, math.sin(t) * self.r)

    def p(self, t: float) -> Vector:  # point on ellipse for set parameter t
        return self.p_no_offset(t) + self.m

    def contains_point(self, p: Vector) -> bool:
        return self.m.d(p) <= self.r + tol

    def contains_points(self, ps: [Vector]) -> bool:
        for p in ps:
            if not self.contains_point(p):
                return False
        return True


class LineSegment:
    def __init__(self, p1: Vector, p2: Vector):

        assert p1 != p2, "Error: Start- and Endpoint are equal. p1: " + str(p1) + " p2: " + str(p2)

        self.p1 = p1  # start point
        self.p2 = p2  # end point
        self.d = p2 - p1  # difference vector
        self.l = self.d.l  # length

        # calculate slope m and y(or x)-intercept n of line defined by the line segment
        if (p2.x - p1.x) != 0:
            self.m = (p2.y - p1.y) / (p2.x - p1.x)
            self.n = p1.y - self.m * p1.x
        else:
            self.m = float("inf")
            self.n = p1.x

    def __str__(self):
        return str(self.p1) + "->" + str(self.p2) + ": l=" + str(self.l) + " m=" + str(self.m) + " n=" + str(self.n)

    # Line Arithmetic

    def fr(self, r: float) -> Vector:  # point for set parameter r = [0,1]
        return self.p1 + self.d * r

    def frl(self, rl: float) -> Vector:  # point for set parameter r = [0,l]
        return self.p1 + self.d * (rl/self.l)

    def fx(self, x: float) -> float:  # y-value for set x-value
        if not math.isinf(self.m):
            return self.m * x + self.n
        else:
            return float("nan")

    def fy(self, y: float) -> float:  # x-value for set y-value
        if math.isinf(self.m):
            return self.n
        elif self.m != 0:
            return (y - self.n) / self.m
        else:
            return float("nan")

    def rx(self, x: float) -> float:  # parameter r for set y-value
        if not math.isinf(self.m):
            return (x - self.p1.x) / (self.p2.x - self.p1.x)
        else:
            return float("nan")

    def ry(self, y: float) -> float:  # parameter r for set x-value
        if self.m != 0:
            return (y - self.p1.y) / (self.p2.y - self.p1.y)
        else:
            return float("nan")

    def contains_point(self, p: Vector) -> bool:  # does point p lie on the line segment
        return (self.m == float("inf") or 0 - tol <= self.rx(p.x) <= 1 + tol) and\
               (self.m == 0 or 0 - tol <= self.ry(p.y) <= 1 + tol)

    def r_point(self, p: Vector) -> float:  # parameter r: [0,1] for set point
        # result only relevant for case p lies on line
        rx = self.rx(p.x)
        ry = self.ry(p.y)
        if math.isinf(self.m):
            rx = ry
        if self.m == 0:
            ry = rx
        return 0.5 * (rx + ry)

    def rl_point(self, p: Vector) -> float:  # parameter r: [0,l] for set point
        # result only relevant for case p lies on line
        return self.r_point(p) * self.l

    def d_l_point(self, p: Vector) -> float:  # closest distance of p to line
        return abs(self.d.cross_product(p - self.p1)) / self.d.l

    def project_p(self, p: Vector) -> Vector:  # projects point onto line
        v = self.d.rotate_90_l()
        ls = LineSegment(p, p + v)
        return self.intersection_p(ls)

    def project_p_rl(self, p: Vector) -> float:  # projects point onto line and returns parameter rl
        v = self.d.rotate_90_l()
        ls = LineSegment(p, p + v)
        return self.intersection_rl(ls)

    def d_ls_point(self, p: Vector) -> float:  # closest distance of p to line segment
        d_l = self.d_l_point(p)  # closest distance of p to line
        projection_p = self.project_p(p)
        if self.contains_point(projection_p):  # if the projection of p lies on the segment return d_l
            return d_l
        else:
            return min(self.p1.d(p), self.p2.d(p))

    def point_above(self, p: Vector) -> bool:  # is given point above line
        return p.y > self.fx(p.x)

    def point_right(self, p: Vector) -> bool:  # is given point on the right side of line
        return p.x > self.fy(p.y)

    def point_on(self, p: Vector) -> bool:  # is given point on line
        return about_equal(p.y, self.fx(p.x)) or about_equal(p.x, self.fy(p.y))

    def intersection_r(self, ls: 'LineSegment') -> float:  # calculates the parameter r of the intersection point
        if not about_equal(self.m, ls.m):
            a = np.array([[-self.d.x, ls.d.x], [-self.d.y, ls.d.y]])
            b = np.array([self.p1.x - ls.p1.x, self.p1.y - ls.p1.y])
            [self_r, ls_r] = np.linalg.solve(a, b)
            return self_r
        else:
            return float('nan')

    def intersection_p(self, ls: 'LineSegment') -> Vector:  # calculates the intersection point of two line segments
        r = self.intersection_r(ls)
        return self.fr(r)

    def intersection_rl(self, ls: 'LineSegment') -> float:  # calculates the parameter r of the intersection point
        r = self.intersection_r(ls)
        return r * self.l

    def r_for_equal_dist_to_points(self, p1: Vector, p2: Vector) -> float:
        # parameter r for point on line that has equal distance to p1 and p2
        if p1 == p2:
            p = p1
            v = self.d.rotate_90_r()
        else:
            p = (p1 + p2) * 0.5
            v = (p2 - p1).rotate_90_l()
        ls = LineSegment(p, p + v)

        return self.intersection_r(ls)

    def rl_for_equal_dist_to_points(self, p1: Vector, p2: Vector) -> float:
        return self.r_for_equal_dist_to_points(p1, p2) * self.l

    def p_for_equal_dist_to_points(self, p1: Vector, p2: Vector) -> Vector:
        return self.fr(self.r_for_equal_dist_to_points(p1, p2))

    def cuts_bounds(self, bounds: ((float, float), (float, float))) -> [Vector]: # points where line cuts given bounds
        x1 = bounds[0][0]
        x2 = bounds[0][1]
        y1 = bounds[1][0]
        y2 = bounds[1][1]

        yx1 = self.fx(x1)
        yx2 = self.fx(x2)
        xy1 = self.fy(y1)
        xy2 = self.fy(y2)

        points = []

        if x1 <= xy1 <= x2:
            points.append(Vector(xy1, y1))
        if x1 <= xy2 <= x2:
            points.append(Vector(xy2, y2))
        if y1 < yx1 < y2:
            points.append(Vector(x1, yx1))
        if y1 < yx2 < y2:
            points.append(Vector(x2, yx2))

        return points

    def hyperbola_with_point(self, point: Vector) -> 'Hyperbola':
        projection_rl = self.project_p_rl(point)
        projection_d = self.d_l_point(point)
        s = Vector(projection_rl, projection_d)
        return Hyperbola(s)

    def hyperbola_with_line(self, line: 'LineSegment') -> 'Hyperbola':  #Todo: rewrite (see hyperbola_with_point)
        l = math.sqrt(math.pow(self.l, 2) + math.pow(line.l, 2))
        p1 = Vector(0, self.p1.d(line.p1))
        p2 = Vector(0.5 * l, self.fr(0.5).d(line.fr(0.5)))
        p3 = Vector(l, self.p2.d(line.p2))
        return Hyperbola(p1, p2, p3)


class Path:
    def __init__(self, points: [Vector]):
        self.points = points
        self.count = len(points) - 1
        self.segments = []
        self.length = 0
        self.lengths = []
        self.offsets = [0]

        for i_segment in range(self.count):
            p1 = points[i_segment]
            p2 = points[i_segment + 1]
            segment = LineSegment(p1, p2)
            self.segments.append(segment)
            length = segment.l
            self.length += length
            self.lengths.append(length)
            self.offsets.append(self.length)

    def __str__(self):
        desc = " (l=" + str(self.length) + "):\n"
        for i in range(self.count):
            desc += "  " + str(i) + ": " + str(self.segments[i]) + '\n'

        return desc

    # Path Arithmetic

    def i_rl(self, rl: float) -> float:  # index for set parameter rl
        if not 0 <= rl <= self.length:
            print("Error: Parameter rl=" + str(rl) + " is not in bounds: [0, " + str(self.length) + "].")
            return math.nan

        i_segment = bisect(self.offsets[1:], rl)
        return min(i_segment, self.count - 1)

    def p_rl(self, rl: float) -> Vector:  # point for set parameter rl
        i_segment = self.i_rl(rl)

        if math.isnan(i_segment):
            return Vector.nan()
        return self.segments[i_segment].frl(rl - self.offsets[i_segment])


class Hyperbola:
    def __init__(self, s: Vector):
        self.s = s

    def __str__(self):
        return "hyperbola: S" + str(self.s) + " func: " + str(self.func_str())

    @staticmethod
    def nan() -> 'Hyperbola':
        return Hyperbola(Vector.nan())

    def is_nan(self) -> bool:
        return math.isnan(self.s.x)

    def func_str(self) -> str:
        return "sqrt( " + str(self.s.y) + "^2 + ( x - " + str(self.s.x) +  " )^2 )"

    def fx(self, x: float) -> float:  # y-value for set x-value
        return math.sqrt(math.pow(self.s.y, 2) + math.pow(x - self.s.x, 2))

    def px(self, x: float) -> Vector:  # point for set x-value
        return Vector(x, self.fx(x))

    def fy(self, y: float) -> [float]:  # x-value(s) for set y-value
        if y < self.s.y:
            return []
        elif y == self.sy:
            return [self.s.x]
        else:
            w = math.sqrt(math.pow(y, 2) - math.pow(self.s.y, 2))
            return [self.s.x - w, self.s.x + w]

    def py(self, y: float) -> [Vector]:  # point(s) for set y-value
        xs = self.fy(y)
        return [Vector(x, y) for x in xs]

    def fax(self, x: float) -> float:  # slope for set x-value
        return (x - self.x) / self.fx(x)

    def orientation(self, x: float) -> float:
        # for set x: returns positive if falling, 0 if constant, negative if rising
        return x - self.s.x

    def reflect_x(self, x: float) -> 'Hyperbola':  # reflects hyperbola over vertical at set x
        return Hyperbola(self.s + Vector(2 * self.orientation(x), 0))

    def move_x(self, d_x: float) -> 'Hyperbola':  # moves hyperbola on x-axis
        return self.move_p(Vector(d_x, 0))

    def move_y(self, d_y: float) -> 'Hyperbola':  # moves hyperbola on y-axis
        return self.move_p(Vector(0, d_y))

    def move_p(self, d_p: Vector) -> 'Hyperbola':  # moves hyperbola on y-axis
        return Hyperbola(self.s + d_p)

    def cuts_bounds_ver(self, bounds: (float, float)) -> [Vector]:
        return [self.px(bounds[0]), self.px(bounds[1])]

    def intersects_hyperbola(self, other: 'Hyperbola') -> Vector:
        s1 = self.s
        s2 = other.s
        if s1 == s2:
            x = float("inf")
        elif about_equal(s1.x, s2.x):
            x = float("nan")
        else:
            x = 0.5 * (s1.y**2 + s1.x**2 - s2.y**2 - s2.x**2) / (s1.x - s2.x)
        return self.px(x)

    def intersects_hyperbola_in_bounds(self, other: 'Hyperbola', bounds: Bounds_2D) -> Vector:
        intersection = self.intersects_hyperbola(other)
        if bounds[0] <= intersection.x <= bounds[1]:
            return intersection
        return Vector.nan()

    def intersects_hyperbola_in_bounds_critical(self, other: 'Hyperbola', bounds: Bounds_2D) -> Vector:
        intersection = self.intersects_hyperbola_in_bounds(other, bounds)
        x = intersection.x
        if self.orientation(x) <= 0 <= other.orientation(x):
            return intersection
        return Vector.nan()

    def sample(self, bounds: (float, float), n: int) -> [Vector]:  # samples hyperbola in bounds with n edges
        sample_points = []
        width = bounds[1] - bounds[0]
        for i in range(n + 1):
            x = width * (i / n) + bounds[0]
            sample_points.append(self.px(x))
        return sample_points


class Ellipse:
    def __init__(self, m: Vector, a: Vector, b: Vector):
        self.m = m  # midpoint
        self.a = a  # axis vector a
        self.b = b  # axis vector b

    def __str__(self):
        return "Ellipsis: M" + str(self.m) + " a:" + str(self.a) + " b:" + str(self.b)

    def __mul__(self, l: float):  # scale ellipsis to l
        a = self.a * l
        b = self.b * l
        return Ellipse(self.m, a, b)

    def p_no_offset(self, t: float) -> Vector:  # point on ellipse for set parameter t without regard for offset
        return self.a * math.cos(t) + self.b * math.sin(t)

    def p(self, t: float) -> Vector:  # point on ellipse for set parameter t
        return self.p_no_offset(t) + self.m

    def tx(self, x: float) -> [float]:  # parameter t for given x-value
        return self.txy(self.a.x, self.b.x, self.m.x, x)

    def ty(self, y: float) -> [float]:  # parameter t for given y-value
        return self.txy(self.a.y, self.b.y, self.m.y, y)

    @staticmethod
    def txy(a: float, b: float, m: float, xy: float) -> [float]:  # helper function for tx and ty
        if a != 0 and not about_equal(m, xy + a):
            w = a**2 + b**2 - m**2 + 2*m*xy - xy**2
            if w < 0:
                return []
            elif w == 0:
                return [2*math.atan2(b, a-m+xy)]
            return [2*math.atan2(b - math.sqrt(w), a-m+xy), 2*math.atan2(b + math.sqrt(w), a-m+xy)]
        else:
            return [math.pi, 2*(math.pi-math.atan2(a, b))]

    def cuts_bounds_t(self, bounds: ((float, float), (float, float))) -> [float]:
        # parameters t where ellipsis cuts given bounds
        tx1 = self.tx(bounds[0][0])
        tx2 = self.tx(bounds[0][1])
        ty1 = self.ty(bounds[1][0])
        ty2 = self.ty(bounds[1][1])

        tx = tx1 + tx2
        ty = ty1 + ty2

        ts = []

        for t in tx:
            t %= 2 * math.pi
            if t < 0:
                t += 2*math.pi
            p = self.p(t)
            if self.p(t).in_bounds_y(bounds[1]):
                ts.append(t)

        for t in ty:
            t %= (2*math.pi)
            if t < 0:
                t += 2*math.pi
            if self.p(t).in_bounds_x(bounds[0]):
                ts.append(t)

        if len(ts) == 0:
            return []

        ts.sort()

        ret_ts = [ts[0]]
        for i in range(len(ts) - 1):
            if not about_equal(ts[i], ts[i+1]):
                ret_ts.append(ts[i+1])

        return ret_ts

    def cuts_bounds_p(self, bounds: ((float, float), (float, float))) -> [Vector]:
        # points where ellipsis cuts given bounds
        points = []
        ts = self.cuts_bounds_t(bounds)

        for t in ts:
            points.append(self.p(t))

        return points


class EllipseInfinite:  # Ellipsis shadow where one axis is infinite (for parallel inputs)
    def __init__(self, m: Vector, a: Vector, l1: float, l2: float = 0):
        self.m = m  # anchor for the middle
        self.m1 = m  # left line anchor
        self.m2 = m  # right line anchor
        if l2 > l1:  # move left and right anchors for given l
            dm = Vector(math.sqrt(l2 ** 2 - l1 ** 2), 0)
            self.m1 -= dm
            self.m2 += dm
        self.l1 = l1
        self.l2 = l2
        self.a = a  # axis vector

    def __str__(self):
        return "EllipseInfinite: M" + str(self.m) + " a:" + str(self.a) + " l1:" + str(self.l1)

    def __mul__(self, l: float):  # scale to l
        return EllipseInfinite(self.m, self.a, self.l1, l)

    def cuts_bounds_p(self, bounds: ((float, float), (float, float))) -> [Vector]:
        # points where the lines cut given bounds
        points = []

        if self.l2 > self.l1:
            points.append(LineSegment(self.m1, self.m1 + self.a).cuts_bounds(bounds))
            points.append(LineSegment(self.m2, self.m2 + self.a).cuts_bounds(bounds))
        else:
            points.append(LineSegment(self.m, self.m + self.a).cuts_bounds(bounds))

        return points
