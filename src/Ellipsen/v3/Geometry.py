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


tol = 1e-13  # global absolute tolerance


def about_equal(f1: float, f2: float) -> bool:
    return math.isclose(f1, f2, abs_tol=tol)


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

    def __hash__(self):
        return hash((self.x, self.y))

    def __lt__(self, other: 'Vector') -> int:
        return self.x <= other.x and self.y <= other.y

    # Vector Arithmetic

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

    def norm(self) -> 'Vector':  # normalized vector (lenght = 1)
        return self * (1 / self.l)

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

    def frl(self, r: float) -> Vector:  # point for set parameter r = [0,l]
        return self.p1 + self.d * (r/self.l)

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
        d_l = self.d_l_point(p)
        pl = p + self.d.rotate_90_l().norm() * d_l
        pr = p + self.d.rotate_90_r().norm() * d_l
        if self.point_on(pl):
            return pl
        else:
            return pr

    def project_p_rl(self, p: Vector) -> float:  # projects point onto line and returns parameter r
        projection = self.project_p(p)
        return self.rl_point(projection)

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

    def intersection(self, b: 'LineSegment'):  # calculates the intersection point of two line segments
        if math.isinf(self.m) and not math.isinf(b.m):
            return Vector(self.n, b.fx(self.n))
        elif math.isinf(b.m) and not math.isinf(self.m):
            return Vector(b.n, self.fx(b.n))
        elif (self.m - b.m) != 0:
            x = (b.n - self.n) / (self.m - b.m)
        else:
            x = float("nan")
        return Vector(x, self.fx(x))

    def p_for_equal_dist_to_points(self, p1: Vector, p2: Vector) -> Vector:
        # point on line that has equal distance to p1 and p2
        if p1 == p2:
            return Vector(math.nan, math.nan)
        m = (p1 + p2) * 0.5
        v = (p2 - p1).rotate_90_l()
        ls = LineSegment(m, m + v)
        if not about_equal(self.m, ls.m):
            return self.intersection(ls)
        elif self.point_on(m):
            return m
        else:
            return Vector(math.nan, math.nan)

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

    def parabola_with_point(self, point: Vector) -> 'Parabola':
        p1 = Vector(0, self.p1.d(point))
        p2 = Vector(0.5 * self.l, self.fr(0.5).d(point))
        p3 = Vector(self.l, self.p2.d(point))
        return Parabola(p1, p2, p3)

    def parabola_with_line(self, line: 'LineSegment') -> 'Parabola':
        l = math.sqrt(math.pow(self.l, 2) + math.pow(line.l, 2))
        p1 = Vector(0, self.p1.d(line.p1))
        p2 = Vector(0.5 * l, self.fr(0.5).d(line.fr(0.5)))
        p3 = Vector(l, self.p2.d(line.p2))
        return Parabola(p1, p2, p3)


class Parabola:
    def __init__(self, p1: Vector, p2: Vector, p3: Vector):
        self.ps = [p1, p2, p3]  # points defining the parabola
        xs = np.array([[math.pow(p.x, 2), p.x, 1] for p in self.ps])
        ys = np.array([p.y for p in self.ps])
        [a, b, c] = np.linalg.solve(xs, ys)
        self.s = Vector(-0.5 * (b / a), c - 0.25 * (math.pow(b, 2) / a))
        self.a = a

    def __str__(self):
        return "Parabola: S" + str(self.s) + " func:" + str(self.func_str())

    def func_str(self) -> str:
        return str(self.a) + "(x - " + str(self.s.x) + ")^2 + " + str(self.s.y)

    def fx(self, x: float) -> float:  # y-value for set x-value
        return self.a * math.pow(x - self.s.x, 2) + self.s.y

    def fy(self, y: float) -> float:  # x-value(s) for set y-value
        if about_equal(y, self.s.y):
            return [self.s.x]
        elif (self.a > 0 and y > self.s.y) or (self.a < 0 and y < self.s.y):
            w = math.sqrt((y - self.s.y) / self.a)
            x1 = w + self.s.x
            x2 = -w + self.s.x
            return [x1, x2]
        return []

    def px(self, x: float) -> Vector:  # point for set x-value
        return Vector(x, self.fx(x))

    def py(self, y: float) -> float:  # points(s) for set y-value
        xs = self.fy(y)
        ps = []
        for x in xs:
            ps.append(self.px(x))
        return ps

    def move_x(self, d_x: float):  # moves parabola on x-axis
        d_p = Vector(d_x, 0)
        self.move_p(d_p)

    def move_y(self, d_y: float):  # moves parabola on y-axis
        d_p = Vector(0, d_y)
        self.move_p(d_p)

    def move_p(self, d_p: Vector):  # moves parabola on y-axis
        self.s += d_p

    def cuts_line(self, line: LineSegment) -> [Vector]:
        cuts = []
        # Todo
        return cuts

    def cuts_linesegment(self, linesegment: LineSegment):
        cuts_line = self.cuts_line(linesegment)
        cuts_linesegment = []
        for cut in cuts_line:
            if linesegment.contains_point(cut):
                cuts_linesegment.append(cut)
        return cuts_linesegment
    
    def cuts_parabola(self, other: 'Parabola'):
        cuts = []
        # Todo: solve self.func - other.func = 0
        return cuts

    def cuts_bounds_vert(self, bounds: (float, float)) -> [Vector]:
        return [self.px(bounds[0]), self.px(bounds[1])]

    def sample(self, bounds: (float, float), n: int) -> [Vector]:  # samples parabola in bounds with n edges
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
