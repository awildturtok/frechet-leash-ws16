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
        return math.isclose(self.x, other.x, rel_tol=1e-13) and math.isclose(self.y, other.y, rel_tol=1e-13)

    def __ne__(self, other) -> bool:
        return not self.__eq__(other)

    def d(self, other) -> float:  # distance to other point
        d = other - self
        return d.l

    def norm(self) -> 'Vector':  # normalized vector (lenght = 1)
        return self * (1 / self.l)

    def cross_product(self, other) -> float:
        return self.x * other.y - self.y * other.x

    def angle(self) -> float:
        return math.atan2(self.y, self.x)

    def rotate(self, a: float):
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
        return bounds[0] <= self.x <= bounds[1]

    def in_bounds_y(self, bounds: (float, float)) -> bool:
        return bounds[0] <= self.y <= bounds[1]


class LineSegment:
    def __init__(self, p1: Vector, p2: Vector):

        assert p1 != p2, "Error: Start- and Endpoint are equal."

        self.p1 = p1  # Startpoint P1
        self.p2 = p2  # Endpoint P2
        self.d = p2 - p1  # Difference Vector
        self.l = self.d.l  # Länge der Strecke

        # calculate slope m and y(or x)-intercept n
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

    def fx(self, x: float) -> float:  # y-Value for set x-Value
        if not math.isinf(self.m):
            return self.m * x + self.n
        else:
            return float("nan")

    def fy(self, y: float) -> float:  # x-Value for set y-Value
        if math.isinf(self.m):
            return self.n
        elif self.m != 0:
            return (y - self.n) / self.m
        else:
            return float("nan")

    def rx(self, x: float) -> float:  # parameter r for set y-Value
        if not math.isinf(self.m):
            return (x - self.p1.x) / (self.p2.x - self.p1.x)
        else:
            return float("nan")

    def ry(self, y: float) -> float:  # parameter r for set x-Value
        if self.m != 0:
            return (y - self.p1.y) / (self.p2.y - self.p1.y)
        else:
            return float("nan")

    def contains_point(self, p: Vector) -> bool:  # does point p lie on the line segment
        return (self.m == float("inf") or 0 <= self.rx(p.x) <= 1) and (self.m == 0 or 0 <= self.ry(p.y) <= 1)

    def r_point(self, p: Vector) -> float:  # parameter r for set point (result only relevant for case p lies on line)
        rx = self.rx(p.x)
        ry = self.ry(p.y)
        if math.isinf(self.m):
            rx = ry
        if self.m == 0:
            ry = rx
        return 0.5 * (rx + ry)

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
        return math.isclose(p.y, self.fx(p.x), rel_tol=1e-13)

    def cuts_bounds(self, bounds: ((float, float), (float, float))) -> [Vector]:
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


class Ellipse:
    def __init__(self, m: Vector, a: Vector, b: Vector):
        self.m = m  # midpoint
        self.a = a  # axis vector a
        self.b = b  # axis vector b

    def __str__(self):
        return "Ellipsis: M" + str(self.m) + " a:" + str(self.a) + " b:" + str(self.b)

    def __mul__(self, l: float):
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
        if a != 0 and not math.isclose(m, xy + a, rel_tol=1e-13):
            w = a**2 + b**2 - m**2 + 2*m*xy - xy**2
            if w < 0:
                return []
            elif w == 0:
                return [2*math.atan2(b, a-m+xy)]
            return [2*math.atan2(b - math.sqrt(w), a-m+xy), 2*math.atan2(b + math.sqrt(w), a-m+xy)]
        else:
            return [math.pi, 2*(math.pi-math.atan2(a, b))]

    def cuts_bounds_t(self, bounds: ((float, float), (float, float))) -> [float]:
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
            if not math.isclose(ts[i], ts[i+1], rel_tol=1e-13):
                ret_ts.append(ts[i+1])

        return ret_ts

    def cuts_bounds_p(self, bounds: ((float, float), (float, float))) -> [Vector]:
        points = []
        ts = self.cuts_bounds_t(bounds)

        for t in ts:
            points.append(self.p(t))

        return points


class EllipseInfinite:
    def __init__(self, m: Vector, a: Vector, l1: float, l2: float = 0):
        self.m = m  # anchor
        self.m1 = m
        self.m2 = m
        if l2 > l1:
            dm = Vector(math.sqrt(l2 ** 2 - l1 ** 2), 0)
            self.m1 -= dm
            self.m2 += dm
        self.l1 = l1
        self.l2 = l2
        self.a = a

    def __str__(self):
        return "EllipseInfinite: M" + str(self.m) + " a:" + str(self.a) + " l1:" + str(self.l1)

    def __mul__(self, l: float):
        return EllipseInfinite(self.m, self.a, self.l1, l)

    def cuts_bounds_p(self, bounds: ((float, float), (float, float))) -> [Vector]:
        points = []

        if self.l2 > self.l1:
            points.append(LineSegment(self.m1, self.m1 + self.a).cuts_bounds(bounds))
            points.append(LineSegment(self.m2, self.m2 + self.a).cuts_bounds(bounds))
        else:
            points.append(LineSegment(self.m, self.m + self.a).cuts_bounds(bounds))

        return points
