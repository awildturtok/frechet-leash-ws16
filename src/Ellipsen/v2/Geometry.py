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


class Vector:
    def __init__(self, x: float, y: float):
        self.x = float(x)
        self.y = float(y)
        self.l = math.sqrt(x ** 2 + y ** 2)

    def __str__(self):
        return "(%s,%s)" % (self.x, self.y)

    # Vector Arithmetic

    def __add__(self, other):
        x = self.x + other.x
        y = self.y + other.y
        return Vector(x, y)

    def __neg__(self):
        return Vector(-self.x, -self.y)

    def __sub__(self, other):
        x = self.x - other.x
        y = self.y - other.y
        return Vector(x, y)

    def __mul__(self, s: float):
        x = self.x * s
        y = self.y * s
        return Vector(x, y)

    def __abs__(self) -> float:
        return self.l

    def __eq__(self, other):
        return self.x == other.x and self.y == other.y

    def __ne__(self, other):
        return not self.__eq__(other)

    def d(self, other) -> float:  # distance to other point
        d = other - self
        return d.l

    def norm(self):  # normalized vector (lenght = 1)
        return self * (1 / self.l)

    def cross_product(self, other):
        return self.x * other.y - self.y * other.x

    def angle(self) -> float:
        return math.atan2(self.y, self.x)

    def rotate(self, a: float):
        x = self.x * math.cos(a) - self.y * math.sin(a)
        y = self.x * math.sin(a) + self.y * math.cos(a)
        return Vector(x, y)

    def rotate_90_l(self):  # rotate 90° to the left
        return Vector(-self.y, self.x)

    def rotate_90_r(self):  # rotate 90° to the right
        return Vector(-self.y, self.x)

    def in_bounds(self, bounds: ((float, float), (float, float))) -> bool:
        return bounds[0][0] <= self.x <= bounds[0][1] and bounds[1][0] <= self.y <= bounds[1][1]


class LineSegment:
    def __init__(self, p1: Vector, p2: Vector):

        if p1 == p2:
            print("Error: Start- and Endpoint are equal.")
            p2 = p2 + Vector(1, 0)

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

    @property
    def __str__(self):
        return str(self.p1) + "->" + str(self.p2) + ": m=" + str(self.m) + " n=" + str(self.n)

    # Line Arithmetic

    def fr(self, r: float) -> Vector:  # point for set parameter r
        return self.p1 + self.d * r

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

    def d_ls_point(self, p: Vector) -> float:  # closest distance of p to line segment
        d_l = self.d_l_point(p)  # closest distance of p to line
        r = self.r_point(p + self.d.rotate_90_l().norm() * d_l)  # parameter r for this distance
        if 0 <= r <= 1:  # if the projection of p lies on the segment return this distance
            return d_l
        else:
            return min(self.p1.d(p), self.p2.d(p))


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

    def p(self, t: float) -> Vector:  # point on ellipse for set parameter t
        return ((self.a * math.cos(t)) + (self.b * math.sin(t))) + self.m
