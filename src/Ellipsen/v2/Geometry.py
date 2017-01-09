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
        self.norm = math.sqrt(x**2 + y**2)

    def __str__(self):
        return "(%s,%s)" % (self.x, self.y)

    # Vector Arithmetic

    def __add__(self, other: Vector) -> Vector:
        x = self.x + other.x
        y = self.y + other.y
        return Vector(x,y)

    def __neg__(self) -> Vector:
        return Vector(-self.x, -self.y)

    def __sub__(self, other: Vector) -> Vector:
        x = self.x - other.x
        y = self.y - other.y
        return Vector(x,y)

    def __mul__(self, s: float) -> Vector:
        x = self.x*s
        y = self.y*s
        return Vector(x,y)

    def __abs__(self) -> float:
        return self.norm

    def distance(self, other: Vector) -> float:
        d = other-self
        return d.norm


class LineSegment:

    def __init__(self, p1: Vector, p2: Vector):

        if p1 == p2:
            print("Error: Start- and Endpoint are equal.")
            p2 = p2 + Vector(1,0)

        self.p1 = p1 # Startpoint P1
        self.p2 = p2 # Endpoint P2
        self.d = p2-p1 # Difference Vector
        self.l = self.d.l # Länge der Strecke

        # Steigung m und Y-Abschnitt n berechnen
        if (p2.x-p1.x) != 0:
            self.m = (p2.y-p1.y)/(p2.x-p1.x)
            self.n = p1.y-self.m*p1.x
        else:
            self.m = float("inf")
            self.n = p1.x

    def __str__(self):
        return str(self.p1)+"->"+str(self.p2)

    #Line Arithmetic

    def fr(self, r: float) -> Vector: # point for set parameter r
        return self.p1+self.d*r

    def fx(self, x: float) -> float: # y-Value for set x-Value
        if not math.isinf(self.m): return self.m*x+self.n
        else: return float("nan")

    def fy(self, y: float) -> float: # x-Value for set y-Value
        if math.isinf(self.m): return self.n
        elif self.m != 0: return (y-self.n)/self.m
        else: return float("nan")

    def rx(self, x: float) -> float: # parameter r for set y-Value
        if not math.isinf(self.m): return (x-self.p1.x)/(self.p2.x-self.p1.x)
        else: return float("nan")

    def ry(self, y: float) -> float: # parameter r for set x-Value
        if self.m != 0: return (y-self.p1.y)/(self.p2.y-self.p1.y)
        else: return float("nan")

    def containsP(self, p: Vector) -> bool: # does point p lie on the line segment
        return (self.m == float("inf") or 0 <= self.rx(p.x) <= 1) and (self.m == 0 or 0 <= self.ry(p.y) <= 1)

    def rP(self, p: Vector) -> float: # parameter r for set point
        if not math.isinf(self.m): rx = self.rx(p.x)
        if not self.m == 0: ry = self.ry(p.y)
        if math.isinf(self.m): rx = ry
        if self.m == 0: ry = rx
        return 0.5*(rx+ry)


class Ellipse:

    def __init__(self, m, a, b, angle):
        self.m = m # midpoint
        self.a = a # axis vector a
        self.b = b # axis vector b

    def p(self, t: float) -> Vector: # point on ellipse for set parameter t
        return self.m + math.cos(t)*self.a + math.sin(t)*self.b

