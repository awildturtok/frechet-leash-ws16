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

from Geometry import *
import plotly.plotly as py
import plotly.graph_objs as go


class OneLineSegment(LineSegment):  # one of the two input line segments
    def __init__(self, ls: LineSegment, s: Vector):
        super().__init__(ls.p1, ls.p2)

        self.intersects = self.contains_point(s)  # intersects S?
        self.dir = (self.r_point(s) <= 1)  # points in direction  of S?
        self.rs = self.r_point(s)  # S in parametric terms of line segment

        '''# calculate offset to S
        if self.dir:
            self.offset = self.p1.d(s)
        else:
            self.offset = self.p2.d(s)
        if self.intersects:
            self.offset = -self.offset'''

    def __str__(self):
        return str(self.p1) + "->" + str(self.p2) + ": m=" + str(self.m) + " n=" + str(self.n) + \
               "\n   Intersects:" + str(self.intersects) + " Direction:" + str(self.dir) + \
               " r(S): " + str(self.rs)  # + " Offset:" + str(self.offset)


def intersection(a: LineSegment, b: LineSegment):  # calculates the intersection point of two line segments
    if math.isinf(a.m) and not math.isinf(b.m):
        x = a.n
    elif math.isinf(b.m) and not math.isinf(a.m):
        x = b.n
    elif (a.m - b.m) != 0:
        x = (b.n - a.n) / (a.m - b.m)
    else:
        x = float("nan")
    return Vector(x, a.fx(x))

class Cell(Ellipse):
    def __init__(self, m, a, b, c, d):
        pass

class TwoLineSegments:  # Input: two line segments
    def __init__(self, a: LineSegment, b: LineSegment):
        # calculate shortest and longest possible line length - FALSCH: ex. 0,0->0,2 und 1,1->2,2
        self.bounds_l = (min(a.p1.d(b.p1), a.p1.d(b.p2), a.p2.d(b.p1), a.p2.d(b.p2)),
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
        return "two Line Segments:\nLine Segment A: " + str(self.a) + "\nLine Segment B: " + str(
            self.b) + "\nIntersect (" + str(self.intersect) + "): " + str(self.s)

    def cell(self) -> Cell:
        angle_a = self.a.d.angle()
        vec_a = self.a.d.rotate(-angle_a)
        vec_b = self.b.d.rotate(-angle_a)
        norm_a = a.norm()
        norm_b = b.norm()
        c = (norm_b + norm_a) * (0.5 / (norm_b - norm_a).l)
        d = (norm_b - norm_a) * (0.5 / (norm_b + norm_a).l)

        print("normA: " + str(norm_a))
        print("normB: " + str(norm_b))
        print("c: " + str(c))
        print("d: " + str(d))


st1 = LineSegment(Vector(0, 0), Vector(0, 1))
st2 = LineSegment(Vector(0, 0), Vector(-1, 1))
eingabe = TwoLineSegments(st1, st2)
eingabe.cell()
print(eingabe)


### Neu schreiben:
# Gibt Ellipsen-Ausgabe für eine Zelle aus (n1: Anzahl der Ellipsen, n2: Punkt-Genauigkeit der Ellipsen)
def zellenAusgabe(eingabe, n1, n2):
    # Punkte in Zelle berechnen
    ellipsen = []
    data = []
    for i1 in range(n1):
        ellipse = eingabe.ellipse(eingabe.minl + (foat(i1) / (n1 - 1)) * (eingabe.maxl - eingabe.minl))
        ellipsen.append(ellipse)
        x1 = []
        x2 = []
        y1 = []
        y2 = []
        for i2 in range(n2 + 1):
            x = (float(i2) / n2)
            aus1 = ellipse.aus1(x)
            aus2 = ellipse.aus2(x)
            if (0 <= aus1 <= 1):
                x1.append(x)
                y1.append(aus1)
            if (0 <= aus2 <= 1):
                x2.append(x)
                y2.append(aus2)
        data.append([[x1, y1], [x2, y2], ellipse.l])
    return data


# Ellipsen-Daten einer Zelle mit Plotly visualisieren
def zelleZuPlotly(ausgabe):
    data = []
    for a in ausgabe:
        aus1 = a[0]
        aus2 = a[1]
        l = a[2]
        trace1 = go.Scatter(x=aus1[0], y=aus1[1], name=str(l) + 'a')
        trace2 = go.Scatter(x=aus2[0], y=aus2[1], name=str(l) + 'b')
        data.append(trace1)
        data.append(trace2)
    fig = dict(data=data)
    py.plot(fig, filename='Ellipsen-Test')
