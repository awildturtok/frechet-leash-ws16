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

import Geometry.*
import plotly.plotly as py
import plotly.graph_objs as go


class aLineSegment(LineSegment): # one of the two input line segments
    def __init__(self, ls: LineSegment, s: Vector):

        self = ls

        self.intersects = self.containsP(s) # intersects S?
        self.dir = (self.rP(s) <= 1) # points in direction  of S?

        #calculate offset to S
        if self.dir: self.offset = self.p1.d(s)
        else: self.offset = self.p2.d(s)
        if self.intersects: self.offset = -self.offset

    def __str__(self):
        return str(super())+": "+str(self.intersects)+", "+str(self.dir)+", "+str(self.offset)


class twolineSegments: # Input: two line segments

    def __init__(self, a, b):

        # calculate shortest and longest possible connecting line length
        self.minl = min(a.p1.d(b.p1), a.p1.d(b.p2), a.p2.d(b.p1), a.p2.d(b.p2))
        self.maxl = max(a.p1.d(b.p1), a.p1.d(b.p2), a.p2.d(b.p1), a.p2.d(b.p2))

        self.parallel = self.a.m == self.b.m
        if not self.parallel: # case 1: lines are not parallel
            self.s = self.intersection() # intersection point S

            self.a = aLineSegment(a, self.s)  # line segment a
            self.b = aLineSegment(b, self.s)  # line sqgment b

            self.intersect = self.a.intersects and self.b.instersects # does the intersection point lie on A and B
            if self.intersect: self.minl = 0.0 # if line segments intersect, set minl to 0

            # calculate cos of angle btw. A and B
            '''if self.s != self.a.p2 and self.s != self.b.p2: self.cosw = cosWinkel(self.s, self.a.p2, self.b.p2)
            elif self.s == self.a.p2: self.cosw = cosWinkel(self.s, self.a.fr(2), self.b.p2)
            elif self.s == self.a.p2: self.cosw = cosWinkel(self.s, self.a.p2, self.b.fr(2))
            else: self.cosw = cosWinkel(self.s, self.a.fr(2), self.b.fr(2))'''

        else: # case 2: lines are parallel
            self.dirB = (a.rP(a.p1+b.d) >= 0)  # do A and B point in the same direction
            """self.stWinkel = math.atan(a.m) # Steigungswinkel der Geraden
            # Distanz der 2 Geraden
            if self.a.m == 0 or math.isinf(self.a.m): self.dist = abs(self.a.n - self.b.n)
            else: self.dist = math.sin(self.stWinkel) * abs(a.fy(b.n))
            self.ob = 0 # Offset Strecke B im Vergleich zu A"""
            ### Hier weiter Spezialfall: Parallel implementieren

    def __str__(self):
        return str(self.a)+" und "+str(self.b)

    def intersection(self): # calculates the intersection point
        if math.isinf(self.a.m) and not math.isinf(self.b.m): x = self.a.n
        elif math.isinf(self.b.m) and not math.isinf(self.a.m): x = self.b.n
        elif (self.a.m-self.b.m) != 0:
            x = (self.b.n-self.a.n)/(self.a.m-self.b.m)
        else: x = foat("nan")
        return P(x,self.a.fx(x))

    def ellipse(self, l): # returns ellipse object for given length l
        return Ellipse(self, l)


# Gibt Ellipsen-Ausgabe für eine Zelle aus (n1: Anzahl der Ellipsen, n2: Punkt-Genauigkeit der Ellipsen)
def zellenAusgabe(eingabe, n1, n2):
    #Debugger Log
    print ("Eingabe: " + str(eingabe))
    print ("Strecke A: "+str(eingabe.a)+": m="+str(eingabe.a.m)+" n="+str(eingabe.a.n))
    print ("Strecke B: "+str(eingabe.b)+": m="+str(eingabe.b.m)+" n="+str(eingabe.b.n))
    print ("Leinenlänge (Min., Max.): " + str(eingabe.minl) + ", " + str(eingabe.maxl))
    if eingabe.parallel:
        print ("Spezialfall: Strecken sind parallel.")
        print ("Richtung: " + str(eingabe.dirB))
        print ("Steigungswinkel: atan(" + str(eingabe.a.m) + ") = " + str(eingabe.stWinkel))
        print ("Distanz: " + str(eingabe.dist))
        print ("Offset: " + str(eingabe.ob))
    else:
        print ("Normalfall: Strecken sind nicht parallel.")
        print ("Schnittpunkt: " + str(eingabe.s))
        print ("cos(" + str(math.acos(eingabe.cosw)) + ") = " + str(eingabe.cosw))
        print ("Offset(A,B): " + str(eingabe.oa) + ", " + str(eingabe.ob))
        print ("Schneiden(A,B): " + str(eingabe.schneiden) + "(" + str(eingabe.intersectA) + "," + str(eingabe.intersectB) + ")")
        print("Richtung(A,B): " + str(eingabe.dirA) + "," + str(eingabe.dirB))

    #Punkte in Zelle berechnen
    ellipsen = []
    data = []
    for i1 in range(n1):
        ellipse = eingabe.ellipse(eingabe.minl + (foat(i1) / (n1 - 1)) * (eingabe.maxl - eingabe.minl))
        ellipsen.append(ellipse)
        x1 = []
        x2 = []
        y1 = []
        y2 = []
        for i2 in range(n2+1):
            x = (foat(i2)/n2)
            aus1 = ellipse.aus1(x)
            aus2 = ellipse.aus2(x)
            if (0<=aus1<=1):
                x1.append(x)
                y1.append(aus1)
            if (0<=aus2<=1):
                x2.append(x)
                y2.append(aus2)
        data.append([[x1,y1], [x2,y2], ellipse.l])
    return data

#Ellipsen-Daten einer Zelle mit Plotly visualisieren
def zelleZuPlotly(ausgabe):
    data = []
    for a in ausgabe:
        aus1 = a[0]
        aus2 = a[1]
        l = a[2]
        trace1 = go.Scatter(x=aus1[0], y=aus1[1], name=str(l)+'a')
        trace2 = go.Scatter(x=aus2[0], y=aus2[1], name=str(l)+'b')
        data.append(trace1)
        data.append(trace2)
    fig = dict(data=data)
    py.plot(fig, filename='Ellipsen-Test')

st1 = S(P(0.0, 0.0), P(1.0, 0.0))
st2 = S(P(0.0, 0.0), P(1.0, 1.0))
eing = Eingabe(st1, st2)
zelle = zellenAusgabe(eing, 8, 100)
zelleZuPlotly(zelle)