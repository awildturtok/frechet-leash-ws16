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
import plotly.plotly as py
import plotly.graph_objs as go

class P: # Punkt Klasse
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.l = self.length()
    def __str__(self):
        return "(%s,%s)"%(self.x, self.y)

    # Vektor Arithmetik
    def __add__(self, other):
        x = self.x + other.x
        y = self.y + other.y
        return P(x,y)
    def __neg__(self):
        return P(-self.x, -self.y)
    def __sub__(self, other):
        x = self.x - other.x
        y = self.y - other.y
        return P(x,y)
    def __mul__(self, s):
        x = self.x*s
        y = self.y*s
        return P(x,y)
    def __abs__(self):
        return self.l

    def length(self): # Berechnet die Distanz (0,0)->P
        return math.sqrt((self.x)**2 + (self.y)**2)

    def d(self, other): # Berechnet die Distanz zwischen 2 Punkten
        d = other-self
        return d.length()

class S: # Strecken Klasse
    def __init__(self, p1: P, p2: P):
        self.p1 = p1 # Anfangspunkt P2
        self.p2 = p2 # Endpunkt P2
        self.d = p2-p1 # Differenz-Vektor
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

    #Strecken Arithmetik
    def fr(self, r): # Punkt auf Geraden für geg. Parameter r
        return self.p1+self.d*r
    def fx(self, x): # y-Wert für geg. x-Wert
        if self.m != float("inf"): return self.m*x+self.n
        else: return float("nan")
    def rx(self, x): # Parameter r für geg. x-Wert
        if self.m != float("inf"): return (x-self.p1.x)/(self.p2.x-self.p1.x)
        else: return float("nan")
    def ry(self, y): # Parameter r für geg. y-Wert
        if self.m != 0: return (y-self.p1.y)/(self.p2.y-self.p1.y)
        else: return float("nan")
    def enthP(self, p): # Liegt P auf der Strecke?
        return (self.m == float("inf") or 0 <= self.rx(p.x) <= 1) and (self.m == 0 or 0 <= self.ry(p.y) <= 1)
    def rP(self, p): # Parameter r für geg. Punkt P (Ausgabe nur sinnvoll wenn P auf S liegt)
        if not math.isinf(self.m): rx = self.rx(p.x)
        if not self.m == 0: ry = self.ry(p.y)
        if math.isinf(self.m): rx = ry
        if self.m == 0: ry = rx
        return 0.5*(rx+ry)


def cosWinkel(s, a, b): # cos-Wert des Winkels des Dreicks ASB
    a = a-s
    b = b-s
    return (a.x*b.x + a.y*b.y)/(a.l*b.l)

class D: # 2-Strecken-Konstellation Klasse
    def __init__(self, a, b):
        self.a = a # Strecke A
        self.b = b # Strecke B
        self.s = self.schnitt() # Schnittpunkt
        self.schnittA = a.enthP(self.s) # Liegt der Schnittpunkt auf A?
        self.schnittB = b.enthP(self.s) # Liegt der Schnittpunkt auf B?
        self.schneiden = self.schnittA and self.schnittB # Liegt der Schnittpunkt auf A und B?
        self.richtA = (a.rP(self.s) <= 1) # Zeigt A in Richtung S?
        self.richtB = (b.rP(self.s) <= 1) # Zeigt B in Richtung S?
        # Berechnet den Offset loa von A
        if self.richtA: self.loa = self.s.d(self.a.p1)
        else: self.loa = self.s.d(self.a.p2)
        # Berechnet den Offset lob von B
        if self.richtB: self.lob = self.s.d(self.b.p1)
        else: self.lob = self.s.d(self.b.p2)
        # Entscheidet das Vorzeichen der Offsets
        if self.schnittA: self.loa = -self.loa
        if self.schnittB: self.lob = -self.lob
        # Berechnet die kürzeste und längste mögliche Leinenlänge
        if self.schneiden: self.minl = 0.0
        else: self.minl = min(a.p1.d(b.p1), a.p1.d(b.p2), a.p2.d(b.p1), a.p2.d(b.p2))
        self.maxl = max(a.p1.d(b.p1), a.p1.d(b.p2), a.p2.d(b.p1), a.p2.d(b.p2))
        # Berechent den Cos-Wert des Winkels zwischen den Geraden
        if self.s != self.a.p2 and self.s != self.b.p2: self.cosw = cosWinkel(self.s, self.a.p2, self.b.p2)
        elif self.s == self.a.p2: self.cosw = cosWinkel(self.s, self.a.fr(2), self.b.p2)
        elif self.s == self.a.p2: self.cosw = cosWinkel(self.s, self.a.p2, self.b.fr(2))
        else: self.cosw = cosWinkel(self.s, self.a.fr(2), self.b.fr(2))
    def __str__(self):
        return str(self.a)+" und "+str(self.b)

    def schnitt(self): # Berechnet den Schnittpunkt der 2 Geraden
        if math.isinf(self.a.m) and not math.isinf(self.b.m): x = self.a.n
        elif math.isinf(self.b.m) and not math.isinf(self.a.m): x = self.b.n
        elif (self.a.m-self.b.m) != 0:
            x = (self.b.n-self.a.n)/(self.a.m-self.b.m)
        else: x = float("nan")
        return P(x,self.a.fx(x))

    def ellipse(self, l): # Gibt eine Ellipsen-Objekt für eine geg. Leinenlänge zurück
        return E(self, l)

class E: # Ellipsen Klasse
    def __init__(self, d, l):
        self.d = d # 2-Strecken-Konstellation
        self.l = l # Leinenlänge

    def aus(self, v, r): # Ausgabe-Parameter s (auf Strecke B) für geg. Vorzeichen der Wurzel und Parameter r (auf Stecke A)
        l = self.l # Leinenlänge
        loa = self.d.loa # Offset loa der Strecke A
        lob = self.d.lob # Offset lob der Strecke B
        la = self.d.a.l # Länge der Strecke A
        lb = self.d.b.l # Länge der Strecke B
        cosw = self.d.cosw # Cos-Wert des Winkels zwischen den Geraden
        s = float("nan") # Ausgabe-Parameter s
        if not self.d.richtA: r = 1 - r # Die Zählrichtung von r ist Abhängig von der Richtung von A

        # Parameter s für geg. Werte berechen
        tmp = lb**2*(loa**2*cosw**2-loa**2+2*loa*la*r*cosw**2-2*loa*la*r+la**2*r**2*cosw**2-la**2*r**2+l**2)
        if tmp >= 0: s = (v*math.sqrt(tmp)+loa*lb*cosw+la*lb*r*cosw-lob*lb)/(lb**2)

        # Checking for direction of the {Strecken}
        if not self.d.richtB: s = 1 - s # Die Zählrichtung von s ist Abhängig von der Richtung von B

        return s

    # Ausgabe-Parameter s für positives (aus1) und negatives (aus2) Vorzeichen der Wurzel
    def aus1(self, r):
        return self.aus(1, r)
    def aus2(self, r):
        return self.aus(-1, r)

def zelle(d, n1, n2):
    print ("Eingabe: "+str(d))
    print ("Strecke 1: "+str(st1)+": m="+str(st1.m)+" n="+str(st1.n))
    print ("Strecke 2: "+str(st2)+": m="+str(st2.m)+" n="+str(st2.n))
    print ("Schnittpunkt: "+str(d.s))
    print ("Winkel: "+str(math.acos(d.cosw)))
    print ("Cos(Winkel): "+str(d.cosw))
    print ("Offset A: "+str(d.loa))
    print ("Offset B: "+str(d.lob))
    print ("Min. Länge: "+str(d.minl))
    print ("Max. Länge: "+str(d.maxl))
    print ("Schneiden(A,B): "+str(d.schneiden)+"("+str(d.schnittA)+","+str(d.schnittB)+")")
    print ("Richtung(A,B): "+str(d.richtA)+","+str(d.richtB))

    ellipsen = []
    data = []
    for i1 in range(n1):
        ellipse = d.ellipse(d.minl+(float(i1)/(n1-1))*(d.maxl-d.minl))
        ellipsen.append(ellipse)
        x1 = []
        x2 = []
        y1 = []
        y2 = []
        for i2 in range(n2+1):
            x = (float(i2)/n2)
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

def plotlyTest(d, n1, n2):
    aus = zelle(d, n1, n2)
    data = []
    for a in aus:
        aus1 = a[0]
        aus2 = a[1]
        l = a[2]
        trace1 = go.Scatter(x=aus1[0], y=aus1[1], name=str(l)+'a')
        trace2 = go.Scatter(x=aus2[0], y=aus2[1], name=str(l)+'b')
        data.append(trace1)
        data.append(trace2)
    fig = dict(data=data)
    py.plot(fig, filename='Ellipsen-Test')

st1 = S(P(1.0, 0.0), P(2.0, 0.0))
st2 = S(P(0.0, 0.0), P(0.0, 1.0))
d = D(st1,st2)
plotlyTest(d, 50, 100)