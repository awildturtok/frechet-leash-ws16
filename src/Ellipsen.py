# -*- coding: utf-8 -*-

import math

import plotly.plotly as py
import plotly.graph_objs as go

class P:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def move(self, dp):
        self.x = self.x + dp.x
        self.y = self.y + dp.y

    def __str__(self):
        return "(%s,%s)"%(self.x, self.y)

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
        return self.l()

    def d(self, other):
        dx = self.x - other.x
        dy = self.y - other.y
        return math.sqrt(dx**2 + dy**2)

    def l(self):
        return math.sqrt((self.x)**2 + (self.y)**2)

class S:
    def __init__(self, p1, p2):
        self.p1 = p1
        self.p2 = p2
        self.d = p2-p1
        self.l = self.d.l()
        if (p2.x-p1.x) != 0:
            self.m = (p2.y-p1.y)/(p2.x-p1.x)
            self.n = p1.y-self.m*p1.x
        else:
            self.m = float("inf")
            self.n = p1.x
        
    def __str__(self):
        return str(self.p1)+"->"+str(self.p2)
        
    def fr(self, r):
        return self.p1+self.d*r
    
    def fx(self, x):
        if self.m != float("inf"): return self.m*x+self.n
        else: return float("nan")
    
    def rx(self, x):
        if self.m != float("inf"): return (x-self.p1.x)/(self.p2.x-self.p1.x)
        else: return float("nan")
    
    def ry(self, y):
        if self.m != 0: return (y-self.p1.y)/(self.p2.y-self.p1.y)
        else: return float("nan")

def winkel(s, a, b):
    a = a-s
    b = b-s
    return math.acos((a.x*b.x + a.y*b.y)/(a.l()*b.l()))

class D:
    def __init__(self, a, b):
        self.a = a
        self.b = b
        self.s = self.schnitt()
        self.schneidenA = (a.m==float("inf") or 0<=a.rx(self.s.x)<=1) and (a.m==0 or 0<=a.ry(self.s.y)<=1)
        self.schneidenB = (b.m==float("inf") or 0<=b.rx(self.s.x)<=1) and (b.m==0 or 0<=b.ry(self.s.y)<=1)
        self.schneiden = self.schneidenA and self.schneidenB
        self.richtungA = (a.m==float("inf") or a.rx(self.s.x)<=1) and (a.m==0 or a.ry(self.s.y)<=1)
        self.richtungB = (b.m==float("inf") or b.rx(self.s.x)<=1) and (b.m==0 or b.ry(self.s.y)<=1)

        if self.richtungA: self.loa = self.s.d(self.a.p1)
        else: self.loa = self.s.d(self.a.p2)
        if self.richtungB: self.lob = self.s.d(self.b.p1)
        else: self.lob = self.s.d(self.b.p2)

        if self.schneidenA: self.loa = -self.loa
        if self.schneidenB: self.lob = -self.lob

        if self.schneiden: self.minl = 0.0
        else: self.minl = min(a.p1.d(b.p1), a.p1.d(b.p2), a.p2.d(b.p1), a.p2.d(b.p2))
        self.maxl = max(a.p1.d(b.p1), a.p1.d(b.p2), a.p2.d(b.p1), a.p2.d(b.p2))

        if self.s != self.a.p2 and self.s != self.b.p2: self.w = winkel(self.s, self.a.p2, self.b.p2)
        elif self.s == self.a.p2: self.w == winkel(self.s, self.a.fr(2), self.b.p2)
        elif self.s == self.a.p2: self.w == winkel(self.s, self.a.p2, self.b.fr(2))
        else: self.w = winkel(self.s, self.a.fr(2), self.b.fr(2))

    def __str__(self):
        return str(self.a)+" und "+str(self.b)
        
    def schnitt(self):
        if (self.a.m-self.b.m) != 0: x = (self.b.n-self.a.n)/(self.a.m-self.b.m)
        else: x = float("nan")
        return P(x,self.a.fx(x))
    
    def ellipse(self, l):
        return E(self, l)

class E:
    def __init__(self, d, l):
        self.d = d
        self.l = l
    
    def aus(self, v, r):
        l = self.l
        loa = self.d.loa
        lob = self.d.lob
        la = self.d.a.l
        lb = self.d.b.l
        w = self.d.w
        r = r
        s = float("nan")

        if (not self.d.richtungA and self.d.richtungB) or (not self.d.richtungA and not self.d.richtungB):
            r = 1 - r

        tmp = lb**2*(loa**2*math.cos(w)**2-loa**2+2*loa*la*r*math.cos(w)**2-2*loa*la*r+la**2*r**2*math.cos(w)**2-la**2*r**2+l**2)
        if tmp >= 0: s = (v*math.sqrt(tmp)+loa*lb*math.cos(w)+la*lb*r*math.cos(w)-lob*lb)/(lb**2)

        if (self.d.richtungA and not self.d.richtungB) or (not self.d.richtungA and not self.d.richtungB):
            s = 1 - s

        return s

    def aus1(self, r):
        return self.aus(1, r)

    def aus2(self, r):
        return self.aus(-1, r)


st1 = S(P(1.0,-1.0),P(2.0,-2.0))
st2 = S(P(0.0,0.0),P(1.0,1.0))
d = D(st1,st2)
print ("Eingabe: "+str(d))
print ("Strecke 1: "+str(st1)+": m="+str(st1.m)+" n="+str(st1.n))
print ("Strecke 2: "+str(st2)+": m="+str(st2.m)+" n="+str(st2.n))
print ("Schnittpunkt: "+str(d.s))
print ("Winkel: "+str(d.w))
print ("Cos(Winkel): "+str(math.cos(d.w)))
print ("Offset A: "+str(d.loa))
print ("Offset B: "+str(d.lob))
print ("Min. Länge: "+str(d.minl))
print ("Max. Länge: "+str(d.maxl))
print ("Schneiden(A,B): "+str(d.schneiden)+"("+str(d.schneidenA)+","+str(d.schneidenB)+")")
print ("Richtung(A,B): "+str(d.richtungA)+","+str(d.richtungB))

ellipsen = []
n1 = 10
n2 = 250
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
    trace1 = go.Scatter(x=x1, y=y1, name=str(ellipse.l)+'a')
    trace2 = go.Scatter(x=x2, y=y2, name=str(ellipse.l)+'b')
    data.append(trace1)
    data.append(trace2)

fig = dict(data=data)
py.plot(fig, filename='Ellipsen-Test')


