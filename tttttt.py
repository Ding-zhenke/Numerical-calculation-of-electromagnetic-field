# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 18:30:16 2020

@author: 31214
"""
from alldcc import *
import sympy as sy
c = CoordSys3D('c')
#直角坐标系
#建立一个矢量函数
def fun(x,y,z):
    v = x*c.i+y*c.j+z*c.k
    return v
#带入具体点（x,y,z)
def sub(a,x,y,z):
    va=a.subs([(c.x,x),
               (c.y,y),
               (c.z,z),
               (sy.pi,3.14)])
    return va
hm,k = sy.symbols("Hm,k",positive=True)
#
E = fun(20*sy.cos(1e9*sy.pi*t),0,0)
D = None
H = None
B = None
J = None
Jd = None
p0 = None
sig1= None
so = solve_method(E,D,H,B,J,Jd,p0,t,sig1)
so.allsub(er=80,sigr=4)
E=so.E
D=so.D
H=so.H
B=so.B
J=so.J
Jd=so.Jd
p0=so.p0
