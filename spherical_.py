# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 16:59:45 2020
球坐标系
例题2.4.1
@author: 31214
"""

import sympy as sy
from sympy.vector import CoordSys3D, Del
from alldcc import *
#建立一个矢量函数(rho,theta,phi)
def fun(i,j,k):
    v = i*c.er+j*c.eth+k*c.eph
    return v
#带入具体点（r,theta,phi)
def sub(a,r,theta,phi):
    va=a.subs([(c.r,r),(c.ph,phi),(c.th,theta),(sy.pi,3.14)])
    return va
c = CoordSys3D('c',
               transformation='spherical',
               vector_names=("er","eth","eph"),
               variable_names=("r", "th", "ph"))
#常量定义
A,a =sy.symbols("A,a")
#已知量
E = None
D = None
H = fun(0,0,a*c.r)
B = None
J = None
Jd = None
p0 = None
sig1= 0
so = solve_method(E,D,H,B,J,Jd,p0,t,sig1)

