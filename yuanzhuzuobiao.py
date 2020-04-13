# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 16:01:58 2020
圆柱坐标系
作业2.14.2
@author: 31214
"""
#圆柱坐标系(i,j,k)=(p,phi,z)
from sympy.vector import CoordSys3D, Del
import sympy as sy
from alldcc import *
c = CoordSys3D('c',
               transformation='cylindrical',
               variable_names=("r", "ph", "z"),
               vector_names=("er","eph","ez"))
#建立一个矢量函数
def fun(i,j,k):
    v = i*c.er+j*c.eph+k*c.ez
    return v
#带入具体点（r,phi,z)
def sub(a,r,phi,z):
    va=a.subs([(c.r,r),(c.ph,phi),(c.z,z),(sy.pi,3.14)])
    return va
'''
创建一个向量(rho,phi,z)
v = fun(F(c.r),F(c.ph),F(c.z))
'''
#常量定义
a =sy.symbols("a")
#
E = None
D = None
H = fun(a*c.r,0,0)
B = None
J = None
Jd = None
p0 = None
sig1= 0
so = solve_method(E,D,H,B,J,Jd,p0,t,sig1)
E=so.E
D=so.D
H=so.H
B=so.B
J=so.J
Jd=so.Jd
p0=so.p0