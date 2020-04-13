# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 22:31:26 2020
例题2.5.3
@author: 31214
"""
#直角坐标系内
from sympy.vector import CoordSys3D, Del
import sympy as sy
from alldcc import *
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
#调用
#常量定义
# em,hm= sy.symbols("Em,Hm")
# #H =C.i*hm*sy.cos(w*t-k*C.z)
# #已知量申明,C.x,C.y
# E = fun(em*sy.cos(2*sy.pi*1e6*t),0,0)

# #计算过程
# D = ep*E
# jd= dt(D)
# j = sig*E
# k = jd.magnitude()/j.magnitude()
# k = sy.cancel(k)
# k = substitude(k,81,1,4)
F =sy.symbols("F",cls = sy.Function)
h = fun(F(c.x,c.y,c.z),F(c.x,c.y,c.z),F(c.x,c.y,c.z))
e = sy.exp(-1/ep*sig*t)*(st(sy.exp(sig/ep*t)*rot(h)))
H_ =st(1/mu * rot(e)) 



