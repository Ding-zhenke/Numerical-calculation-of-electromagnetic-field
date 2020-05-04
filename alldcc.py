# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 17:55:05 2020
This is dcc's surpports  
@author: 31214

from sympy.vector import CoordSys3D, Del
C = CoordSys3D('C')
delop = Del()
v = C.x*C.y*C.z * (C.i + C.j + C.k)
a=delop.cross(v, doit = True)
#v ^ C.i  点乘

magnitude()幅度
Returns the magnitude of this vector.
normalize()单位向量
Returns the normalized version of this vector.
outer(other)
Returns the outer product of this vector with another, in the form of a Dyadic instance.
projection(other, scalar=False)投影
Returns the vector or scalar projection of the ‘other’ on ‘self’.
    from sympy.vector.coordsysrect import CoordSys3D
    from sympy.vector.vector import Vector, BaseVector
    C = CoordSys3D('C')
    i, j, k = C.base_vectors()
    v1 = i + j + k
    v2 = 3*i + 4*j
    v1.projection(v2)
    7/3*C.i + 7/3*C.j + 7/3*C.k
    v1.projection(v2, scalar=True)
    7/3
to_matrix(system)
Returns the matrix form of this vector with respect to the specified coordinate system.
v.to_matrix(哪一个坐标系)
"""
from sympy.vector import CoordSys3D, Del
import sympy as sy
delop = Del()
#基础函数
def rot(v):#旋度
    delop = Del()
    rot_v=delop.cross(v, doit = True)
    rot_v=rot_v.doit()
    return rot_v
def div(v):#散度
    delop = Del()
    div_v=delop.dot(v,doit=True)
    div_v=div_v.doit()
    return div_v
def grad(s):#梯度
    delop = Del()
    g=delop.gradient(s,doit=True)
    g=g.doit()
    return g
#公式调用
def dt(x):#对T偏导
    dx = sy.Derivative(x,t)
    dx=dx.doit()
    return sy.vector.vector.VectorAdd(dx)
def st(x):#对T积分
    sx = sy.integrate(x,t)
    sx = sx.doit()
    return sy.vector.vector.VectorAdd(sx)
def i_conv(j,t):
    p = div(j)
    p = st(p)
    return p #电流连续性方程
def e_div(e,e0):
    p = div(e)*e0#电场散度
    return p#电荷密度
def e_rot(e,t):#0
    e = rot(e)
    e = st(-e)
    return e
def b_div(b):#0
    return div(b)
def b_rot(b,mu0):#总电流密度，包含磁化电流密度
    return rot(b)/mu0
def d_div(d):#返回的事电荷密度
    return div(d)
def h_rot(h):#电流密度
    return rot(h)

'''
电磁场基本常量，pi就是派
t:时间
p0:电荷密度
ep:介电常数
si:电导率
mu:磁导率
'''
t,w,mu,ep,ep0,mu0,sig=sy.symbols("t,omega,mu,epsilon,epsilon0,mu0,sigma",
                               positive=True)
#电磁场变量 
E,H,B,D,J=sy.symbols("E,H,B,D,J")
q,p0=sy.symbols("q,rho0",positive=True)
'''
E:电场强度
H：磁场强度
B：磁感应强度
D：电通量
q：电荷量>=0
J: 电流密度
p0:电荷密度>=0
'''
class solve_method():
    def __init__(self,E,D,H,B,J,Jd,p0,t,sig1):
        self.sig_value = sig1
        if E != None:
            print("已知E")
            self.E=E
            #先求D
            self.D  = ep*E 
            #再求Jd
            self.Jd = dt(self.D)
            #求p0
            self.p0 = div(self.D)
            #推出B
            self.B  = st(-rot(E))
            #推出H
            self.H  = 1/mu*self.B
            #自然推出 J
            if sig1 !=0:#就用第一方程算，说明此时无其他未知量了
                if self.H != 0*self.H:                
                    self.J  = rot(self.H)-self.Jd
                else:
                    self.J = sig*self.E
            else:
                self.J = 0 #可用第一方程推出未知数   
        elif D != None:
            print("已知D")
            self.D = D
            self.E  = D*(1/ep)
            self.Jd = dt(D)
            self.p0 = div(self.D)
            self.B  = st(-rot(self.E))
            self.H  = 1/mu*self.B1
            #自然推出 J
            if sig1 !=0:#就用第一方程算，说明此时无其他未知量了                
                if self.H != 0*self.H:                
                    self.J  = rot(self.H)-self.Jd
                else:
                    self.J = sig*self.E
            else:
                self.J = 0 #可用第一方程推出未知数
        elif H != None:
            self.H = H
            print("已知H")
            self.B=H*mu
            if div(self.B)!=0:
                print('divB doesn’t equal to 0!')
            if sig1 == 0 or J == 0:#能不能把传导电流消除
                self.D  = st(rot(H))
                self.E  = 1/ep*self.D
                self.Jd = dt(self.D)
                self.p0 = div(self.D)
                self.J = 0
            else:#看来是消除不了
                print("当作是均匀导体来解")
                #对第一方程解方程
                self.E  = sy.exp(-1/ep*sig*t)
                #这一步我省略的常数，不知道有影响不
                self.E  = self.E*(st(sy.exp(1/ep*sig*t)*rot(H)))
                self.D  = ep*self.E
                self.p0 = div(self.D)
                self.J  = self.E*sig
                self.Jd = dt(self.D)
        elif B!= None:
            print("已知B")
            self.B=B
            self.H = 1/ep*B
            if div(self.B)!=0:
                print('divB doesn’t equal to 0!')
            if sig1 == 0 or J == 0:#能不能把传导电流消除
                self.D  = st(div(self.H))
                self.E  = 1/ep*self.D
                self.Jd = dt(self.D)
                self.p0 = div(self.D)
            else:#看来是消除不了
                print("当作是均匀导体来解")
                #对第一方程解方程
                self.E  = sy.exp(-1/ep*sig*t)
                #这一步我省略的常数，不知道有影响不
                self.E  = self.E*(st(sy.exp(1/ep*sig*t)*rot(self.H)))
                self.D  = ep*self.E
                self.p0 = div(self.D)
                self.J  = self.E*sig
            self.Jd = dt(self.D)
        elif J != None:
            self.J=J
            print("已知传导电流J")
            self.E  = 1/sig*J#此时不需要sig了
            self.D  = self.E*ep
            self.Jd = dt(self.D)
            self.p0 = div(self.D)
            self.B  = st(-rot(self.E))
            self.H  = 1/mu * self.B
        elif Jd !=None:
            self.Jd=Jd
            self.D = st(Jd)
            self.E = 1/ep*self.D
            self.p0 = div(self.D)
            self.B  = st(-rot(self.E))
            self.H  = 1/mu * self.B
            self.J = rot(self.H)-Jd
        else:
            print("error!")
        self.E = self.E.subs(sig,sig1)
        self.D = self.D.subs(sig,sig1)
        self.H = self.H.subs(sig,sig1)
        self.B = self.B.subs(sig,sig1)
        if self.J !=0:     
            self.J = self.J.subs(sig,sig1)
        self.Jd = self.Jd.subs(sig,sig1)
        self.p0 = self.p0.subs(sig,sig1)
    #如果需要数值替换的话
    def substitude(self,a,er=1,mur=1,sigr=1):
        #ep = ep0*epr
        a = a.subs(ep,er*8.85*1e-12)
        #mu = mu0*mur
        a = a.subs(mu,mur*4*sy.pi*1e-7)
        #pi
        a = a.subs(sy.pi,3.14)
        #sigma
        a = a.subs(sig,sigr)
        return a
    def allsub(self,er=1,mur=1,sigr=1):
        self.E = self.substitude(self.E,er=1,mur=1,sigr=1)
        self.D = self.substitude(self.D,er=1,mur=1,sigr=1)
        self.H = self.substitude(self.H,er=1,mur=1,sigr=1)
        self.B = self.substitude(self.B,er=1,mur=1,sigr=1)
        if self.J !=0:     
            self.J = self.substitude(self.J,er=1,mur=1,sigr=1)
        self.Jd = self.substitude(self.Jd,er=1,mur=1,sigr=1)
        self.p0 = self.substitude(self.p0,er=1,mur=1,sigr=1)
            
        
            

























