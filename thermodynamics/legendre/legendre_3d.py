import numpy as np
import sympy as smp

from manim import *

# Define some constants
a = 0.01
b = 0.005
c = 150

# Define thermodynamic variables and potentials
s, v = smp.symbols("s v", real=True)
# symbolic U(S,V)
u = a * s ** 2 + b * (v - c) ** 2
# use u_f internal as a function that takes entropy and v as arguments
u_f = smp.lambdify(args=[s, v], expr=u)
# get the partial derivatives and their physical variables
duds = smp.diff(u, s)
t = duds
dudv = smp.diff(u, v)
p = - dudv
# use legendre transform to obtain G(T,P)
g = u - t * s + p * v

u_s = u.subs(v, 10)
u_v = u.subs(s, 10)
