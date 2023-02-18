import numpy as np
import plotly.graph_objects as go
import sympy as smp

# a = 0.01
# x = np.linspace(0, 100, num=101)
# s = 2 * a * x
# X, S = np.meshgrid(x, s)
# U = a * np.power(X, 2)
# G = U - X * S
#
# fig = go.Figure(data=[go.Surface(z=G, x=x, y=s)])
# fig.show()

a = 0.01
b = 0.005
c = 150
# s = np.linspace(0, 100, num=101)
# v = np.linspace(0, c, num=101)
# S, V = np.meshgrid(s, v)
# U = a * np.power(S, 2) + b * np.power(V - c, 2)
# G = U - X * S

s, v = smp.symbols("s v", real=True)
u = a * s ** 2 + b * (v - c) ** 2
u_f = smp.lambdify(args=[s, v], expr=u)
duds = smp.diff(u, s)
t = duds
dudv = smp.diff(u, v)
p = - dudv
g = u - t * s + p * v

u_s = u.subs(v, 10)
u_v = u.subs(s, 10)
print(g)

# fig = go.Figure(data=[go.Surface(z=U, x=s, y=v)])
# fig.update_layout(
#     scene=dict(
#         xaxis_title="S",
#         yaxis_title="V",
#         zaxis_title="U"
#     )
# )
# fig.show()
