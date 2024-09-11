from ngsolve import *
import numpy as np
import matplotlib.pyplot as plt

import sys
sys.argv = ["fun"]

mesh = Mesh (unit_square.GenerateMesh(maxh=0.2))
myfunc = x*(1-x)-3
Draw(myfunc, mesh, "cfu");

X = np.linspace(0, 1, num=11)
Y = np.ones_like(X) * 0.2
myfunc(mesh(X, Y))
plt.plot(X, myfunc(mesh(X, Y)))
plt.xlabel('x')
#plt.show()

fun=x*(1-x)-3
fes = H1(mesh, order=1)
gf_fun = GridFunction(fes)
gf_fun.Set(fun)
a = gf_fun
Draw(a, mesh, 'a')


