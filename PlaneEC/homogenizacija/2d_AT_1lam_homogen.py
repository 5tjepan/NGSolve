from ngsolve import *
from netgen.occ import *


size= 0.00008 #0.00039 #0.0006
d=0.001
h=0.01
core = MoveTo(-d/2, -h/2).Rectangle(d,h).Face()
#core.edges.name="Gama"
#core.faces.maxh=0.0002
core.edges.Max(X).name = "r"
core.edges.Min(X).name = "l"
core.edges.Min(Y).name = "b"
core.edges.Max(Y).name = "t"
#air.edges.Min(X).maxh=0.1

air = MoveTo(-d/1.5, -h/1.8).Rectangle(2*d/1.5, 2*h/1.8).Face()
air.edges.name="rub"
air-=core

air.faces.name="air"
core.faces.name="core"
core.faces.col = (0.3, 1, 0)  #colour

geo = Glue([air, core])

#####################################

mesh = Mesh(OCCGeometry(geo, dim=2).GenerateMesh(maxh=size, quad_dominated=True))

fsU = HCurl(mesh, order=0, dirichlet="rub", complex=True, nograds = False)
fsV = H1(mesh, order=1, dirichlet="l|r", complex=True)
fes=fsU*fsV
mvp, csp = fes.TrialFunction()
alpha, theta = fes.TestFunction()

print('mesh.nv',mesh.nv)
print('mesh.nedge', mesh.nedge)
print('mesh.ne', mesh.ne)

sol = GridFunction(fes)
Apot, Tpot = sol.components

#::::::::::::::::::::::::::::::::::::::

#sig = mesh.MaterialCF({ "core" : 1 }, default=0.0000012)
#B0=200*x+0.3
B0 = 0.4
#J0=62700*1j
omega=314
mu0 = 1.257e-6

#rel = 0#1000 
hgrel=1000.5441573734199+52.257037207118685j
#rho= 5e-7

rel = mesh.MaterialCF({ "core" : 1000 }, default=795775)
rho = mesh.MaterialCF({ "core" : 5e-7 }, default=1)

term1 = rho * grad(csp)*grad(theta)*dx + 1j*omega*curl(mvp)*theta*dx 
term2 = -1j*omega* rel *curl(mvp)*curl(alpha)*dx + 1j*omega*csp*curl(alpha)*dx
term3 = 0.1*mvp*alpha*dx
#termHG= -1j*omega* hgrel *curl(mvp)*curl(alpha)*dx
#term2 = - csp*curl(alpha)*dx + rel*curl(mvp)*curl(alpha)*dx 

a = BilinearForm(term1+term2+term3)
#a = BilinearForm(termHG+term3)
a.Assemble()

force = -1j *omega*B0*theta*dx
f=LinearForm(force)
f.Assemble()

r = f.vec - a.mat * sol.vec
sol.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs()) * r

#...................
delta=1/sqrt(omega/(rho*rel*2))
gama = (1+1j)/delta
Bfun = B0*gama*d/2 *(cosh(gama*x)/sinh(gama*d/2))
Jfun = -1j*omega*d/(2*rho)* B0 *(sinh(gama*x)/sinh(gama*d/2))

Draw(Bfun, mesh, 'Bfun')

Pfun=(h*d)*(B0*omega)**2 *d*delta/(8*rho) * ( (sinh(d/delta)-sin(d/delta))/(cosh(d/delta)-cos(d/delta)) )
#print('Pfun=',Pfun)
#print('delta=',delta)

#PPPPPPPPPPPPPPPPPPPP

rot=CF( (0 , 1,  -1, 0), dims=(2,2) )

B=curl(Apot)
J = rot*grad(Tpot) #+ B0*(Jfun/B0-J0)
rho=5e-7
Power= 0.5*rho*J*Conj(J) + 8200*rho*0.5*(B0+B)**2
print('Power losses =', Integrate(Power,mesh,order=5).real)


#-------------------------------------

f_Km = IfPos( (x)*(x-d/2), 0, IfPos( (y-h/20)*(y+h/20), 0, 1))
f_core = mesh.MaterialCF({ "core" : 1 }, default=0.0)
reg_Km =f_Km
#pozivom metode mesh.Materials() dobiva se objekt klase ngsolve.comp.Region
defon = mesh.Materials('Km')
#int_Km = Integrate(CF(1)*reg_Km, mesh, definedon=defon)
int_Km = Integrate(CF(1)*reg_Km, mesh)

B_bar= Integrate(B.Norm()*reg_Km, mesh, definedon=defon)/int_Km
J_bar= Integrate(J.Norm()*reg_Km, mesh, definedon=defon)/int_Km

print('B_bar = ',B_bar)
print('J_bar = ', J_bar)

#b_abs = B.Norm()
#b_sqr=(B*Conj(B)).Norm()

Draw (f_core, mesh, "f_core")
Draw (reg_Km, mesh, "reg_Km")

Draw(Tpot, mesh, "T")
Draw(Apot, mesh, "A")
Draw (B, mesh, "B")
Draw (J, mesh, "J")


""" Jtilda=J+J_bar*CF((0,1))
ggint= Integrate(Jtilda*Jtilda*reg_Km, mesh, definedon=defon)/int_Km/B0**2
print('ggint=',ggint)
print('rho*ggint=',rho*ggint) """


import numpy as np
import matplotlib.pyplot as plt
import sys
sys.argv = ["dummy"]

X = np.linspace(-0.5*d, 0.5*d, num=21)
Y = np.zeros_like(X)

#plt.plot(X, J.imag(mesh(X, Y)))
plt.plot(X, Jfun.real(mesh(X, Y)))
plt.xlabel('x')
#plt.show()



