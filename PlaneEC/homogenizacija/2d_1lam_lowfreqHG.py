from ngsolve import *
from netgen.occ import *

import netgen.gui

size= 0.00082 #0.00145 #0.00046  #0.00156 #0.0008
d=0.001*3
h=0.006
core = MoveTo(-d/2, -h/2).Rectangle(d,h).Face()
#core.edges.name="Gama"
#core.faces.maxh=0.0002
core.edges.Max(X).name = "r"
core.edges.Min(X).name = "l"
core.edges.Min(Y).name = "b"
core.edges.Max(Y).name = "t"
#air.edges.Min(X).maxh=0.1


core.faces.name="core"
core.faces.col = (1, 1, 0)  #colour

geo = Glue([core])

#####################################

mesh = Mesh(OCCGeometry(geo, dim=2).GenerateMesh(maxh=size, quad_dominated=True))

fsU = HCurl(mesh, order=0, dirichlet="b|l|t|r", complex=True, nograds = False)
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

#sig = mesh.MaterialCF({ "core|Km" : 1 }, default=0.000000001)
#B0=200*x+0.3
B0 = 0.4
omega=314*2
mu0 = 1.257e-6

rel = 1000
rho= 5e-7 

kh= d/4 #d/2/sqrt(5) #size # d/4

term1 = rho*grad(csp)*grad(theta)*dx + 1j*omega*curl(mvp)*theta*dx
term1hg= 1j*omega*(1)/(12*rel)*(kh**2) * grad(csp)*grad(theta)*dx
term2 = - 1j*omega* rel*curl(mvp)*curl(alpha)*dx + 1j*omega*csp*curl(alpha)*dx
term2hg= 1/12*(kh*omega)**2 /rho *curl(mvp)*curl(alpha)*dx
#term2 = - csp*curl(alpha)*dx + rel*curl(mvp)*curl(alpha)*dx 

a = BilinearForm(term1+term2+term1hg+term2hg)
a.Assemble()

force = -1j * omega*B0*theta*dx
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
print('Pfun=',Pfun)
#print('delta=',delta)

#PPPPPPPPPPPPPPPPPPPP

rot=CF( (0 , 1,  -1, 0), dims=(2,2) )
J = rot*grad(Tpot)
B=curl(Apot)

Power=0.5*rho*J*Conj(J) + 0.5* 1/12*(kh*omega)**2 /rho * (B0+B)*Conj(B0+B)
print('PowLoss =', Integrate(Power,mesh,order=5).real)
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#print(Integrate(1,mesh,element_wise = True))

#-------------------------------------
f_Km = IfPos( (x-d/2)*(x), 0, IfPos( (y-4*h/12)*(y-h/2), 0, 1))
f_core = mesh.MaterialCF({ "core" : 1 }, default=0.0)
reg_Km =f_Km
#pozivom metode mesh.Materials() dobiva se objekt klase ngsolve.comp.Region
defon = mesh.Materials('Km')
#int_Km = Integrate(CF(1)*reg_Km, mesh, definedon=defon)
int_Km = Integrate(CF(1)*reg_Km, mesh)

#B_bar= Integrate((B+B0)*reg_Km, mesh, definedon=defon)/int_Km
B_bar= Integrate((B+B0)*reg_Km, mesh)/int_Km
J_bar= Integrate(J*CF((0,1))*reg_Km, mesh)/int_Km

print('B_bar = ',B_bar)
print('J_bar = ', J_bar)



Draw (f_core, mesh, "f_core")
Draw (reg_Km, mesh, "reg_Km")


Draw(Tpot, mesh, "T")
Draw(Apot, mesh, "A")
Draw (B, mesh, "B")
Draw (J, mesh, "J")
Draw(grad(Tpot),mesh,"gradT")


Jfuntilda=Jfun-J_bar#*CF((0,1))
Jtilda=J-J_bar*CF((0,1))
ggint= Integrate(Jfuntilda*Jfuntilda*reg_Km, mesh)/int_Km/B0**2
print('ggint=',ggint)
print('rho*ggint=',rho*ggint)

Pow_ggint=Integrate(Jfuntilda*Conj(Jfuntilda)*reg_Km, mesh)/int_Km/B0**2
print('Pow_ggint=',Pow_ggint)

import numpy as np
import matplotlib.pyplot as plt
import sys
sys.argv = ["dummy"]

X = np.linspace(-0.5*d, 0.5*d, num=21)
Y = np.zeros_like(X)+0.00065

plt.plot(X, (J).imag(mesh(X, Y)))
plt.plot(X, (Jfun).imag(mesh(X, Y)))
plt.xlabel('x')
#plt.show()



