from ngsolve import *
from netgen.occ import *


size=0.00004 #0.0008
d=0.001#*3
h=0.01
core = MoveTo(-d/2, -h/2).Rectangle(d,h).Face()
#core.edges.name="Gama"
#core.faces.maxh=0.0002
core.edges.Max(X).name = "r"
core.edges.Min(X).name = "l"
core.edges.Min(Y).name = "b"
core.edges.Max(Y).name = "t"
#air.edges.Min(X).maxh=0.1


Km = MoveTo(0.0, -h/20).Rectangle(d/2,h/10).Face()
Km.faces.name = "Km"
Km.edges.Max(X).name = "kmr"
Km.edges.Min(X).name = "kml"
Km.edges.Min(Y).name = "kmb"
Km.edges.Max(Y).name = "kmt"
Km.faces.maxh=size
#core -= Km

core.faces.name="core"
core.faces.col = (1, 1, 0)  #colour

geo = Glue([core])

#####################################

mesh = Mesh(OCCGeometry(geo, dim=2).GenerateMesh(maxh=size, quad_dominated=True))

fsU = HCurl(mesh, order=0, dirichlet="b|l|t|r|kmr", complex=True, nograds = False)
fsV = H1(mesh, order=1, dirichlet="l|r|kmr", complex=True)
fes=fsU*fsV
mvp, csp = fes.TrialFunction()
alpha, theta = fes.TestFunction()

print('mesh.nv',mesh.nv)
print('mesh.nedge', mesh.nedge)
print('mesh.ne', mesh.ne)

sol = GridFunction(fes)
Apot, Tpot = sol.components

#::::::::::::::::::::::::::::::::::::::

sig = mesh.MaterialCF({ "core|Km" : 1 }, default=0.000000001)
#B0=200*x+0.3
B0 = 0.4
omega=314
mu0 = 1.257e-6

rel = 1000 
rho= 5e-7 


term1 = rho * grad(csp)*grad(theta)*dx + 1j * omega*curl(mvp)*theta*dx #+ 0.1*mvp*alpha*dx
term2 = - 1j*omega* rel*curl(mvp)*curl(alpha)*dx + 1j*omega* csp*curl(alpha)*dx 
#term2 = - csp*curl(alpha)*dx + rel*curl(mvp)*curl(alpha)*dx 

a = BilinearForm(term1+term2)
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

Power=0.5*rho*J*Conj(J)
print('PowLoss =', Integrate(Power,mesh,order=5).real)


#-------------------------------------
f_Km = IfPos( (x-d/2)*(x+d/2), 0, IfPos( (y-h/8)*(y+h/8), 0, 1))
f_core = mesh.MaterialCF({ "core" : 1 }, default=0.0)
reg_Km =f_Km
#pozivom metode mesh.Materials() dobiva se objekt klase ngsolve.comp.Region
#defon = mesh.Materials('Km')
#int_Km = Integrate(CF(1)*reg_Km, mesh, definedon=defon)
int_Km = Integrate(CF(1)*reg_Km, mesh)

#B_bar= Integrate((B+B0)*reg_Km, mesh, definedon=defon)/int_Km
B_bar= Integrate((B+B0)*reg_Km, mesh)/int_Km
J_bar= Integrate(J*CF((0,1))*reg_Km, mesh)/int_Km

print('B_bar = ',B_bar)
print('J_bar = ', J_bar)

#-------------------------------
#b_abs = B.Norm()
#b_sqr=(B*Conj(B)).Norm()
b_abs=B+B0
b_sqr=(B+B0)*(B+B0)
b_bar= Integrate(b_abs*reg_Km, mesh)/int_Km
int_bsqr = Integrate(b_sqr*reg_Km, mesh)
FF_avg=(1/int_Km * int_bsqr/(b_bar**2))
print('FF=1/Km * int_bsqr/bbar^2= rel* =',FF_avg)

#j_sqr=(J*Conj(J)).Norm()
j_sqr=J*J
int_jsqr = Integrate(j_sqr*reg_Km, mesh)
GG_avg =(1/int_Km * int_jsqr/(b_bar**2))
print('GG=1/Km * int_jsqr/jbar^2= rho* =',GG_avg)

nuzz=rel*FF_avg - 1j/omega*rho * GG_avg
print('nuzz',nuzz)


babs_bar= Integrate((B+B0).Norm()*reg_Km, mesh)/int_Km
jabs_sqr=(J.Norm())**2
int_jabs_sqr = Integrate(jabs_sqr*reg_Km, mesh)
GGabs_avg =(1/int_Km * int_jabs_sqr/babs_bar**2)

Pcoeff=0.5*rho * GGabs_avg #(1/int_Km * abs(int_jsqr)/((abs(b_bar))**2))
print('Pcoeff =', Pcoeff)
print('P_Km = Pcoeff*b_bar^2 * int_Km=', Pcoeff*babs_bar**2 * int_Km)

tok= Integrate((B+B0).Norm()*f_core, mesh)
print('tok',tok)

print('b_bar(Km)=',b_bar)


Draw (f_core, mesh, "f_core")
Draw (reg_Km, mesh, "reg_Km")


#Draw(rel,mesh, "rel")

Draw(Tpot, mesh, "T")
Draw(Apot, mesh, "A")
Draw (B, mesh, "B")
Draw (J, mesh, "J")


""" Jfuntilda=Jfun-J_bar#*CF((0,1))
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
Y = np.zeros_like(X)

plt.plot(X, (J).real(mesh(X, Y)))
plt.plot(X, (Jfun).real(mesh(X, Y)))
plt.xlabel('x')
plt.show() """



