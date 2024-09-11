from ngsolve import *
from netgen.occ import *

import netgen.gui

size= 0.0009 #0.0118 #0.0059 #0.00145 #0.00046  #0.00156 #0.0008
d=0.001*3
h=0.006
core1 = MoveTo(-d/2, -h/2).Rectangle(d/4,h).Face()
core1.edges.Min(X).name = "l1"
core1.edges.Min(Y).name = "b1"
core1.edges.Max(Y).name = "t1"
#air.edges.Min(X).maxh=0.1
core2 = MoveTo(-d/4, -h/2).Rectangle(d/4,h).Face()
core2.edges.Min(Y).name = "b2"
core2.edges.Max(Y).name = "t2"

core3 = MoveTo(0, -h/2).Rectangle(d/4,h).Face()
core3.edges.Min(Y).name = "b3"
core3.edges.Max(Y).name = "t3"
#air.edges.Min(X).maxh=0.1

core4 = MoveTo(d/4, -h/2).Rectangle(d/4,h).Face()
core4.edges.Max(X).name = "r4"
core4.edges.Min(Y).name = "b4"
core4.edges.Max(Y).name = "t4"

core1.faces.name="core1"
core2.faces.name="core2"
core3.faces.name="core3"
core4.faces.name="core4"
core1.faces.col = (1, 1, 0)  #colour
core3.faces.col = (1, 0.3, 0)  #colour

geo = Glue([core1,core2,core3, core4])

#####################################

mesh = Mesh(OCCGeometry(geo, dim=2).GenerateMesh(maxh=size, quad_dominated=False))

fsU = HCurl(mesh, order=0, dirichlet="l1|t1|t2|t3|t4|b4|b3|b2|b1|r4", complex=True, nograds = False)
fsV = H1(mesh, order=1, dirichlet="l1|r4", complex=True)
#fsV = H1(mesh, order=1, dirichlet="l1|r4|t1|t2|t3|t4|b4|b3|b2|b1", complex=True)
fes=fsU*fsV
mvp, csp = fes.TrialFunction()
alpha, theta = fes.TestFunction()

print('mesh.nv',mesh.nv)
print('mesh.nedge', mesh.nedge)
print('mesh.ne', mesh.ne)

sol = GridFunction(fes)
Apot, Tpot = sol.components


#sig = mesh.MaterialCF({ "core|Km" : 1 }, default=0.000000001)
#B0=6e2*x+0.4
#B0=6e5*(x**2-d**2/12)+0.4
B0 = 0.4
Draw(B0,mesh,'B0')

omega=314*8
mu0 = 1.257e-6

rel = 1000 
rho= 5e-7 

kh= d/4 #d/2/sqrt(5) #size # d/4

rot=CF( (0 , 1,  -1, 0), dims=(2,2) )
diry=CF((0,1))

#::::::::::::::::::::::::::::::::::::::
fesl2=L2(mesh, order=0)
xavg=GridFunction(fesl2)

arealist=Integrate(1,mesh,element_wise=True)
intlist=Integrate(x, mesh, element_wise=True)

#for i in range(len(intlist)):   xavg.vec[i]=intlist[i]/arealist[i]
xavg.Set(x)

#xavg = CF({"core1" : -0.001125}, {"core2" : -0.000375}, {"core3" : 0.000375},{"core4" : 0.001125})

#xavg=IfPos(x, IfPos((x-d/4),0.001125,0.000375),IfPos((x+d/4),-0.000375,-0.001125))

#xbump=IfPos((y-5*d/9)*(y+d/3),1,0)
bump=IfPos((x-1*d/7)*(x+d/3),0,1) * IfPos((0.5*y-2*d/7)*(y+d/3),0,1) *sin(100*x)
Draw(bump,mesh,'bump')

ay=CF((0,1))
gfun= (x-1.00*xavg)*ay *(-1j)*omega/rho * bump #*10*Norm(x)**0.5  
ffun= -(x-1.00*xavg)/rel * bump  #vjerojatno je minus!
#::::::::::::::::::::::::::::::::::::::


#PAZI NA 24 UMJESTO 12
#term1 = rho*grad(csp)*grad(theta)*dx + 1j*omega*curl(mvp)*theta*dx
term1 = rho*grad(csp)*grad(theta)*dx + 1j*omega*mvp*(rot*grad(theta))*dx
#term2 = - 1j*omega* rel*curl(mvp)*curl(alpha)*dx + 1j*omega*csp*curl(alpha)*dx
term2 = - 1j*omega*rel*curl(mvp)*curl(alpha)*dx + 1j*omega*alpha*(rot*grad(csp))*dx
#term2 = - 1j*omega* rel*curl(mvp)*curl(alpha)*dx + 1j*omega*csp*curl(alpha)*dx

term1hg= 1j*omega*(1)/(12*rel)*(kh**2) * grad(csp)*grad(theta)*dx
term2hg= 1/12*(kh*omega)**2 /rho *curl(mvp)*curl(alpha)*dx

term1a= 1j*omega*theta*ffun* rot*grad(csp)*CF((0,1))*dx #('core1|core2|core3')
#term1a= 1j*omega*ffun*theta*Norm(rot*grad(csp))*dx
term2a= 1j*omega*gfun*alpha*curl(mvp)*dx   #('core1|core2|core3')

term1b= 1j*(rho*(gfun*gfun))*curl(alpha)*curl(mvp)*dx   #('core1|core2|core3')
term2b= -(1j)*omega*ffun*ffun*(rot*(grad(theta)))*(rot*(grad(csp)))*dx  #('core1|core2|core3')

term4= 1j*omega*gfun*gfun*curl(alpha)*(mvp)*dx - 1j*omega*gfun*(alpha)*curl(mvp)*dx
#term1a= -rho*gfun*gfun*curl(alpha)*curl(mvp)*dx
#term2a= (1j)*omega*rel*ffun*ffun*grad(theta)*(grad(csp))*dx

#a = BilinearForm(term1 - term2 + 0*term1a + 0*term2a + term1b + 5*term2b )
a = BilinearForm(term1 + term2 + 12*term1a + 2*term2a + 15*term1b + term2b)
#a = BilinearForm(term1 + term2 + 1*term1a + 1*term2a )
#a = BilinearForm(term1 - term2 + term1hg - term2hg)
#a = BilinearForm(term1 + term2 + term1hg + term2hg)
a.Assemble()

force = -1j * omega*B0*theta*dx
#force = -1j * omega*B0*theta*dx - rho*gfun*gfun*curl(alpha)*B0*dx - 2*1j*omega*gfun*alpha*B0*dx
#force = -1j * omega*B0*theta*dx - rho*gfun*gfun*curl(alpha)*B0*dx
#force = -1j * omega*B0*theta*dx - 1*1j*omega*gfun*alpha*B0*dx
###force = -1j * omega*B0*theta*dx - 1*1j*omega*Norm(gfun*alpha)*B0*dx
f=LinearForm(force)
f.Assemble()

r = f.vec - a.mat * sol.vec
sol.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs()) * r

#...................
delta=1/sqrt(omega/(rho*rel*2))
gama = (1+1j)/delta
Bfun = B0*gama*d/2 *(cosh(gama*x)/sinh(gama*d/2))
Jfun = -1j*omega*d/(2*rho)* B0 *(sinh(gama*x)/sinh(gama*d/2))

#Draw(Bfun, mesh, 'Bfun')

Pfun=(h*d)*(B0*omega)**2 *d*delta/(8*rho) * ( (sinh(d/delta)-sin(d/delta))/(cosh(d/delta)-cos(d/delta)) )
#print('Pfun=',Pfun)
print('Pfun =', Integrate(Pfun/(h*d),mesh,order=5).real)
#print('delta=',delta)

#PPPPPPPPPPPPPPPPPPPP

J = rot*grad(Tpot)
B=curl(Apot)

#Power=0.5*rho*(J + gfun*(B0+B)) * Conj(J + gfun*(B0+B))


Power=0.5*rho*J*Conj(J) #+ 0.5*rho*(gfun*(B0+B))*Conj((B0+B)*gfun)
#Power=0.5*rho*J*Conj(J) + 0.5* 1/12*(kh*omega)**2 /rho * (B0+B)*Conj(B0+B)


#Power=0.5*rho*J*Conj(J) + 0.5* 1/12*(kh*omega)**2 /rho * (B0+B)*Conj(B0+B) \
#    + 0*0.5*omega*kh/3 * ((J*diry)*Conj(B)).imag #0.5*omega*kh/3 * (Norm(J.imag)*Norm(B.real) - Norm(J.real)*Norm(B.imag))

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


Draw(Tpot, mesh, "T")
Draw(Apot, mesh, "A")
Draw (B, mesh, "B")
Draw (J, mesh, "J")
#Draw(grad(Tpot),mesh,"gradT")


Jfuntilda=Jfun-J_bar#*CF((0,1))
Jtilda=J-J_bar*CF((0,1))
ggint= Integrate(Jfuntilda*Jfuntilda*reg_Km, mesh)/int_Km/B0**2
#print('ggint=',ggint)
#print('rho*ggint=',rho*ggint)

Pow_ggint=Integrate(Jfuntilda*Conj(Jfuntilda)*reg_Km, mesh)/int_Km/B0**2
#print('Pow_ggint=',Pow_ggint)


""" import numpy as np
import matplotlib.pyplot as plt
import sys
sys.argv = ["dummy"]

X = np.linspace(-0.5*d, 0.5*d, num=100)
Y = np.zeros_like(X)+0.00065

plt.plot(X, gfun.imag(mesh(X, Y)))
plt.plot(X, gfun.imag(mesh(X, Y)))
plt.xlabel('x')
#plt.show() """


""" fesl2=L2(mesh, order=0)
gfun=GridFunction(fesl2)

intlist=Integrate(x, mesh, element_wise=True)

for i in range(len(intlist)): gfun.vec[i]=intlist[i] """

Draw(gfun, mesh, "gfun")

