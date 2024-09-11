from ngsolve import *
from netgen.occ import *

d=0.001*10
h=0.01
core = MoveTo(-d/2, -h/2).Rectangle(d,h).Face()
#core = MoveTo(-0.08,-0.08).Line(0.16,0.0).Line(0,0.16).Line(-0.16,0).Close().Face()
#core.edges.name="Gama"
#core.faces.maxh=0.0002
core.edges.Max(X).name = "r"
core.edges.Min(X).name = "l"
core.edges.Min(Y).name = "b"
core.edges.Max(Y).name = "t"
#air.edges.Min(X).maxh=0.1

core.faces.name="core"
core.faces.col = (1, 1, 0)  #colour

geo = core
#Draw(OCCGeometry(geo));

mesh = Mesh(OCCGeometry(geo, dim=2).GenerateMesh(maxh=0.001, quad_dominated=False))

fsU = HCurl(mesh, order=0, dirichlet="b|l|t|r", complex=True, nograds = False)
fsV = H1(mesh, order=2, dirichlet="b|l|t|r", complex=True)
#fsV = H1(mesh, order=1, dirichlet="interface", definedon="core", complex=True)
fes=fsU*fsV
mvp, csp = fes.TrialFunction()
alpha, theta = fes.TestFunction()

print('mesh.nv',mesh.nv)
print('mesh.nedge', mesh.nedge)
print('mesh.ne', mesh.ne)

sol = GridFunction(fes)
Apot, Tpot = sol.components

sig = mesh.MaterialCF({ "core" : 1 }, default=0.000000001)
#B0=200*x+0.3
B0 = 0.4
omega=1
mu0 = 1.257e-6
#rho = 0.00001

#rel = mesh.MaterialCF({ "core|Km" : 32000 }, default=(1/mu0))
#rho = mesh.MaterialCF({ "core|Km" : 6e-7 }, default=10)
rel = 100000 
rho= 5e-7 

#sigma=CoefficientFunction( (0, 0,  0, 10), dims=(2,2) )
#rel=CoefficientFunction( ( 1200, 0, 0, 1200), dims=(2,2) )

""" term1 = rho * grad(csp)*grad(theta)*dx('core') + 1j*omega*curl(mvp)*theta*dx('core')
term2 = rel*curl(mvp)*curl(alpha)*dx - csp*curl(alpha)*dx('core')
term3 = 1*mvp*alpha*dx('air') """


term1 = rho * grad(csp)*grad(theta)*dx + 1j * omega*curl(mvp)*theta*dx
term2 = - csp*curl(alpha)*dx + rel*curl(mvp)*curl(alpha)*dx 

a = BilinearForm(term1+term2)
a.Assemble()

force = -1j * omega*B0*theta*dx
f=LinearForm(force)
f.Assemble()

r = f.vec - a.mat * sol.vec
sol.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs()) * r

rot=CF( (0 , 1,  -1, 0), dims=(2,2) )
J = rot*grad(Tpot)
B=curl(Apot)

Power=0.5*rho*J*Conj(J)
print('Power losses =', Integrate(Power,mesh,order=5).real)


###############

f_Km = IfPos( (x)*(x-d/2), 0, IfPos( (y-h/10)*(y+h/10), 0, 1))
f_core = mesh.MaterialCF({ "core" : 1 }, default=0.0)
reg_Km =f_Km
#pozivom metode mesh.Materials() dobiva se objekt klase ngsolve.comp.Region
defon = mesh.Materials('core')
#int_Km = Integrate(CF(1)*reg_Km, mesh)
int_Km = Integrate(CF(1)*reg_Km, mesh, definedon=defon)

#b_abs = B.Norm()
#b_sqr=(B*Conj(B)).Norm()
""" b_abs=B
b_sqr=B*B
b_bar= Integrate(b_abs*reg_Km, mesh, definedon=defon)/int_Km
int_bsqr = Integrate(b_sqr*reg_Km, mesh, definedon=defon)
FF_avg=(1/int_Km * int_bsqr/(b_bar**2))
print('FF=1/Km * int_bsqr/bbar^2= rel* =',FF_avg)

#j_sqr=(J*Conj(J)).Norm()
j_sqr=J*J
int_jsqr = Integrate(j_sqr*reg_Km, mesh, definedon=defon)
GG_avg =(1/int_Km * int_jsqr/(b_bar**2))
print('GG=1/Km * int_jsqr/jbar^2= rho* =',GG_avg)

nuzz=rel*FF_avg - 1j/omega*rho * GG_avg
print('nuzz',nuzz)


babs_bar= Integrate(B.Norm()*reg_Km, mesh, definedon=defon)/int_Km
jabs_sqr=(J.Norm())**2
int_jabs_sqr = Integrate(jabs_sqr*reg_Km, mesh, definedon=defon)
GGabs_avg =(1/int_Km * int_jabs_sqr/babs_bar**2)

Pcoeff=0.5*rho * GGabs_avg #(1/int_Km * abs(int_jsqr)/((abs(b_bar))**2))
print('Pcoeff =', Pcoeff)
print('P_Km = Pcoeff*b_bar^2 * int_Km=', Pcoeff*babs_bar**2 * int_Km)

tok= Integrate(B.Norm()*f_core, mesh, definedon=defon)
print('tok',tok)

print('b_bar(Km)=',b_bar) """


Draw (f_core, mesh, "f_core")
Draw (reg_Km, mesh, "reg_Km")


#Draw(rel,mesh, "rel")

Draw(Tpot, mesh, "T")
Draw(Apot, mesh, "A")
Draw (B, mesh, "B")
Draw (J, mesh, "J")



import numpy as np
import matplotlib.pyplot as plt
import sys
sys.argv = ["dummy"]

X = np.linspace(-0.5*d, 0.5*d, num=21)
Y = np.zeros_like(X)

plt.plot(X, J.Norm()(mesh(X, Y)))
plt.xlabel('x')
plt.show()

