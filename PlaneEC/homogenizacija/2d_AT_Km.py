from ngsolve import *
from netgen.occ import *


size=0.00002
d=0.001
h=0.01

Km = MoveTo(0.0, -h/20).Rectangle(d/2,h/10).Face()
Km.faces.name = "Km"
Km.edges.Max(X).name = "kmr"
Km.edges.Min(X).name = "kml"
Km.edges.Min(Y).name = "kmb"
Km.edges.Max(Y).name = "kmt"
Km.faces.maxh=size

geo=Km

mesh = Mesh(OCCGeometry(geo, dim=2).GenerateMesh(maxh=size, quad_dominated=False))

fsU = HCurl(mesh, order=0, dirichlet="kmb|kml|kmt|kmr", complex=True, nograds = False)
fsV = H1(mesh, order=1, dirichlet="kml|kmr", complex=True)
fes=fsU*fsV
mvp, csp = fes.TrialFunction()
alpha, theta = fes.TestFunction()

print('mesh.nv',mesh.nv)
print('mesh.nedge', mesh.nedge)
print('mesh.ne', mesh.ne)

sol = GridFunction(fes)
Apot, Tpot = sol.components

#B0=200*x+0.3
B0 = 0.4716
J0 = 4203317*CF((0,1))
omega=314
mu0 = 1.257e-6
rel = 100
rho= 5e-7 

rot=CF( (0 , 1,  -1, 0), dims=(2,2) )

term1 = rho * grad(csp)*grad(theta)*dx + 1j * omega*curl(mvp)*theta*dx 
term2 = - 1j*omega* rel*curl(mvp)*curl(alpha)*dx + 1j*omega* csp*curl(alpha)*dx #+ 0.1*mvp*alpha*dx 
#term1 = rho * grad(csp)*grad(theta)*dx + 1j * omega*(mvp)*(rot*grad(theta))*dx 
#term2 = - 1j*omega* rel*curl(mvp)*curl(alpha)*dx + 1j*omega* (rot*grad(csp))*(alpha)*dx #+ 0.1*mvp*alpha*dx 

a = BilinearForm(term1+term2)
a.Assemble()

force =-1j * omega*B0*theta*dx -1j*omega* J0*alpha*dx 
f=LinearForm(force)
f.Assemble()

r = f.vec - a.mat * sol.vec
sol.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs()) * r

J = rot*grad(Tpot)
B=curl(Apot)
A=Apot*CF((1,0))


Power=0.5*rho*J*Conj(J)
print('Power losses =', Integrate(Power,mesh,order=5).real)


Draw (B, mesh, "B")
Draw (J, mesh, "J")
Draw(Tpot, mesh, "T")
Draw(Apot, mesh, "A")

###############

f_Km = IfPos( (x)*(x-d/2), 0, IfPos( (y-h/20)*(y+h/20), 0, 1))
f_core = mesh.MaterialCF({ "core" : 1 }, default=0.0)
reg_Km =f_Km
#pozivom metode mesh.Materials() dobiva se objekt klase ngsolve.comp.Region
defon = mesh.Materials('Km')
int_Km = Integrate(CF(1)*reg_Km, mesh, definedon=defon)
#int_Km = Integrate(CF(1)*reg_Km, mesh)

B_bar= Integrate(B*reg_Km, mesh, definedon=defon)/int_Km
J_bar= Integrate(Apot*CF((0,1))*reg_Km, mesh, definedon=defon)/int_Km

print('B_bar = ',B_bar)
print('J_bar = ', J_bar)


#Draw (f_core, mesh, "f_core")
#Draw (reg_Km, mesh, "reg_Km")





import numpy as np
import matplotlib.pyplot as plt
import sys
sys.argv = ["dummy"]

X = np.linspace(0.0, 0.5*d, num=21)
Y = np.zeros_like(X)


plt.plot(X, Tpot.imag(mesh(X, Y)))
plt.plot(X, B.real(mesh(X, Y)))
plt.xlabel('x')
plt.show()



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



