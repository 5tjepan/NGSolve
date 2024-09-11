

from ngsolve import *
from netgen.occ import *
import netgen.gui
import matplotlib.pyplot as plt
import sys
sys.argv = ["fun"]

outer= Circle((0,0), 0.2).Face()
outer.edges.name = 'rub'

inner = MoveTo(-0.08,-0.08).Line(0.16,0.0).Line(0,0.16).Line(-0.16,0).Close().Face()
#inner = MoveTo(-0.8,-0.8).Line(1.6,0).Line(0,1.6).Line(-1.6,0).Close().Face()
#inner = MoveTo(-0.8,-0.8).Line(1.6,-0.1).Line(0,1.7).Line(-1.6,-0.2).Close().Face()
#inner = MoveTo(0.0,-1.0).Line(1,1).Line(-1,1).Line(-1,-1).Close().Face()

inner.edges.name="interface"
inner.faces.maxh=0.01
outer = outer - inner

inner.faces.name="inner"
inner.faces.col = (1, 1, 0)  #colour
outer.faces.name="outer"

geo = Glue([inner, outer])
#Draw(OCCGeometry(geo));

mesh = Mesh(OCCGeometry(geo, dim=2).GenerateMesh(maxh=0.05, quad_dominated=False))
print('mesh.nv',mesh.nv)
print('mesh.nedge', mesh.nedge)
print('mesh.ne', mesh.ne)


fes = HCurl(mesh, order=0, dirichlet="rub",  complex=True, nograds = False) #CMPLX
u, v = fes.TnT()

bb=CF((0.1*y,-0.1*x))
#bb=CF((10,10))

gfu = GridFunction(fes)
old = GridFunction(fes)
#dirich = GridFunction(fes)

#gfu.Set(bb, VOL_or_BND = BND)
gfu.Set(bb, definedon=mesh.Boundaries('rub'))
#old.Set(bb, definedon=mesh.Boundaries('rub'))
#dirich.Set(bb, definedon=mesh.Boundaries('rub'))

omega=314
#rel = 1200
sigma = 1e2
sig = mesh.MaterialCF({ "inner" : 1 }, default=None)
#sigma=CoefficientFunction( (10000, 0,  0, 10000), dims=(2,2) )

mu0=4*pi*1e-7
Href=[0,42,53,62,70,79,88,100,113,132,157,193,255,376,677,1624,1e9]
Bref=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1e9*mu0]
nu=[Href[i]/Bref[i] if Bref[i]!=0 else 0 for i in range(len(Bref))]
#HBcurve = BSpline(2, [0]+list(Href), list(Bref)) #ovo bi trebala biti instanca klase BSpline
#diffHBcurve= HBcurve.Differentiate() #munonlin.Differentiate() metoda daje objekt klase BSpline
nucurve = BSpline(2, [0]+list(Bref), list(nu)) #ovo bi trebala biti instanca klase BSpline

rel = (100*curl(gfu)+ 1e-5).Norm()

from vh import myNewtonSolver


term1 = (1/mu0)*curl(u)*curl(v)*dx('outer') + rel*curl(u)*curl(v)*dx('inner') + \
    +1j*omega*sigma*u*v*dx('inner') + 1j*1e-4*u*v*dx('outer') 

a = BilinearForm(term1)

dummy=CF((1e-17,1e-17))
rhs1 = dummy*v*dx
f = LinearForm(rhs1)

myNewtonSolver(gfu,a,f, dampfactor=0.5, Nmax=100)

#####POSTPROCESING
A=gfu
E = - 1j * omega * A
J = - 1j * omega * sigma * A *sig
B = curl(gfu)

Pow=0.5*E*Conj(J) 
Peddy=Integrate(Pow, mesh, order=5)
print('Peddy',Peddy)

Draw(A, mesh, "A")
Draw (B, mesh, "B")
Draw (J, mesh, "J")


