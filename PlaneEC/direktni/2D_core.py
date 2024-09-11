from ngsolve import *
from netgen.occ import *

outer= Circle((0,0), 0.2).Face()
outer.edges.name = 'rub'

""" outer = MoveTo(0, -1).Rectangle(3,2).Face()
outer.edges.name="rub"
outer.edges.Max(X).name = "r"
outer.edges.Min(X).name = "l"
outer.edges.Min(Y).name = "b"
outer.edges.Max(Y).name = "t"
#outer.edges.Min(X).maxh=0.1 """

inner = MoveTo(-0.08,-0.08).Line(0.16,0.0).Line(0,0.16).Line(-0.16,0).Close().Face()
#inner = MoveTo(-0.8,-0.8).Line(1.6,-0.1).Line(0,1.7).Line(-1.6,-0.2).Close().Face()
#inner = MoveTo(0.0,-1.0).Line(1,1).Line(-1,1).Line(-1,-1).Close().Face()
#inner = MoveTo(-1,-0.5).Line(1.5,-0.5).Line(0.5,1.5).Line(-1.5,0.5).Close().Face()
#inner = MoveTo(0.0,-0.9).Line(1.1,0.9).Line(-1.1,0.9).Line(-1.0,-0.9).Close().Face()
#inner = MoveTo(0.0,-0.9).Line(1.1,0.9).Line(-1.1,0.9).Line(-1.1,-0.9).Close().Face()

inner.edges.name="interface"
inner.faces.maxh=0.01
outer = outer - inner

inner.faces.name="inner"
inner.faces.col = (1, 1, 0)  #colour
outer.faces.name="outer"


geo = Glue([inner, outer])
#Draw(OCCGeometry(geo));


mesh = Mesh(OCCGeometry(geo, dim=2).GenerateMesh(maxh=0.05, quad_dominated=False))
fes = HCurl(mesh, order=0, dirichlet="rub",  complex=True, nograds = False) #CMPLX
gfu = GridFunction(fes)
print('mesh.nv',mesh.nv)
print('mesh.nedge', mesh.nedge)
print('mesh.ne', mesh.ne)

u, v = fes.TnT()

bb=CF((0.1*y,-0.1*x))
#bb=CF((10,10))

#gfu.Set(bb, VOL_or_BND = BND)
gfu.Set(bb, definedon=mesh.Boundaries('rub'))


omega=314
mu0 = 1.257e-6
rel = 200
sigma = 2e3
sig = mesh.MaterialCF({ "inner" : 1 }, default=None)
#sigma=CoefficientFunction( (10, 0,  0, 10), dims=(2,2) )

#rel=CoefficientFunction( ( 1200, 0, 0, 1200), dims=(2,2) )

a = BilinearForm(fes)
a += (1/mu0)*curl(u)*curl(v)*dx('outer') + rel*curl(u)*curl(v)*dx('inner') + \
    +1j*omega*sigma*u*v*dx('inner') + 1j*1e-3*u*v*dx('outer') 

a.Assemble()


sila=CF((0,0))
f = LinearForm(fes)
f += sila*v*dx
f.Assemble()

#solvers.BVP(bf=a, lf=f, gf=gfu, pre=None, maxsteps=200, print=True)

r = - a.mat * gfu.vec
gfu.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs())*r

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


testu = GridFunction(fes)
testu.vec[:] = 0
k=43
testu.vec[k]=gfu.vec[k]
print('k=',k,'testu.vec[k]=',testu.vec[k])
Draw(testu)



""" import numpy as np
dirich=GridFunction(fes)
dirich.Set(bb, definedon=mesh.Boundaries('rub'))

rows,cols,vals = a.mat.COO()
import scipy.sparse as sp
AA = sp.csr_matrix((vals,(rows,cols)))

S=AA.todense() #/795544.94828958

maska=fes.FreeDofs()
print(maska)
for dof in range(fes.ndof):
    if (not maska[dof]):
        S[dof,:]=0
        S[dof,dof]=1

print(S)
print(dirich.vec)
print(gfu.vec)

#determ=np.linalg.det(S)
Seigenvalues=np.linalg.eigvals(S)
print('eigenvaules = ',Seigenvalues)
 """
