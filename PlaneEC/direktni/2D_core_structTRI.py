from ngsolve import *
from netgen.occ import *

outer= Circle((0,0), 3).Face()
outer.edges.name = 'rub'

""" outer = MoveTo(0, -1).Rectangle(3,2).Face()
outer.edges.name="rub"
outer.edges.Max(X).name = "r"
outer.edges.Min(X).name = "l"
outer.edges.Min(Y).name = "b"
outer.edges.Max(Y).name = "t"
#outer.edges.Min(X).maxh=0.1 """


""" inner1 = MoveTo(0.0,-1.0).Line(1,1).Line(-0.5,0.5).Line(-1,-1).Close().Face()
inner1.edges.name="interface1"
inner1.faces.maxh=0.05

inner2 = MoveTo(0.1,-1.0).Line(1,1).Line(0.5,-0.5).Line(-1,-1).Close().Face()
#inner = MoveTo(1,0).Line(1,1).Line(-1,1).Line(-1,-1).Close().Face()
inner2.edges.name="interface2"
inner2.faces.maxh=0.05 """

w=1
inner1 = MoveTo(-w,-1.5).Rectangle(w, 3).Face()
#inner = MoveTo(1,0).Line(1,1).Line(-1,1).Line(-1,-1).Close().Face()
inner1.edges.name="interface1"
inner1.faces.maxh=1

inner2 = MoveTo(0.0, -1.5).Rectangle(w, 3).Face()
#inner = MoveTo(1,0).Line(1,1).Line(-1,1).Line(-1,-1).Close().Face()
inner2.edges.name="interface2"
inner2.faces.maxh=1


outer = outer - inner1 - inner2

inner1.faces.name="inner1"
inner1.faces.col = (1, 1, 0)  #colour

inner2.faces.name="inner2"
inner2.faces.col = (0.5, 1, 1)  #colour
outer.faces.name="outer"


geo = Glue([inner1,inner2, outer])
#Draw(OCCGeometry(geo));


mesh = Mesh(OCCGeometry(geo, dim=2).GenerateMesh(maxh=1, quad_dominated=False))
fes = HCurl(mesh, order=0, dirichlet="rub",  complex=True, nograds = False) #CMPLX
gfu = GridFunction(fes)
print('mesh.nv',mesh.nv)
print('mesh.nedge', mesh.nedge)
print('mesh.ne', mesh.ne)

u, v = fes.TnT()

bb=CF((1*y,-1*x))
#bb=CF((10,10))

#gfu.Set(bb, VOL_or_BND = BND)
gfu.Set(bb, definedon=mesh.Boundaries('rub'))

Draw(gfu)

omega=314
mu0 = 1.257e-6
rel = 1200
#sigma = 2e3
sig = mesh.MaterialCF({ "inner1|inner2" : 1 }, default=None)
sigma=CoefficientFunction( (0, 0,  0, 100), dims=(2,2) )

#rel=CoefficientFunction( ( 1200, 0, 0, 1200), dims=(2,2) )

a = BilinearForm(fes)
a += (1/mu0)*curl(u)*curl(v)*dx('outer') + rel*curl(u)*curl(v)*dx('inner1|inner2') + \
    + 1j*1e-3*u*v*dx('outer') +1j*omega*sigma*u*v*dx('inner1|inner2')

a.Assemble()
r = - a.mat * gfu.vec
gfu.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs())*r


J = - 1j * omega * sigma * gfu *sig
B = curl(gfu)
A=gfu

Draw(A, mesh, "A")
Draw (B, mesh, "B")
Draw (J, mesh, "J")

#print(gfu.vec)



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
