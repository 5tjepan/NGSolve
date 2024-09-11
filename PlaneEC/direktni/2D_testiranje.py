from ngsolve import *
from netgen.occ import *

#outer= Circle((1,1), 2).Face()
#outer.edges.name = 'rub'

outer = MoveTo(0, -1).Rectangle(3,2).Face()
outer.edges.name="rub"
outer.edges.Max(X).name = "r"
outer.edges.Min(X).name = "l"
outer.edges.Min(Y).name = "b"
outer.edges.Max(Y).name = "t"
#outer.edges.Min(X).maxh=0.1

inner = MoveTo(1, 0).Rectangle(1, 1).Face()
#inner = MoveTo(1.5,0.5).Line(1,1).Line(-1,1).Line(-1,-1).Close().Face()
inner.edges.name="interface"
inner.faces.maxh=1
outer = outer - inner

inner.faces.name="inner"
inner.faces.col = (1, 1, 0)  #colour
outer.faces.name="outer"

geo = Glue([inner, outer])
#Draw(OCCGeometry(geo));

mesh = Mesh(OCCGeometry(geo, dim=2).GenerateMesh(maxh=0.9, quad_dominated=True))
fes = HCurl(mesh, order=0, dirichlet="r|l|b|t",  complex=True, nograds = False) #CMPLX
gfu = GridFunction(fes)
print('mesh.nv',mesh.nv)
print('mesh.nedge', mesh.nedge)
print('mesh.ne', mesh.ne)

u, v = fes.TnT()

#bb=CF((10*y,-10*x))
bb=CF((10,10))

#gfu.Set(bb, VOL_or_BND = BND)
#gfu.Set(bb, definedon=mesh.Boundaries('rub'))
gfu.Set(bb, definedon=mesh.Boundaries('r'))

Draw(gfu)

omega=314
mu0 = 1.257e-6
rel = 1200
#sigma = 2e3
sig = mesh.MaterialCF({ "inner" : 1 }, default=None)
sigma=CoefficientFunction( (10, 0,  0, 10), dims=(2,2) )

#rel=CoefficientFunction( ( 1200, 0, 0, 1200), dims=(2,2) )

a = BilinearForm(fes)
a += (1/mu0)*curl(u)*curl(v)*dx #('outer') + rel*curl(u)*curl(v)*dx('inner') #+ \
    #+ 1j*1e-1*u*v*dx('inner') #+1j*omega*sigma*u*v*dx('inner')

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



""" rows,cols,vals = a.mat.COO()
import scipy.sparse as sp
AA = sp.csr_matrix((vals,(rows,cols)))
import numpy as np
AD=AA.todense()
#print(AD)
invAD=np.linalg.inv(AD)
print(invAD)
 """
