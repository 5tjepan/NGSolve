from ngsolve import *
from netgen.occ import *

import sys
sys.argv = ["asdf"]

#outer= Circle((1,1), 2).Face()
#outer.edges.name = 'rub'

outer = MoveTo(0, -1).Rectangle(2,2).Face()
#outer.edges.name="rub"
outer.edges.Max(X).name = "r"
outer.edges.Min(X).name = "l"
outer.edges.Min(Y).name = "b"
outer.edges.Max(Y).name = "t"
#outer.edges.Min(X).maxh=0.1

outer.faces.name="outer"

#Draw(OCCGeometry(geo));

mesh = Mesh(OCCGeometry(outer, dim=2).GenerateMesh(maxh=1, quad_dominated=True))
fes = HCurl(mesh, order=0, dirichlet="r|l|b|t",  complex=False, nograds = False) #CMPLX
gfu = GridFunction(fes)
print('mesh.nv',mesh.nv)
print('mesh.nedge', mesh.nedge)
print('mesh.ne', mesh.ne)

u, v = fes.TnT()

import numpy as np

#bb=CF((10*y,-10*x))
bb=CF((10,10))

#gfu.Set(bb, VOL_or_BND = BND)
#gfu.Set(bb, definedon=mesh.Boundaries('rub'))

gfu.Set(bb, definedon=mesh.Boundaries('r'))
dirich=GridFunction(fes)
dirich.Set(bb, definedon=mesh.Boundaries('r'))
#dirich.Set(bb,definedon=mesh.Boundaries('r'))

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


#J = - 1j * omega * sigma * gfu
B = curl(gfu)
A=gfu

Draw(A, mesh, "A")
Draw (B, mesh, "B")
#Draw (J, mesh, "J")

#print(gfu.vec)

# it has only 7 free dofs but only six lineraly independent curl(alpha) as test functions
print("free dofs: \n", fes.FreeDofs())

#print(a.mat)

rows,cols,vals = a.mat.COO()
import scipy.sparse as sp
AA = sp.csr_matrix((vals,(rows,cols)))

S=AA.todense()[-4:,-4:] #/795544.94828958

maska=fes.FreeDofs()

for dof in range(fes.ndof):
    if (not maska[dof]):
        S[dof,:]=0
        S[dof,dof]=1

print(S)
print(dirich.vec)
print(gfu.vec)

determ=np.linalg.det(S)
Seigenvalues=np.linalg.eigvals(S)
print('eigenvaules = ',Seigenvalues)
print(np.linalg.inv(S))
print(determ)

print(S)

import matplotlib.pyplot as plt
plt.spy(S)
plt.show()

""" qq=np.array([[2,-1,0,-1],[-1,2,-1,0],[0,-1,2,-1],[-1,0,-1,2]])
eigenv=np.linalg.eigvals(qq)
print('eigenvaules = ',eigenv) """

#pp=np.linalg.inv(qq)
#print('det', np.linalg.det(qq))


""" I=np.identity(10)
I[7:9,:]=AD[7:9,:]

brid=NodeId(EDGE,2)
print(brid)
print(mesh[brid])
print(fes.GetDofNrs(brid))"""

