from ngsolve import *
from netgen.occ import *

outer= Circle((1,1), 2).Face()
outer.edges.name = 'rub'

""" outer = MoveTo(0, -1).Rectangle(3,2).Face()
outer.edges.name="rub"
outer.edges.Max(X).name = "r"
outer.edges.Min(X).name = "l"
outer.edges.Min(Y).name = "b"
outer.edges.Max(Y).name = "t"
#outer.edges.Min(X).maxh=0.1 """

inner = MoveTo(0, 0).Rectangle(2,2).Face()
inner.edges.name="rub"
inner.edges.Max(X).name = "r"
inner.edges.Min(X).name = "l"
inner.edges.Min(Y).name = "b"
inner.edges.Max(Y).name = "t"
#outer.edges.Min(X).maxh=0.1

#inner = MoveTo(0, 0).Rectangle(2,2).Face()
#inner = MoveTo(1.5,0.5).Line(1,1).Line(-1,1).Line(-1,-1).Close().Face()
#inner.edges.name="interface"

inner.faces.maxh=1
outer = outer - inner

inner.faces.name="inner"
inner.faces.col = (1, 1, 0)  #colour
outer.faces.name="outer"

geo = Glue([inner, outer])
#Draw(OCCGeometry(geo));

mesh = Mesh(OCCGeometry(geo, dim=2).GenerateMesh(maxh=0.2, quad_dominated=False))

fsU = HCurl(mesh, order=0, dirichlet="rub", complex=True, nograds = False)
fsV = H1(mesh, order=1, definedon="inner", complex=True)
#fsV = H1(mesh, order=1, dirichlet="r", definedon="inner", complex=True)
fes=fsU*fsV
mvp, esp = fes.TrialFunction()
alpha, phi = fes.TestFunction()

print('mesh.nv',mesh.nv)
print('mesh.nedge', mesh.nedge)
print('mesh.ne', mesh.ne)

u, v = fes.TnT()

bb=CF((10*y,-10*x))
#bb=CF((10,10))

#gfu.Set(bb, VOL_or_BND = BND)

sol = GridFunction(fes)
Apot, Vpot = sol.components
Apot.Set(bb, definedon=mesh.Boundaries('rub'))

#Draw(gfu)

omega=314
mu0 = 1.257e-6
rel = 1200
#sigma = 2e3
sig = mesh.MaterialCF({ "inner" : 1 }, default=None)
sigma=CoefficientFunction( (10, 0,  0, 10), dims=(2,2) )

#rel=CoefficientFunction( ( 1200, 0, 0, 1200), dims=(2,2) )

term1=(1/mu0)*curl(mvp)*curl(alpha)*dx('outer') + rel*curl(mvp)*curl(alpha)*dx('inner') + \
    1j*omega*sigma*mvp*alpha*dx('inner') + 1j*omega*sigma*grad(esp)*alpha*dx('inner')
term2=1j*omega*sigma*mvp*grad(phi)*dx('inner') + 1j*omega*sigma*grad(esp)*grad(phi)*dx('inner')
term3=1j*1e0*mvp*alpha*dx('outer')+1j*1e0*(esp)*(phi)*dx('inner')

term4=1j*1e0*mvp*alpha*dx('outer') + 1j*1e0*grad(esp)*alpha*dx('outer') + \
    1j*1e0*mvp*grad(phi)*dx('outer') + 1j*1e0*grad(esp)*grad(phi)*dx('outer')

a = BilinearForm(term1+term2+term3)

a.Assemble()
r = - a.mat * sol.vec

#gfu.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs())*r


sol.vec.data += a.mat.Inverse(fes.FreeDofs()) * r

J = - 1j * omega*sig * sigma * Apot - 1j*omega*sig*sigma*grad(Vpot)
B=curl(Apot)

Draw(Vpot, mesh, "V")
Draw(Apot, mesh, "A")
Draw (B, mesh, "B")
Draw (J, mesh, "J")

#print(gfu.vec)

#for i in range(fes.ndof):
#    print(fes.CouplingType(i))

import numpy as np
dirich=GridFunction(fes)

Apotdir, Vpotdir = dirich.components
Apotdir.Set(bb, definedon=mesh.Boundaries('rub'))

rows,cols,vals = a.mat.COO()
import scipy.sparse as sp
AA = sp.csr_matrix((vals,(rows,cols)))

S=AA.todense() #/795544.94828958
print('shape of S', np.shape(S))

print('fsU.ndof',fsU.ndof)
print('fsV.ndof',fsV.ndof)
print('fes.ndof',fes.ndof)

maskaU=fsU.FreeDofs()
print('maskaU',maskaU)

maskaV=fsV.FreeDofs()
print('maskaV',maskaV)

maska=fes.FreeDofs()
print(maska)
for dof in range(max(np.shape(S))):
    if (not maska[dof]):
        S[dof,:]=0
        S[dof,dof]=1

#print(S)
#print(dirich.vec)
#print(sol.vec)

#determ=np.linalg.det(S)
Seigenvalues=np.linalg.eigvals(S)
print('eigenvaules = ',Seigenvalues)


