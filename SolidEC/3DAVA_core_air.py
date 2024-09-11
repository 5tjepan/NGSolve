from netgen.csg import *
from ngsolve import *


def MakeGeometry():
    geometry = CSGeometry()
    #box = OrthoBrick(Pnt(-1,-1,-1),Pnt(1,1,1)).bc("outer")
    box = OrthoBrick(Pnt(-0.2,-0.2,-0.2),Pnt(0.2,0.2,0.2)).bc("outer")

    core = OrthoBrick(Pnt(-0.05,-0.05,-0.1),Pnt(0.05,0.05,0.1))
    core.maxh(0.99)
    core.mat("core")
    
      
    geometry.Add ((box-core).mat("air"))
    geometry.Add (core)
  

    return geometry


ngmesh = MakeGeometry().GenerateMesh(maxh=0.99)
#ngmesh.Save("coil.vol")
mesh = Mesh(ngmesh)
# curve elements for geometry approximation
mesh.Curve(5)  #broj 5 je stupanj krivulje? 

Draw(mesh)

fsU = HCurl(mesh, order=0, dirichlet="outer", complex=True, nograds = False)
fsV = H1(mesh, order=1, definedon="core", complex=True)
fes=fsU*fsV
mvp, esp = fes.TrialFunction()
alpha, phi = fes.TestFunction()

mur = mesh.MaterialCF({ "core" : 1000 }, default=1)
mu0 = 1.257e-6
nu = 1/(mu0*mur)

omega=314

sig = mesh.MaterialCF({ "core" : 1 }, default=None)
sigma=CoefficientFunction( (1 , 0, 0,   0, 2e3, 0,  0, 0, 2e3), dims=(3,3) )
rel=CoefficientFunction( ( 32000, 0, 0,   0, 1200, 0,  0, 0, 1200), dims=(3,3) )


#a = BilinearForm(fes, symmetric=True)

term1=(1/mu0)*curl(mvp)*curl(alpha)*dx('air') + rel*curl(mvp)*curl(alpha)*dx('core') + \
    1j*omega*sigma*mvp*alpha*dx('core') + 1j*omega*sigma*grad(esp)*alpha*dx('core')
term2=1j*omega*sigma*mvp*grad(phi)*dx('core') + 1j*omega*sigma*grad(esp)*grad(phi)*dx('core')
term3=1j*1e0*mvp*alpha*dx('air') + 1j*1e0*(esp)*(phi)*dx('core')

a = BilinearForm(term1+term2+term3)


#c = Preconditioner(a, type="bddcc")
#c = Preconditioner(a, type="multigrid")
#c= Preconditioner(a, type='direct', inverse='umfpack')

#Js=CF((y,-x,0))
#Draw(Js, mesh, 'Js')

f = LinearForm(fes)

""" Ts_air= IfPos(z*z-0.01, (0,0,0), IfPos((x*x+y*y-0.01),CF((0,0,0.015)),(0,0,0)))
Ts_coil= CF((0,0,0.5*(x*x+y*y)-0.005))
f += Ts_coil *1e7* curl(alpha) * dx("coil") + Ts_air *1e7* curl(alpha) * dx("air|core") """ 

#f += CoefficientFunction((y,-x,0)) *1e7* alpha * dx("coil")

a.Assemble()
f.Assemble()

sol = GridFunction(fes)
Apot, Vpot = sol.components

bb=CF((10*y,-10*x, 0))
Apot.Set(bb, definedon=mesh.Boundaries('outer'))

#solvers.BVP(bf=a, lf=f, gf=sol, pre=None, maxsteps=200, print=True)
#solver = CGSolver(mat=a.mat, pre=None, maxsteps=50000)
#sol.vec.data = solver * f.vec

#sol.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec

J = - 1j * omega*sig * sigma * Apot - 1j*omega*sig*sigma*grad(Vpot)
B=curl(Apot)

Draw (Apot, mesh, "A")
Draw (B, mesh, "B")
Draw (J, mesh, "J")



import numpy as np

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
#print(sol.vec)

#determ=np.linalg.det(S)
Seigenvalues=np.linalg.eigvals(S)
print('eigenvaules = ',Seigenvalues)

