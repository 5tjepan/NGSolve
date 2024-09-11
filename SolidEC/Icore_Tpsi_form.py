from netgen.csg import *
from ngsolve import *


def MakeGeometry():
    geometry = CSGeometry()
    box = OrthoBrick(Pnt(-1,-1,-1),Pnt(1,1,1)).bc("outer")

    core = OrthoBrick(Pnt(-0.05,-0.05,-0.2),Pnt(0.05,0.05,0.2)).bc("interface")
    core.maxh(0.02)
    core.mat("core")
    
    coil = (Cylinder(Pnt(0.0,0,-0.5), Pnt(0.0,0,0.5), 0.2) - \
            Cylinder(Pnt(0.0,0,-0.5), Pnt(0.0,0,0.5), 0.1)) * \
            OrthoBrick (Pnt(-1,-1,-0.1),Pnt(1,1,0.1))
    coil.maxh(0.04)
    coil.mat("coil")
    
    geometry.Add ((box-core-coil).mat("air"))
    geometry.Add (core)
    geometry.Add (coil)

    return geometry


ngmesh = MakeGeometry().GenerateMesh(maxh=0.3)
#ngmesh.Save("coil.vol")
mesh = Mesh(ngmesh)
# curve elements for geometry approximation
mesh.Curve(5)  #broj 5 je stupanj krivulje? 

Draw(mesh)

fsU = HCurl(mesh, order=0, definedon='core', dirichlet="interface", complex=True, nograds = False)
fsV = H1(mesh, order=1, complex=True)
fes=fsU*fsV
cvp, msp = fes.TrialFunction()
alpha, phi = fes.TestFunction()

mur = mesh.MaterialCF({ "core" : 1000 }, default=1)
mu0 = 1.257e-6
mu= mu0*mur
nu = 1/(mu)
rho= 0.001
omega=314

yes_coil = mesh.MaterialCF({ "core|air" : 0 }, default=1)
not_coil = mesh.MaterialCF({ "coil" : 0 }, default=1)
#sigma=CoefficientFunction( (1 , 0, 0,   0, 2e3, 0,  0, 0, 2e3), dims=(3,3) )
#rho = CF((1 , 0, 0,   0, 0.5e-3, 0,  0, 0, 0.5e-3), dims=(3,3))
#rel=CoefficientFunction( ( 32000, 0, 0,   0, 1200, 0,  0, 0, 1200), dims=(3,3) )

#a = BilinearForm(fes, symmetric=True)
term1 = -rho*curl(cvp)*curl(alpha)*dx('core') + 1j*omega*mu*cvp*alpha*dx('core') \
    + 1j*omega*mu*grad(msp)*alpha*dx('core')
term2 = -mu*cvp*grad(phi)*dx('core') -mu*grad(msp)*grad(phi)*dx
term3 = 0.1*msp*phi*dx #regularization
a = BilinearForm(term1+term2+term3)

Ts_air= 10e6*IfPos(z*z-0.01, (0,0,0), IfPos((x*x+y*y-0.01),CF((0,0,0.015)),(0,0,0)))
Ts_coil= 10e6*CF((0,0,0.5*(x*x+y*y)-0.005))
force1 =-1j*omega*mu*Ts_coil*alpha*dx("coil") - 1j*omega*mu*Ts_air*alpha*dx("air|core") 
force2 =mu*Ts_coil*grad(phi)*dx("coil") +mu*Ts_air*grad(phi)*dx("air|core")
f = LinearForm(force1 + force2)
Ts= yes_coil*Ts_coil + not_coil*Ts_air

a.Assemble()
f.Assemble()

#c = Preconditioner(a, type="bddcc")
#c = Preconditioner(a, type="multigrid")
#c= Preconditioner(a, type='direct', inverse='umfpack')

sol = GridFunction(fes)
Tpot, Ppot = sol.components

#solvers.BVP(bf=a, lf=f, gf=sol, pre=c, maxsteps=200, print=True)
#solver = CGSolver(mat=a.mat, pre=None, maxsteps=50000)
#sol.vec.data = solver * f.vec

sol.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec

J = curl(Tpot)
B = mu * Tpot + mu*grad(Ppot)# + mu*Ts 

Draw (Ts, mesh, "Ts")
""" Draw (Ts_coil, mesh, "Tcoil")
Draw (Ts_air, mesh, "Tair")
Draw(yes_coil,mesh,"yescoil")
Draw(not_coil,mesh,"notcoil") """

Draw (J, mesh, "J")
Draw (Ppot, mesh, "psi")
Draw (Tpot, mesh, "T")
Draw (B, mesh, "B")




""" import numpy as np

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
print('eigenvaules = ',Seigenvalues) """

