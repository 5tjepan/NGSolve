from netgen.csg import *
from ngsolve import *


def MakeGeometry():
    geometry = CSGeometry()
    #box = OrthoBrick(Pnt(-1,-1,-1),Pnt(1,1,1)).bc("outer")
    box = OrthoBrick(Pnt(-0.2,-0.3,-0.3),Pnt(0.2,0.3,0.3)).bc("outer")

    #core = OrthoBrick(Pnt(-0.05,-0.05,-0.1),Pnt(0.05,0.05,0.1))
    
    front= Plane (Pnt(-0.1,-0.1,-0.1), Vec(-1,0,0) ).bc("front")
    right= Plane (Pnt(-0.1,-0.2,-0.15), Vec(0,-1,0) ).bc("right")
    bot  = Plane (Pnt(-0.1,-0.1,-0.2), Vec(0,0,-1) ).bc("bot")
    back = Plane (Pnt(0.1,0.1,0.1), Vec(1, 0,0) ).bc("back")
    left = Plane (Pnt(0.1,0.2,0.15), Vec(0,1,0) ).bc("left")
    top  = Plane (Pnt(0.1,0.1,0.2), Vec(0,0, 1) ).bc("top")
    core = left * right * front * back * bot * top
        
    core.maxh(0.1)
    core.mat("core")
          
    geometry.Add ((box-core).mat("air"))
    geometry.Add (core)
  

    return geometry


ngmesh = MakeGeometry().GenerateMesh(maxh=0.9)
#ngmesh.Save("coil.vol")
mesh = Mesh(ngmesh)
# curve elements for geometry approximation
mesh.Curve(5)  #broj 5 je stupanj krivulje? 

Draw(mesh)

fsU = HCurl(mesh, order=0, dirichlet="outer", complex=True, nograds = False)
fsV = H1(mesh, order=1, dirichlet='left|right|bot|top', definedon='core', complex=True)
fes=fsU*fsV

mvp, csp = fes.TrialFunction()
alpha, tau = fes.TestFunction()

mur = mesh.MaterialCF({ "core" : 1000 }, default=1)
mu0 = 1.257e-6
nu = 1/(mu0*mur)

omega=314
rot=CF( (0 , 0, 0,   0, 0, 1,  0, -1, 0), dims=(3,3) )

#sig = mesh.MaterialCF({ "core" : 1 }, default=None)
sigma=CoefficientFunction( (1 , 0, 0,   0, 2e3, 0,  0, 0, 2e3), dims=(3,3) )
rho=CF( (0.1 , 0, 0,   0, 0.0005, 0,  0, 0, 0.0005), dims=(3,3) )
rel=CF( ( 32000, 0, 0,   0, 1200, 0,  0, 0, 1200), dims=(3,3) )


#a = BilinearForm(fes)

term1=(1/mu0)*curl(mvp)*curl(alpha)*dx('air') + rel*curl(mvp)*curl(alpha)*dx('core') + \
1j*omega*(rot*grad(csp))*alpha*dx('core') #+ 0.1*mvp*alpha*dx('air|coil')
term2= -1j*omega*rho*(rot*grad(csp))*(rot*grad(tau))*dx('core') + 1j*omega*(rot*grad(tau))*mvp*dx('core')

term3=1*1e0*mvp*alpha*dx #('air') + 1*1e0*(csp)*(tau)*dx('core')

a = BilinearForm(term1+term2+term3)


f = LinearForm(fes)

""" Ts_air= IfPos(z*z-0.01, (0,0,0), IfPos((x*x+y*y-0.01),CF((0,0,0.015)),(0,0,0)))
Ts_coil= CF((0,0,0.5*(x*x+y*y)-0.005))
f += Ts_coil *1e7* curl(alpha) * dx("coil") + Ts_air *1e7* curl(alpha) * dx("air|core") """ 

#f += CoefficientFunction((y,-x,0)) *1e7* alpha * dx("coil")

a.Assemble()
f.Assemble()

sol = GridFunction(fes)
Apot, Tpot = sol.components

bb=CF((10*y,-10*x, 0))
Apot.Set(bb, definedon=mesh.Boundaries('outer'))

#c = Preconditioner(a, type="bddcc")
#c = Preconditioner(a, type="multigrid")
#c= Preconditioner(a, type='direct', inverse='umfpack')

#solvers.BVP(bf=a, lf=f, gf=sol, pre=c.mat, maxsteps=200, print=True)
#solver = CGSolver(mat=a.mat, pre=None, maxsteps=50000)
#sol.vec.data = solver * f.vec

#r=a.mat*sol.vec
#sol.vec.data += a.mat.Inverse(fes.FreeDofs()) * r

Tx= - 1j * omega*Tpot
J = - 1j * omega*(rot*grad(Tpot))
B=curl(Apot)
Pow=0.5*rho*J*Conj(J) 
Peddy=Integrate(Pow, mesh, order=5)
print(Peddy) 

Draw (Apot, mesh, "A")
Draw (J, mesh, "J")
Draw (Tx, mesh, "Tx")
Draw (B, mesh, "B")



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
print('maska', maska)
for dof in range(max(np.shape(S))):
    if (not maska[dof]):
        S[dof,:]=0
        S[dof,dof]=1

#print(S)
#print(sol.vec)

print(np.linalg.inv(S))
#determ=np.linalg.det(S)
Seigenvalues=np.linalg.eigvals(S)
print('eigenvaules = ',Seigenvalues)

