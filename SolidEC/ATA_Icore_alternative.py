from netgen.csg import *
from ngsolve import *


def MakeGeometry():
    geometry = CSGeometry()
    box = OrthoBrick(Pnt(-1,-1,-1),Pnt(1,1,1)).bc("outer")

    front= Plane (Pnt(-0.05,-0.05,-0.2), Vec(-1,0,0) ).bc("front").maxh(0.015)
    right= Plane (Pnt(-0.05,-0.05,-0.2), Vec(0,-1,0) ).bc("right")
    bot  = Plane (Pnt(-0.05,-0.05,-0.2), Vec(0,0,-1) ).bc("bot")
    back = Plane (Pnt(0.05,0.05,0.2), Vec(1, 0,0) ).bc("back").maxh(0.015)
    left = Plane (Pnt(0.05,0.05,0.2), Vec(0,1,0) ).bc("left")
    top  = Plane (Pnt(0.05,0.05,0.2), Vec(0,0, 1) ).bc("top")
    core = left * right * front * back * bot * top
    #core = OrthoBrick(Pnt(-0.05,-0.05,-0.2),Pnt(0.05,0.05,0.2))
    core.maxh(0.06)
    #front.maxh(0.01)
    core.mat("core")
    
    coil = (Cylinder(Pnt(0.0,0,-0.5), Pnt(0.0,0,0.5), 0.2) - \
            Cylinder(Pnt(0.0,0,-0.5), Pnt(0.0,0,0.5), 0.1)) * \
            OrthoBrick (Pnt(-1,-1,-0.1),Pnt(1,1,0.1))
    coil.maxh(0.12)
    coil.mat("coil")
    
    geometry.Add ((box-core-coil).mat("air"))
    geometry.Add (core)
    geometry.Add (coil)

    return geometry


ngmesh = MakeGeometry().GenerateMesh(maxh=0.4)
#ngmesh.Save("coil.vol")
mesh = Mesh(ngmesh)
# curve elements for geometry approximation
mesh.Curve(5)  #broj 5 je stupanj krivulje? 

Draw(mesh)
print(mesh.GetBoundaries())

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

version=4

if version==1:   
    term1=(1/mu0)*curl(mvp)*curl(alpha)*dx('air|coil') + rel*curl(mvp)*curl(alpha)*dx('core') + \
    1j*omega*(rot*grad(csp))*alpha*dx('core') #+ 0.1*mvp*alpha*dx('air|coil')
    term2= -1j*omega*rho*(rot*grad(csp))*(rot*grad(tau))*dx('core') + 1j*omega*(rot*grad(tau))*mvp*dx('core')
    term3=1*mvp*alpha*dx

elif version==2: #ALTERNATIVNA VERZIJA::::::::::::
    term1=(1/mu0)*curl(mvp)*curl(alpha)*dx('air|coil') + rel*curl(mvp)*curl(alpha)*dx('core') \
    - (rot*grad(csp))*alpha*dx('core') #+ 0.1*mvp*alpha*dx('air|coil')
    term2= -1j/omega*rho*(rot*grad(csp))*(rot*grad(tau))*dx('core') + (rot*grad(tau))*mvp*dx('core')
    term3=1*mvp*alpha*dx

elif version==3: #ALTERNATIVNA VERZIJA::::::::::::kao version 2 ali zamjena: curlT*A => T*curlA
    term1=(1/mu0)*curl(mvp)*curl(alpha)*dx('air|coil') + rel*curl(mvp)*curl(alpha)*dx('core') \
    - csp*CF((1,0,0))*curl(alpha)*dx('core') #+ 0.1*mvp*alpha*dx('air|coil')
    term2= -1j/omega*rho*(rot*grad(csp))*(rot*grad(tau))*dx('core') + tau*CF((1,0,0))*curl(mvp)*dx('core')
    term3=1*mvp*alpha*dx

elif version==4: #ALTERNATIVNA VERZIJA::::::::::::
    term1=-1j*omega*(1/mu0)*curl(mvp)*curl(alpha)*dx('air|coil') -1j*omega* rel*curl(mvp)*curl(alpha)*dx('core') \
    +1j*omega* (rot*grad(csp))*alpha*dx('core') #+ 0.1*mvp*alpha*dx('air|coil')
    term2= rho*(rot*grad(csp))*(rot*grad(tau))*dx('core') + 1j*omega*(rot*grad(tau))*mvp*dx('core')
    term3=1e0*mvp*alpha*dx

a = BilinearForm(term1+term2+term3)

#c = Preconditioner(a, type="bddcc")
#c = Preconditioner(a, type="multigrid")
c= Preconditioner(a, type='direct', inverse='umfpack')

#Js=CF((y,-x,0))
#Draw(Js, mesh, 'Js')

f = LinearForm(fes)

Ts_air= IfPos(z*z-0.01, (0,0,0), IfPos((x*x+y*y-0.01),CF((0,0,0.015)),(0,0,0)))
Ts_coil= CF((0,0,0.5*(x*x+y*y)-0.005))
if version==4:
    f += (-1j*omega*Ts_coil *1e7* curl(alpha)) * dx("coil") - 1j*omega*Ts_air *1e7* curl(alpha) * dx("air|core") 
else:
    f += Ts_coil *1e7* curl(alpha) * dx("coil") + Ts_air *1e7* curl(alpha) * dx("air|core") 

#f += CoefficientFunction((y,-x,0)) *1e7* v * dx("coil")

#import time
#start = time.time()

a.Assemble()
f.Assemble()

sol = GridFunction(fes)
Apot, Tpot = sol.components

#Tpot.Set(CF(0),BND)

solvers.BVP(bf=a, lf=f, gf=sol, pre=c.mat, maxsteps=2000, print=True, needsassembling=False)

#end = time.time()
#print(end-start)

#solver = CGSolver(mat=a.mat, pre=None, maxsteps=50000)
#sol.vec.data = solver * f.vec
#solver.Inverse(freedofs=fes.FreeDofs())

#sol.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec

if version==1:
    Tx= - 1j * omega*Tpot
    J = - 1j * omega*(rot*grad(Tpot))
else:
    Tx = Tpot
    J= rot*grad(Tpot)

B=curl(Apot)
Pow=0.5*rho*J*Conj(J) 
Peddy=Integrate(Pow, mesh, order=5)
print(Peddy) 

Draw (J, mesh, "J")
Draw (Tx, mesh, "Tx")
Draw (B, mesh, "B")

