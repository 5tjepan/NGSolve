from netgen.csg import *
from ngsolve import *

csize=0.0105
def MakeGeometry(csize):
    geometry = CSGeometry()
    w=0.8
    box = OrthoBrick(Pnt(-w,-w,-w),Pnt(w,w,w)).bc("outer")
    h=0.1
    front= Plane (Pnt(-0.03,-0.05,-h), Vec(-1,0,0) ).bc("front").maxh(csize)
    right= Plane (Pnt(-0.03,-0.05,-h), Vec(0,-1,0) ).bc("right")
    bot  = Plane (Pnt(-0.03,-0.05,-h), Vec(0,0,-1) ).bc("bot")
    back = Plane (Pnt(0.03,0.05,h), Vec(1, 0,0) ).bc("back").maxh(csize)
    left = Plane (Pnt(0.03,0.05,h), Vec(0,1,0) ).bc("left")
    top  = Plane (Pnt(0.03,0.05,h), Vec(0,0, 1) ).bc("top")
    core = left * right * front * back * bot * top
    #core = OrthoBrick(Pnt(-0.05,-0.05,-0.2),Pnt(0.05,0.05,0.2))
    core.maxh(0.0105)
    #front.maxh(0.01)
    core.mat("core")
    
    coil = (Cylinder(Pnt(0.0,0,-0.5), Pnt(0.0,0,0.5), 0.2) - \
            Cylinder(Pnt(0.0,0,-0.5), Pnt(0.0,0,0.5), 0.1)) * \
            OrthoBrick (Pnt(-1,-1,-0.075),Pnt(1,1,0.075))
    coil.maxh(0.09)
    coil.mat("coil")
    
    geometry.Add ((box-core-coil).mat("air"))
    geometry.Add (core)
    geometry.Add (coil)

    return geometry


ngmesh = MakeGeometry(csize).GenerateMesh(maxh=0.27)
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


rel_yz=1200
rel_x=32000
rho_yz=5e-7
omega=314
rot=CF( (0 , 0, 0,   0, 0, 1,  0, -1, 0), dims=(3,3) )

#sig = mesh.MaterialCF({ "core" : 1 }, default=None)
#sigma=CoefficientFunction( (1 , 0, 0,   0, 2e3, 0,  0, 0, 2e3), dims=(3,3) )
rho=CF( (0.1 , 0, 0,   0, rho_yz, 0,  0, 0, rho_yz), dims=(3,3) )
rel=CF( ( rel_x, 0, 0,   0, rel_yz, 0,  0, 0, rel_yz), dims=(3,3) )


d=0.0003
kh= 0.1/11/9 #csize/2
khx=kh #0.06/7/5
xdir=CF((1,0,0))

#term1=(1/mu0)*curl(mvp)*curl(alpha)*dx('air|coil') + rel*curl(mvp)*curl(alpha)*dx('core') + \
#1j*omega*(rot*grad(csp))*alpha*dx('core') #+ 0.1*mvp*alpha*dx('air|coil')
#term2= -1j*omega*rho*(rot*grad(csp))*(rot*grad(tau))*dx('core') + 1j*omega*(rot*grad(tau))*mvp*dx('core')
#term3=1*1e0*mvp*alpha*dx

""" 
term1 = rho*grad(csp)*grad(theta)*dx + 1j*omega*curl(mvp)*theta*dx
term1hg= 1j*omega*(1)/(12*rel)*(kh**2) * grad(csp)*grad(theta)*dx
term2 = - 1j*omega* rel*curl(mvp)*curl(alpha)*dx + 1j*omega*csp*curl(alpha)*dx
term2hg= 1/12*(kh*omega)**2 /rho *curl(mvp)*curl(alpha)*dx """

term1=-1j*omega*(1/mu0)*curl(mvp)*curl(alpha)*dx('air|coil')-1j*omega*rel*curl(mvp)*curl(alpha)*dx('core') \
    +1j*omega* (rot*grad(csp))*alpha*dx('core') #+ 0.1*mvp*alpha*dx('air|coil')
term2= rho*(rot*grad(csp))*(rot*grad(tau))*dx('core') + 1j*omega*(rot*grad(tau))*mvp*dx('core')
term3=1*1e0*mvp*alpha*dx
term4 = 1j*omega*(1)/(12*rel_x)*(kh**2) * (rot*grad(csp))*(rot*grad(tau))*dx('core')
term5 = 1/12*(kh*omega)**2 /rho_yz *(curl(mvp)*xdir)*(curl(alpha)*xdir)*dx('core')
term6= -1j*omega*1/(12*rel_yz)*(khx**2) * (rot*grad(csp))*(rot*grad(tau))*dx('core') #doprinos Byz zbog Jyz
term7 = 1/12*(d*omega)**2 /rho_yz *(curl(mvp)*curl(alpha)-(curl(mvp)*xdir)*(curl(alpha)*xdir))*dx('core')

a = BilinearForm(term1+term2+term3)#+term4+term5 + term6)
#a = BilinearForm(term1+term2+term3+term5)
#a += (term1+term2)

#c = Preconditioner(a, type="bddcc")
#c = Preconditioner(a, type="multigrid")
#c= Preconditioner(a, type='direct', inverse='umfpack')

#Js=CF((y,-x,0))
#Draw(Js, mesh, 'Js')

f = LinearForm(fes)

Ts_air= IfPos(z*z-(0.075)**2, (0,0,0), IfPos((x*x+y*y-0.01),CF((0,0,0.015)),(0,0,0)))
Ts_coil= CF((0,0,0.5*(x*x+y*y)-0.005))
f += Ts_coil *1e7* curl(alpha) * dx("coil") + Ts_air *1e7* curl(alpha) * dx("air|core") 

#f += CoefficientFunction((y,-x,0)) *1e7* v * dx("coil")

#import time
#start = time.time()

a.Assemble()
f.Assemble()

sol = GridFunction(fes)
Apot, Tpot = sol.components

#Tpot.Set(CF(0),BND)

#solvers.BVP(bf=a, lf=f, gf=sol, pre=None, maxsteps=2, print=True, needsassembling=False)

#end = time.time()
#print(end-start)

#solver = CGSolver(mat=a.mat, pre=None, maxsteps=50000)
#sol.vec.data = solver * f.vec
#solver.Inverse(freedofs=fes.FreeDofs())

sol.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec

Tx= Tpot   #- 1j * omega*Tpot
J = rot*grad(Tpot) # - 1j * omega*(rot*grad(Tpot))
B=curl(Apot)
Pow=0.5*rho*J*Conj(J)# + 0.5* 1/12*(kh*omega)**2 /rho_yz *(B*xdir)*(Conj(B)*xdir)
Peddy=Integrate(Pow.real, mesh, order=5)
print(Peddy) 


Draw (J, mesh, "J")
Draw (Tx, mesh, "Tx")
Draw (B, mesh, "B")
#Draw(Ts_air,mesh,"Ts_air")