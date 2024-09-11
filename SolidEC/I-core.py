from netgen.csg import *
from ngsolve import *






def MakeGeometry():
    geometry = CSGeometry()
    box = OrthoBrick(Pnt(-1,-1,-1),Pnt(1,1,1)).bc("outer")
  
    front= Plane (Pnt(-0.05,-0.05,-0.2), Vec(-1,0,0) ).bc("front").maxh(0.03)
    right= Plane (Pnt(-0.05,-0.05,-0.2), Vec(0,-1,0) ).bc("right")
    bot  = Plane (Pnt(-0.05,-0.05,-0.2), Vec(0,0,-1) ).bc("bot")
    back = Plane (Pnt(0.05,0.05,0.2), Vec(1, 0,0) ).bc("back").maxh(0.03)
    left = Plane (Pnt(0.05,0.05,0.2), Vec(0,1,0) ).bc("left")
    top  = Plane (Pnt(0.05,0.05,0.2), Vec(0,0, 1) ).bc("top")
    core = left * right * front * back * bot * top
    #core = OrthoBrick(Pnt(-0.05,-0.05,-0.2),Pnt(0.05,0.05,0.2))
    core.maxh(0.03)
    #front.maxh(0.01)
    core.mat("core")

    coil = (Cylinder(Pnt(0.0,0,-0.5), Pnt(0.0,0,0.5), 0.2) - \
            Cylinder(Pnt(0.0,0,-0.5), Pnt(0.0,0,0.5), 0.1)) * \
            OrthoBrick (Pnt(-1,-1,-0.1),Pnt(1,1,0.1))
    coil.maxh(0.1)
    coil.mat("coil")
    
    geometry.Add ((box-core-coil).mat("air"))
    geometry.Add (core)
    geometry.Add (coil)

    return geometry


ngmesh = MakeGeometry().GenerateMesh(maxh=0.35)
#ngmesh.Save("coil.vol")
mesh = Mesh(ngmesh)
# curve elements for geometry approximation
mesh.Curve(5)  #broj 5 je stupanj krivulje? 

Draw(mesh)

fes = HCurl(mesh, order=1, dirichlet="outer", complex=True, nograds = False)
u = fes.TrialFunction()
v = fes.TestFunction()

mur = mesh.MaterialCF({ "core" : 1000 }, default=1)
mu0 = 1.257e-6
nu = 1/(mu0*mur)

omega=314

sig = mesh.MaterialCF({ "core" : 1 }, default=None)
sigma=CoefficientFunction( (0.01 , 0, 0,   0, 2e3, 0,  0, 0, 2e3), dims=(3,3) )
rel=CoefficientFunction( ( 32000, 0, 0,   0, 1200, 0,  0, 0, 1200), dims=(3,3) )

rho=CF( (100 , 0, 0,   0, 0.0005, 0,  0, 0, 0.0005), dims=(3,3) )

#a = BilinearForm(fes, symmetric=True)

term1 = (1/mu0)*curl(u)*curl(v)*dx('air|coil') + rel*curl(u)*curl(v)*dx('core') + \
    1j*omega*sigma*u*v*dx('core')
term2 = 0.01*1e0*u*v*dx('air|coil')

a = BilinearForm(term1+term2)

#c = Preconditioner(a, type="bddcc")
#c = Preconditioner(a, type="multigrid")

#Js=CF((y,-x,0))
#Draw(Js, mesh, 'Js')

f = LinearForm(fes)

Ts_air= IfPos(z*z-0.01, (0,0,0), IfPos((x*x+y*y-0.01),CF((0,0,0.015)),(0,0,0)))
Ts_coil= CF((0,0,0.5*(x*x+y*y)-0.005))
f += Ts_coil *1e7* curl(v) * dx("coil") + Ts_air *1e7* curl(v) * dx("air|core") 

#f += CoefficientFunction((y,-x,0)) *1e7* v * dx("coil")

sol = GridFunction(fes)

a.Assemble()
f.Assemble()

#solvers.BVP(bf=a, lf=f, gf=sol, pre=c, maxsteps=200, print=True)

#needsassembling=False
#solver = CGSolver(mat=a.mat, pre=None, maxsteps=50000)
#sol.vec.data = solver * f.vec

sol.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec

J = - 1j * omega*sig * sigma * sol
B = curl(sol)

Pow=0.5*rho*J*Conj(J) 
Peddy=Integrate(Pow, mesh, order=5)
print(Peddy) 


Draw (B, mesh, "B", draw_surf=False)
Draw (J, mesh, "J", draw_surf=False)
Draw (sol, mesh, "A", draw_surf=True)



""" B=curl(sol)
Bz=B[1]
feh=H1(mesh, order=1)
H1B=GridFunction(feh)
H1B.Set(Bz)
Draw(H1B) """
