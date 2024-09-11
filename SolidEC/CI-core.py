from netgen.csg import *
from ngsolve import *


def MakeGeometry():
    geometry = CSGeometry()
    box = OrthoBrick(Pnt(-1,-1,-1),Pnt(2,1,2)).bc("outer")

    core = OrthoBrick(Pnt(0,-0.05,0),Pnt(0.8,0.05,1))- \
           OrthoBrick(Pnt(0.1,-1,0.1),Pnt(0.7,1,0.9))
    core.maxh(0.05)
    core.mat("core")

    coil = (Cylinder(Pnt(0.05,0,0), Pnt(0.05,0,1), 0.3) - \
            Cylinder(Pnt(0.05,0,0), Pnt(0.05,0,1), 0.15)) * \
            OrthoBrick (Pnt(-1,-1,0.3),Pnt(1,1,0.7))
    coil.maxh(0.2)
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

fes = HCurl(mesh, order=1, dirichlet="outer",  complex=True, nograds = False)
u = fes.TrialFunction()
v = fes.TestFunction()

mur = mesh.MaterialCF({ "core" : 1000 }, default=1)
mu0 = 1.257e-6
nu = 1/(mu0*mur)

omega=314

sig = mesh.MaterialCF({ "core" : 1 }, default=0e-6)
sigma=CoefficientFunction( ( 2e3, 0, 0,   0, 0, 0,  0, 0, 2e3), dims=(3,3) )
rel=CoefficientFunction( ( 1200, 0, 0,   0, 32000, 0,  0, 0, 1200), dims=(3,3) )


#a = BilinearForm(fes, symmetric=True)
a = BilinearForm(fes)
a += (1/mu0)*curl(u)*curl(v)*dx('air|coil') + rel*curl(u)*curl(v)*dx('core') + \
    1j*omega*sigma*u*v*dx('core') + 1j*1e-3*u*v*dx('air|coil')
#a += nu*curl(u)*curl(v)*dx + 1e-6*nu*u*v*dx

#c = Preconditioner(a, type="bddcc")
c = Preconditioner(a, type="multigrid")

#Js=CF((y,0.05-x,0))
#Draw(Js, mesh, 'Js')

f = LinearForm(fes)
f += CoefficientFunction((y,0.05-x,0)) *1e2* v * dx("coil")

sol = GridFunction(fes)

a.Assemble()
f.Assemble()

solver = CGSolver(mat=a.mat, pre=c, maxsteps=20000)
sol.vec.data = solver * f.vec
#sol.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec
J = - 1j * omega*sig * sigma * sol

Draw (curl(sol), mesh, "B", draw_surf=False)
Draw (J, mesh, "J", draw_surf=False)