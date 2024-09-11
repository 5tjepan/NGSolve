from netgen.csg import *
from netgen.occ import *

from ngsolve import *

""" geometry = CSGeometry()
w=0.8
box = OrthoBrick(Pnt(-w,-w,-w),Pnt(w,w,w)).bc("outer")
geometry.Add(box) 
ngmesh = geometry.GenerateMesh(maxh=0.27)
mesh = Mesh(ngmesh) """

w=0.8
geom = MoveTo(-w, -w).Rectangle(2*w,w).Face()
mesh = Mesh(OCCGeometry(geom, dim=2).GenerateMesh(maxh=0.1, quad_dominated=True))


Draw(mesh)
print(mesh.GetBoundaries())

fes = HCurl(mesh, order=0, dirichlet="outer", complex=True, nograds = False)
mvp = fes.TrialFunction()
alpha = fes.TestFunction()

#func=CF(((x*y)**3+x-z**2+4, x*y*z,x**3))  OVO JE ZA 3D
func=CF(((x*y)**3+x-y**2+4, x*y*2)) #OVO JE ZA 2D
sol = GridFunction(fes)
sol.Set(func)

rot=curl(sol)
#Tpot.Set(CF(0),BND)

#solvers.BVP(bf=a, lf=f, gf=sol, pre=None, maxsteps=2, print=True, needsassembling=False)

#end = time.time()
#print(end-start)


Draw (sol, mesh, "sol")
Draw (func, mesh, "f")
Draw (rot, mesh, "rot")

