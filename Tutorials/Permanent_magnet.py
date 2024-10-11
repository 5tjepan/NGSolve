from ngsolve import *
from netgen.occ import *

# box = OrthoBrick(Pnt(-3,-3,-3),Pnt(3,3,3)).bc("outer")
box = Box((-3,-3,-3), (3,3,3))
box.faces.name = 'outer'

#magnet = Cylinder(Pnt(-1,0,0),Pnt(1,0,0), 0.3)
magnet = Cylinder((-1,0,0),X, r=0.3, h=2)
magnet.mat('magnet') #.mat znaci material
magnet.faces.col =(1,0,0)
magnet.faces.name = 'magBND'
magnet.faces.maxh=0.1

air = box - magnet
air.mat('air')
shape = Glue([air,magnet]) #ovo je 'netgen.libngpy._NgOCC.TopoDS_Shape' objekt
geo=OCCGeometry(shape) #ovo je netgen.libngpy._NgOCC.OCCGeometry objekt

Draw(geo)

#mesh = Mesh(geo.GenerateMesh(maxh=2))
mesh = Mesh(geo.GenerateMesh(maxh=2, curvaturesafety=1))
mesh.Curve(3);

print(mesh.GetMaterials())
print(mesh.GetBoundaries())

fes = HCurl(mesh, order=1, dirichlet ='outer', nograds=True)
print('ndof', fes.ndof)
u,v = fes.TnT()

from math import pi
mu0 = 4*pi*1e-7
mur= mesh.MaterialCF({'magnet': 1000}, default=1)

a = BilinearForm(fes)
a+= 1/(mu0*mur)*curl(u)*curl(v)*dx + 1e-8/(mu0*mur)*u*v*dx
a.Assemble()

f=LinearForm(fes)
mag = mesh.MaterialCF({'magnet': (10,0,0)}, default=(0,0,0))
f+= mag*curl(v)*dx('magnet')
f.Assemble()

freedofss=fes.FreeDofs()
gfu= GridFunction(fes)
#gfu.vec.data = a.mat.Inverse(freedofs=fes.FreeDofs()) * f.vec
#solvers.CG(sol=gfu.vec, rhs=f.vec, mat=a.mat, maxsteps=200, tol=1e-5, printrates=False, freedofs=freedofss)
matrica = CGSolver(mat=a.mat, pre=None ,maxsteps=1000, printrates=True, precision=3e-05)
gfu.vec.data = matrica * f.vec 

Draw(gfu)
rot=curl(gfu);
Draw(rot, mesh,'B')
