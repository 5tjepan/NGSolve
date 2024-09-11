from ngsolve import *
from netgen.occ import *

outer= Circle((1,1), 2).Face()
outer.edges.name = 'rub'
""" outer = Rectangle(2, 2).Face()
outer.edges.name="outer"
outer.edges.Max(X).name = "r"
outer.edges.Min(X).name = "l"
outer.edges.Min(Y).name = "b"
outer.edges.Max(Y).name = "t"
outer.edges.Min(X).maxh=0.1 """

inner1 = MoveTo(0.4, 0.5).Rectangle(0.4, 1).Face()
#inner = MoveTo(1,0).Line(1,1).Line(-1,1).Line(-1,-1).Close().Face()
inner1.edges.name="interface"
inner1.faces.maxh=0.05

inner2 = MoveTo(1.1, 0.5).Rectangle(0.4, 1).Face()
#inner = MoveTo(1,0).Line(1,1).Line(-1,1).Line(-1,-1).Close().Face()
inner2.edges.name="interface"
inner2.faces.maxh=0.05

inner = inner1 + inner2
outer = outer - inner

inner.faces.name="inner"
inner.faces.col = (1, 1, 0)  #colour
outer.faces.name="outer"

geo = Glue([inner, outer])
#Draw(OCCGeometry(geo));

mesh = Mesh(OCCGeometry(geo, dim=2).GenerateMesh(maxh=0.3, quad_dominated=True))
fes = HCurl(mesh, order=0, dirichlet="rub",  complex=True, nograds = False)
gfu = GridFunction(fes)
print('mesh.nv',mesh.nv)
print('mesh.nedge', mesh.nedge)

u, v = fes.TnT()

bb=CF((10*y,-10*x))
#gfu.Set(bb, VOL_or_BND = BND)
gfu.Set(bb, definedon=mesh.Boundaries('rub'))
#gfu.Set(bb, definedon=mesh.Boundaries('b'))


omega=314
mu0 = 1.257e-6
rel = 1200
#sigma = 2e3
sig = mesh.MaterialCF({ "inner" : 1 }, default=None)
sigma=CoefficientFunction( (10, 0,  0, 10), dims=(2,2) )

#rel=CoefficientFunction( ( 1200, 0, 0, 1200), dims=(2,2) )

a = BilinearForm(fes)
a += (1/mu0)*curl(u)*curl(v)*dx('outer') + rel*curl(u)*curl(v)*dx('inner') + \
    1j*omega*sigma*u*v*dx('inner') + 1j*1e-3*u*v*dx('outer')

f=LinearForm(fes)
f += CF((0,0))*v*dx

a.Assemble()
f.Assemble()

r = f.vec - a.mat * gfu.vec
#gfu.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs())*r

solvers.BVP(bf=a, lf=f, gf=gfu, maxsteps=200, print=True)

#solver = CGSolver(mat=a.mat, maxsteps=80000, precision=1e-14)
#gfu.vec.data += solver * f.vec

Draw(gfu)

J = - 1j * omega*sig * sigma * gfu
B = curl(gfu)

Draw (B, mesh, "B")
Draw (J, mesh, "J")
