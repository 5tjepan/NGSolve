from ngsolve import *
from netgen.occ import *

outer= Circle((0,0), 0.2).Face()
outer.edges.name = 'rub'


inner = MoveTo(-0.08,-0.08).Line(0.16,0.0).Line(0,0.16).Line(-0.16,0).Close().Face()
#inner = MoveTo(-0.8,-0.8).Line(1.6,-0.1).Line(0,1.7).Line(-1.6,-0.2).Close().Face()

inner.edges.name="interface"
inner.faces.maxh=0.02
outer = outer - inner

inner.faces.name="inner"
inner.faces.col = (1, 1, 0)  #colour
outer.faces.name="outer"

geo = Glue([inner, outer])
#Draw(OCCGeometry(geo));

mesh = Mesh(OCCGeometry(geo, dim=2).GenerateMesh(maxh=0.1, quad_dominated=False))
fes = HCurl(mesh, order=0, dirichlet="rub",  complex=True, nograds = False) #CMPLX
gfu = GridFunction(fes)
print('mesh.nv',mesh.nv)
print('mesh.nedge', mesh.nedge)
print('mesh.ne', mesh.ne)

u, v = fes.TnT()

bb=CF((0.1*y,-0.1*x))
gfu.Set(bb, definedon=mesh.Boundaries('rub'))

omega=3100
mu0 = 1.257e-6
rel = (0.2-x)*1000
sigma = 2e1
sig = mesh.MaterialCF({ "inner" : 1 }, default=None)
#sigma=CoefficientFunction( (10, 0,  0, 10), dims=(2,2) )

#rel=CoefficientFunction( ( 1200, 0, 0, 1200), dims=(2,2) )

a = BilinearForm(fes)
a += (1/mu0)*curl(u)*curl(v)*dx('outer') + rel*curl(u)*curl(v)*dx('inner') + \
    +0.1*u*v*dx('outer') + 1j*omega*sigma*u*v*dx('inner') 

a.Assemble()

sila=CF((0,0))
f = LinearForm(fes)
f += sila*v*dx
f.Assemble()

#solvers.BVP(bf=a, lf=f, gf=gfu, pre=None, maxsteps=200, print=True)

r = - a.mat * gfu.vec
gfu.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs())*r

A=gfu

B = curl(gfu)
E = - 1j * omega * A
J = - 1j * omega * sigma * A *sig

Draw(A, mesh, "A")
Draw (B, mesh, "B")
Draw (J, mesh, "J")



""" testu = GridFunction(fes)
testu.vec[:] = 0
k=43
testu.vec[k]=gfu.vec[k]
print('k=',k,'testu.vec[k]=',testu.vec[k])
Draw(testu)

 """
