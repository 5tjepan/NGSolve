from ngsolve import *
from netgen.occ import *

outer = Rectangle(2, 2).Face()
outer.edges.name="outer"
outer.edges.Max(X).name = "r"
outer.edges.Min(X).name = "l"
outer.edges.Min(Y).name = "b"
outer.edges.Max(Y).name = "t"
outer.edges.Min(X).maxh=0.1

inner = MoveTo(1, 0.9).Rectangle(0.4, 0.5).Face()
inner.edges.name="interface"
inner.faces.maxh=0.1
outer = outer - inner

inner.faces.name="inner"
inner.faces.col = (1, 1, 0)  #colour
outer.faces.name="outer"

geo = Glue([inner, outer])
#Draw(OCCGeometry(geo));

mesh = Mesh(OCCGeometry(geo, dim=2).GenerateMesh(maxh=0.3))

print("material names", mesh.GetMaterials())
print("boundary names", mesh.GetBoundaries())
Draw(mesh);

#fes1 = L2(mesh, definedon="inner", order=0)
fes1 = H1(mesh, definedon="inner", order=2)
u1 = GridFunction(fes1, "u1")
u1.Set(x*y)
#
Draw(u1)

#fes1comp = Compress(fes1)
#print('compress',fes1comp.ndof)

fes = H1(mesh, order=2, dirichlet='b|l|r')
u,v =fes.TnT()
sol = GridFunction(fes)

#term1= u1*v*dx('inner')
term1= u1*v*dx(definedon=mesh.Materials('inner'))
#term2 = 0.1*v*ds(definedon=mesh.Material('t'))
term2 = 0.1*v*ds('t')
f=LinearForm(term1+term2).Assemble()

a=BilinearForm(fes)
a+= grad(u)*grad(v)*dx
a.Assemble()
sol.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec
Draw(sol)

print(mesh.Materials('inner|outer').Mask())
print(mesh.Materials("[a-z]*").Mask())  # can use regexp
print(~mesh.Boundaries('t|l').Mask())

io = mesh.Materials("inner") + mesh.Materials("outer")
print(io.Mask())

domain_values = {'inner': 3.7,  'outer': 1}
cf = mesh.MaterialCF(domain_values)
Draw(cf, mesh,'cf')

bdry_values = {'b': x, 'r': 2-y}
cf = mesh.BoundaryCF(bdry_values, default=0)
Draw(cf, mesh,'cf2')

g = GridFunction(H1(mesh), name='bdry')
g.Set(cf, definedon=mesh.Boundaries('b|r'))
Draw(g);

