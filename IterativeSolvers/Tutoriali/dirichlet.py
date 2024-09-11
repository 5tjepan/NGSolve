from ngsolve import *
from netgen.geom2d import unit_square

mesh = Mesh(unit_square.GenerateMesh(maxh=0.3))
mesh.GetBoundaries()

fes = H1(mesh, order=2, dirichlet="left|right")
#print("free dofs of fes:\n", fes.FreeDofs())
gfu = GridFunction(fes)

g = sin(y)

#gfu.Set(coefficient=g, VOL_or_BND = VOL)
gfu.Set(coefficient=g, VOL_or_BND = BND)

###
# u = fes.TrialFunction()  # symbolic object
# v = fes.TestFunction()  
u, v =fes.TnT()
a = BilinearForm(grad(u)*grad(v)*dx).Assemble()
f = LinearForm(1*v*dx).Assemble()
r = f.vec - a.mat * gfu.vec

gfu.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs())*r

Draw(gfu)
