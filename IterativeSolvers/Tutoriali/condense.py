from ngsolve import *
mesh = Mesh(unit_square.GenerateMesh(maxh=0.4))
Draw(mesh)

fes = H1(mesh, order=3, dirichlet='bottom|right')
#fes = H1(mesh, order=4)
u, v = fes.TnT()

a=BilinearForm(grad(u)*grad(v)*dx, condense=False).Assemble()
#a=BilinearForm(grad(u)*grad(v)*dx, condense=True).Assemble()
f=LinearForm(x*v*dx).Assemble()
u=GridFunction(fes)

u.vec.data = a.mat.Inverse(freedofs=fes.FreeDofs()) * f.vec
#u.vec.data = a.mat.Inverse(freedofs=fes.FreeDofs(coupling=True)) * f.vec

""" invS = a.mat.Inverse(freedofs=fes.FreeDofs(coupling=True))+a.inner_solve
ext = IdentityMatrix() + a.harmonic_extension
extT = IdentityMatrix() + a.harmonic_extension_trans
#invA =  ext @ invS @ extT + a.inner_solve
invA =  ext @ invS @ extT
u.vec.data = invA * f.vec """

Draw(u)