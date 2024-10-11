from ngsolve import *

mesh=Mesh(unit_square.GenerateMesh(maxh=0.3))
Draw(mesh)
fes = H1(mesh, order=1, dirichlet='bottom|right|top|left') #treba bi cijela granica ici?

U = x*x*(1-y)*(1-y)          # U = manufactured solution
DeltaU = 2*((1-y)*(1-y)+x*x) # Source: DeltaU = ∆U

u, v =fes.TnT()
a=BilinearForm(grad(u)*grad(v)*dx).Assemble()

f=LinearForm(-DeltaU*v*dx).Assemble()

gfu=GridFunction(fes)  #vektor rjesenja s nul vrijednostima
gfu.Set(U, BND)  #vektor rjesenja s tocnim rubnim vrijednostima
print(gfu.vec)

invA=a.mat.Inverse(freedofs=fes.FreeDofs())
gfu.vec.data += invA*(f.vec-a.mat*gfu.vec)

Draw(gfu)

print('error=', sqrt(Integrate((U-gfu)*(U-gfu),mesh)))  # Compute error
