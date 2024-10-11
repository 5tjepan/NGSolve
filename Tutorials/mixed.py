from ngsolve import *

mesh = Mesh(unit_square.GenerateMesh(maxh=0.1))

source = sin(3.14*x)
ud = CF(5)
g = mesh.BoundaryCF( {"left" : y*(1-y)}, default=0)
lam = 10

#Draw(g,mesh,'g')

order_flux=1
V = HDiv(mesh, order=order_flux, dirichlet="right|top|left")
Q = L2(mesh, order=order_flux-1)
fesm = V*Q

sigma, u = fesm.TrialFunction()
tau, v = fesm.TestFunction()

normal = specialcf.normal(mesh.dim)
print (normal)

am = BilinearForm((1/lam*sigma*tau + div(sigma)*v + div(tau)*u)*dx).Assemble()
fm = LinearForm(-source*v*dx + ud*(tau.Trace()*normal)*ds).Assemble()

gfm = GridFunction(fesm)

gfsigma, gfu = gfm.components

gfsigma.Set(g*normal, BND)


res = fm.vec.data - am.mat * gfm.vec
gfm.vec.data += am.mat.Inverse(freedofs=fesm.FreeDofs(), inverse="umfpack") * res
#matrica = CGSolver(mat=am.mat, pre=None ,maxsteps=1000, printrates=True, precision=3e-05)
#gfm.vec.data += matrica *res
# solvers.BVP(bf=am, lf=fm, gf=gfm)

Draw (gfsigma, mesh, "flux-mixed")
Draw (gfu, mesh, "u-mixed");



