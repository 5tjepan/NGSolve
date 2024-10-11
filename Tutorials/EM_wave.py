from ngsolve import *
from netgen.occ import *
from netgen.geom2d import SplineGeometry

air= Circle((0.5,0.5), 0.8).Face()
air.edges.name = 'outer'
scatter= MoveTo(0.7,0.3).Rectangle(0.05,0.4).Face()
scatter.edges.name = 'scat'
geo = OCCGeometry(air -scatter, dim=2)
mesh = Mesh(geo.GenerateMesh(maxh=0.05))
Draw(mesh)

fes = H1(mesh, order=5, complex=True) #COMPLEX True
u, v = fes.TnT()

omega = 100
pulse = 5e4*exp(-(40**2)*((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5))) #vrlo strma CF funkcija jer je puls
Draw(pulse, mesh, 'pulse',order=3)

""" a1=grad(u)*grad(v)*dx - omega**2*u*v*dx
a2= - omega *1j*u*v *ds('outer')
a = BilinearForm(a1 + a2).Assemble() """
a = BilinearForm(fes)
a += grad(u)*grad(v)*dx - omega**2*u*v*dx
a += -omega*1j*u*v * ds("outer")
a.Assemble()

f= LinearForm(pulse * v *dx).Assemble()

#sol = GridFunction(fes, name="sol")
#print(fes.FreeDofs()) #printa kaze da su svi Dofs slobodni
#sol.vec.data = a.mat.Inverse() * f.vec
#Draw(sol, mesh, min=-1, max=1, order=3, animate_complex=True)
#Draw(sol, mesh, autoscale= False, min=-1, max=1, order=3)
#Draw(sol, mesh)

gfu = GridFunction(fes, name="u")
gfu.vec.data = a.mat.Inverse() * f.vec
Draw(gfu, min=-1, max=1, order=3, animate_complex=True);

