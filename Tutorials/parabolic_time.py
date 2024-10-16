#imports
from ngsolve import *
from netgen.geom2d import SplineGeometry
#from ngsolve.webgui import Draw

from netgen.occ import *
#from netgen.webgui import Draw as DrawGeo

shape = Rectangle(2,2).Face().Move((-1,-1,0))
shape.edges.Min(X).name="left"
shape.edges.Max(X).name="right"
shape.edges.Min(Y).name="bottom"
shape.edges.Max(Y).name="top"
mesh = Mesh(OCCGeometry(shape, dim=2).GenerateMesh(maxh=0.25))

fes = H1(mesh, order=3, dirichlet="bottom|right|left|top")
u,v = fes.TnT()
time = 0.0
dt = 0.001

b = CoefficientFunction((2*y*(1-x*x),-2*x*(1-y*y)))
Draw(b,mesh,"wind", vectors={"grid_size": 32}, order=3)

a = BilinearForm(fes, symmetric=False)
a += 0.01*grad(u)*grad(v)*dx + b*grad(u)*v*dx
a.Assemble()

m = BilinearForm(fes, symmetric=False)
m += u*v*dx
m.Assemble()

mstar = m.mat.CreateMatrix()
print(f"m.mat.nze = {m.mat.nze}, a.mat.nze={a.mat.nze}, mstar.nze={mstar.nze}")


print(f"mstar.nze={mstar.nze}, len(mstar.AsVector())={len(mstar.AsVector())}")
#print(mstar.AsVector())
mstar.AsVector().data = m.mat.AsVector() + dt * a.mat.AsVector() # M* = M + dt * A
invmstar = mstar.Inverse(freedofs=fes.FreeDofs())

f = LinearForm(fes)
gaussp = exp(-6*((x+0.5)*(x+0.5)+y*y))-exp(-6*((x-0.5)*(x-0.5)+y*y))
Draw(gaussp,mesh,"f", deformation=True)
f += gaussp*v*dx
f.Assemble()

gfu = GridFunction(fes)
gfu.Set((1-y*y)*x) # note that boundary conditions remain
scene = Draw(gfu,mesh,"u")

def TimeStepping(invmstar, initial_cond = None, t0 = 0, tend = 2,
                 nsamples = 10):
    if initial_cond:
        gfu.Set(initial_cond)
    cnt = 0; time = t0
    sample_int = int(floor(tend / dt / nsamples)+1)
    gfut = GridFunction(gfu.space,multidim=0)
    gfut.AddMultiDimComponent(gfu.vec)
    while time < tend - 0.5 * dt:
        res = dt * f.vec - dt * a.mat * gfu.vec
        gfu.vec.data += invmstar * res #uoci operaciju +=
        print("\r",time,end="")
        #scene.Redraw()
        if cnt % sample_int == 0:
            gfut.AddMultiDimComponent(gfu.vec)
        cnt += 1; time = cnt * dt
    return gfut

gfut = TimeStepping(invmstar, (1-y*y)*x)

Draw(gfut, mesh, "gfut", interpolate_multidim=True, animate=True)

print(' ')
print(len(gfut.vecs)) #broj multidim komponenata
gfu.vec.data=gfut.vecs[7]
Draw(gfu, mesh, "gfut_7")



#from helper import ShowPattern
#print("sparsity pattern a.mat:")
#ShowPattern(a.mat,binarize=True)

#print("sparsity pattern mstar:")
#ShowPattern(mstar,binarize=True)