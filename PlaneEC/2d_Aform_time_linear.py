from ngsolve import *
from netgen.occ import *

outer= Circle((0,0), 0.2).Face()
outer.edges.name = 'rub'

inner = MoveTo(-0.08,-0.08).Line(0.16,0.0).Line(0,0.16).Line(-0.16,0).Close().Face()
#inner = MoveTo(-0.8,-0.8).Line(1.6,-0.1).Line(0,1.7).Line(-1.6,-0.2).Close().Face()


inner.edges.name="interface"
inner.faces.maxh=0.01
outer = outer - inner

inner.faces.name="inner"
inner.faces.col = (1, 1, 0) 
outer.faces.name="outer"


geo = Glue([inner, outer])
#Draw(OCCGeometry(geo));

mesh = Mesh(OCCGeometry(geo, dim=2).GenerateMesh(maxh=0.05, quad_dominated=False))

fes = HCurl(mesh, order=0, dirichlet="rub",  complex=False, nograds = False)
print('mesh.nv',mesh.nv)
print('mesh.nedge', mesh.nedge)
print('mesh.ne', mesh.ne)

u, v = fes.TnT()
time = 0.0
omega=2*pi*1
dt = 2*pi/omega /40


mu0 = 1.257e-6
rel = 200
sigma = 2e3
sig = mesh.MaterialCF({ "inner" : 1 }, default=None)

a = BilinearForm(fes, symmetric=False)
a += (1/mu0)*curl(u)*curl(v)*dx('outer') + rel*curl(u)*curl(v)*dx('inner')
a.Assemble()

m= BilinearForm(fes, symmetric=False)
m += sigma*u*v*dx('inner') + 1e-3*u*v*dx('outer') 
m.Assemble() 

#+1j*omega*sigma*u*v*dx('inner') + 1j*1e-3*u*v*dx('outer') 

mstar = m.mat.CreateMatrix()
print(f"m.mat.nze = {m.mat.nze}, a.mat.nze={a.mat.nze}, mstar.nze={mstar.nze}")
print(f"mstar.nze={mstar.nze}, len(mstar.AsVector())={len(mstar.AsVector())}")

mstar.AsVector().data = m.mat.AsVector() + dt * a.mat.AsVector()
invmstar = mstar.Inverse(freedofs=fes.FreeDofs())

sila=CF((0,0))
f = LinearForm(fes)
f += sila*v*dx
f.Assemble()


gfu = GridFunction(fes)

t = Parameter(0.0)

bc = CF((0.1*y*cos(omega*t),-0.1*x*cos(omega*t)))  #
time = 0.0
t.Set(0.0)
gfu.Set(bc, definedon=mesh.Boundaries('rub'))
#Draw(gfu, mesh, "gfu_bnd")


def TimeSteppingSin(invmstar, initial_cond = None, t0 = 0, tend = 2):
    if initial_cond:
        gfu.Set(initial_cond)
    cnt = 0; time = t0
    #sample_int = int(floor(tend / dt / nsamples)+1)
    gfuD = GridFunction(gfu.space)
    DgfuDt= GridFunction(gfu.space) 
    gfuPrev=GridFunction(gfu.space)
    
    gfut = GridFunction(gfu.space,multidim=0)
    Jgfut = GridFunction(gfu.space,multidim=0)
    
    gfut.AddMultiDimComponent(gfu.vec)
    while time < tend - 0.5 * dt:
        t.Set(time)
        #print(bc(mesh(0.198, 0.0)))
        gfuD.Set(bc,definedon=mesh.Boundaries('rub'))
        res = m.mat * gfu.vec - mstar *gfuD.vec
        gfu.vec.data = gfuD.vec + invmstar * res #uoci znak "=" umjesto "+="
        DgfuDt.vec.data = (gfu.vec - gfuPrev.vec)
        gfuPrev.vec.data=gfu.vec
        
        period=2*pi/omega
        if time>1.3*period and time<1.4*period :
            gfut.AddMultiDimComponent(gfu.vec)
            Jgfut.AddMultiDimComponent(DgfuDt.vec)
        cnt += 1; time = cnt * dt
        #print(cnt)
    return gfut, Jgfut


gfut, Jgfut = TimeSteppingSin(invmstar, initial_cond=CF((0,0)),tend=2)


####....... P O S T P R O C E S I N G........
##
comp=2

gfu.vec.data=gfut.vecs[comp]
B=curl(gfu)
Draw(B, mesh, f"B_{comp}")

Jgfu=GridFunction(fes)
Jgfu.vec.data=Jgfut.vecs[comp] 
J= Jgfu*sigma*(-1) * sig
Draw(J, mesh, f"J_{comp}")


#Pow=0.5*E*Conj(J) 
#Peddy=Integrate(Pow, mesh, order=5)
#print('Peddy',Peddy) 


#alternativa:
def TimeStepping(invmstar, initial_cond = None, t0 = 0, tend = 2,
                      nsamples = 10):
    if initial_cond:
        gfu.Set(initial_cond)
    cnt = 0; time = t0
    sample_int = int(floor(tend / dt / nsamples)+1)
    gfuD = GridFunction(gfu.space)
    gfut = GridFunction(gfu.space,multidim=0)
    gfut.AddMultiDimComponent(gfu.vec)
    while time < tend - 0.5 * dt:
        t.Set(time)
        #print(bc(mesh(0.198, 0.0)))
        gfuD.Set(bc,definedon=mesh.Boundaries('rub'))
        res = m.mat * gfu.vec - mstar *gfuD.vec
        gfu.vec.data = gfuD.vec + invmstar * res #uoci znak "=" umjesto "+="
        #print("\r",time,end="")
        if cnt % sample_int == 0:
            gfut.AddMultiDimComponent(gfu.vec)
        cnt += 1; time = cnt * dt
        #print(cnt)
    return gfut


#gfut = TimeStepping(invmstar, initial_cond=CF((0,0)),tend=1)