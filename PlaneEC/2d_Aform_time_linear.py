from ngsolve import *
from netgen.occ import *

outer= Circle((0,0), 0.2).Face()
outer.edges.name = 'rub'
brid=0.16
core = MoveTo(-0.08,-0.08).Line(brid,0.0).Line(0,brid).Line(-brid,0).Close().Face()
core.edges.name="interface"
core.faces.maxh=0.01
outer = outer - core
core.faces.name="core"
core.faces.col = (1, 1, 0) 
outer.faces.name="outer"

geo = Glue([core, outer])
#Draw(OCCGeometry(geo));

mesh = Mesh(OCCGeometry(geo, dim=2).GenerateMesh(maxh=0.05, quad_dominated=False))

fes = HCurl(mesh, order=0, dirichlet="rub",  complex=False, nograds = False)
print('mesh.nv',mesh.nv)
print('mesh.nedge', mesh.nedge)
print('mesh.ne', mesh.ne)

u, v = fes.TnT()
time = 0.0
omega=2*pi*1
dt = 2*pi/omega /80

mu0 = 1.257e-6
rel = 200
sigma = 2e4
sig = mesh.MaterialCF({ "core" : 1 }, default=None)

a = BilinearForm(fes, symmetric=False)
a += (1/mu0)*curl(u)*curl(v)*dx('outer') + rel*curl(u)*curl(v)*dx('core')
a.Assemble()

m= BilinearForm(fes, symmetric=False)
m += sigma*u*v*dx('core') + 1e-3*u*v*dx('outer') 
m.Assemble() 

#+1j*omega*sigma*u*v*dx('core') + 1j*1e-3*u*v*dx('outer') 

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

bc = CF((0.1*y*cos(omega*t),-0.1*x*cos(omega*t))) #*sqrt(2)
time = 0.0
t.Set(0.0)
gfu.Set(bc, definedon=mesh.Boundaries('rub'))
#Draw(gfu, mesh, "gfu_bnd")

time_axis=[]
Blist=[]
Jlistoftuple=[]
Plist=[]
def TimeStepping(invmstar, initial_cond = None, t0 = 0, tend = 2,
                     range=(0.8,0.9), samplestep=2):
    if initial_cond:
        gfu.Set(initial_cond)
    cnt = 0; time = t0
    gfuD = GridFunction(gfu.space)
    DgfuDt= GridFunction(gfu.space) 
    gfuPrev=GridFunction(gfu.space)
    gfuPrev.vec.data=gfu.vec
    
    gfut = GridFunction(gfu.space,multidim=0)
    Jgfut = GridFunction(gfu.space,multidim=0)
    gfut.AddMultiDimComponent(gfu.vec)
    Jgfut.AddMultiDimComponent(gfu.vec)
    while time < tend - 0.5 * dt:
        t.Set(time)
        #print(bc(mesh(0.198, 0.0)))
        gfuD.Set(bc,definedon=mesh.Boundaries('rub')) #vrem. ovisan bc (t.Set(time))
        res = m.mat * gfu.vec - mstar *gfuD.vec
        gfu.vec.data = gfuD.vec + invmstar * res #uoci znak "=" umjesto "+="
        DgfuDt.vec.data = (gfu.vec - gfuPrev.vec)
        gfuPrev.vec.data=gfu.vec
        Blist.append(curl(gfu)[0](mesh(0,0.495*brid)))
        Jlistoftuple.append(DgfuDt(mesh(0,0.495*brid)))
        #print(DgfuDt.Norm()(mesh(0,0.49*brid)))
        time_axis.append(time)
        #current power:
        currPow=Integrate(sigma*DgfuDt*DgfuDt/dt**2, mesh, 
                      definedon=mesh.Materials('core'))
        Plist.append(currPow)
        if time>range[0]*(tend-t0) and time<range[1]*(tend-t0) and cnt%samplestep==0 :
            gfut.AddMultiDimComponent(gfu.vec)
            Jgfut.AddMultiDimComponent(DgfuDt.vec)
        cnt += 1; time = cnt * dt
        #print(cnt)
    return gfut, Jgfut

gfut, Jgfut = TimeStepping(invmstar, initial_cond=CF((0,0)),tend=3,samplestep=3)


####....... P O S T P R O C E S I N G........
##

#Blist=[elem[0] for elem in Blistoftuple]
Jlist=[-1*sigma/dt*elem[0] for elem in Jlistoftuple] #elem[0] je Jx komponenta

period_samples=int(round(2*pi/omega/dt))
print('period_samples', period_samples)
Pavg=sum(Plist[-(period_samples+1):-1])/period_samples
print('Pavg=',Pavg)

#...V I Z U A L I Z A C I J A:::
import matplotlib.pyplot as plt
import sys
sys.argv = ["fun"]

plt.figure()

plt.subplot(311)
plt.plot(time_axis[10:],Blist[10:],'b')
plt.grid(True)
#...
plt.subplot(312)
plt.plot(time_axis[10:], Jlist[10:],'r--')
plt.grid(True)

plt.subplot(313)
plt.plot(time_axis[10:], Plist[10:],'g')
plt.grid(True)
plt.show()

#...
comp=1
gfu.vec.data=gfut.vecs[comp]
B=curl(gfu)
Draw(B, mesh, f"B_{comp}")

Jgfu=GridFunction(fes)
Jgfu.vec.data=Jgfut.vecs[comp] 
J= Jgfu*sigma*(-1) * sig
Draw(J, mesh, f"J_{comp}")


#alternativa samo za B:
def TimeSteppingOld(invmstar, initial_cond = None, t0 = 0, tend = 2,
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

#gfut = TimeSteppingOld(invmstar, initial_cond=CF((0,0)),tend=1)