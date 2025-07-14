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
print('mesh.nv',mesh.nv)
print('mesh.nedge', mesh.nedge)
print('mesh.ne', mesh.ne)

#.............
fes = HCurl(mesh, order=0, dirichlet="rub",  complex=False, nograds = False)
u, v = fes.TnT()

omega=2*pi*1
dt = 2*pi/omega /100
cnt = 0; t0=0.0 ; time = t0; tend=2.1
interval=(0.4,1.0); samplestep=2

mu0 = 1.257e-6
#rel = 200
sigma = 2e3
sig = mesh.MaterialCF({ "core" : 1 }, default=None)
Kf=1

#-------------
B_ref=[0.0, 0.04, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.83, 1000]
H_ref=[0.0, 3.5, 5.0, 7.6, 10.0, 12.07, 14.08, 16.0, 17.75, 19.42, 21.05, 22.64, 24.16, 25.76, 27.98, 29.66, 32.05, 35.33, 40.06, 47.17, 58.35, 77.19, 113.0, 196.0, 395.0, 800, 1000/mu0]

HBcurve = BSpline(2, [0]+list(B_ref), list(H_ref)) #ovo bi trebala biti instanca klase BSpline
diffHB= HBcurve.Differentiate() #HBcurve.Differentiate() metoda daje objekt klase BSpline
#-------------

gfu = GridFunction(fes)
oldgfu = GridFunction(fes)
t = Parameter(0.0)
#...
ic=CF((0,0)) #initial condition
gfu.Set(ic)
oldgfu.Set(CF((0,0)))
#..
bc = 1*CF((0.1*y*sin(omega*t),-0.1*x*sin(omega*t))) #boundary condition
gfu.Set(bc, definedon=mesh.Boundaries('rub'))
#Draw(gfu, mesh, "gfu_bnd")

#inicijalizacija listi:::
time_axis=[]
Blist=[]
Jlistoftuple=[]
Plist=[]
Bavglist=[]


#inicijalizacija gridfunkcija:::
gfuD = GridFunction(gfu.space) #dirichlet
DgfuDt= GridFunction(gfu.space) #dA/dt
gfuPrev=GridFunction(gfu.space)
gfuPrev.vec.data=gfu.vec

gfut = GridFunction(gfu.space,multidim=0) #objekt za spremanje gfu u razlicitim trenucima; slovo t znaci time
Jgfut = GridFunction(gfu.space,multidim=0) #objekt za spremanje Jgfu u razlicitim trenucima; t znaci time
gfut.AddMultiDimComponent(gfu.vec) #u prvom stupcu je inicijalno stanje gfu(t=0)
Jgfut.AddMultiDimComponent(gfu.vec)

#............................
while time < tend - 0.5 * dt:  #Euler time-stepping petlja
#for iter in range(1,14):  #Euler time-stepping petlja

    t.Set(time)
    gfuD.Set(bc,definedon=mesh.Boundaries('rub')) #vrem. ovisan bc (t.Set(time))
    errorlist=[]
    
    Babs= curl(gfu).Norm()
    
    for it in range(1, 6):  #NewtonRaphson petlja

        #rel=200
        #dHdB=200
        rel = (HBcurve(Babs+1e-6))/(Babs+1e-6) #+ 1j*omega*sigma*d**2*1/12 #Babs+1e-6
        dHdB = diffHB(Babs+1e-6) #+ 1j*omega*sigma*d**2*1/12
        print('rel=',rel(mesh(0,0.495*brid)), 'dHdB=',dHdB(mesh(0,0.495*brid)))

        a = BilinearForm(fes, symmetric=False)
        #a += (1/mu0)*curl(u)*curl(v)*dx('outer') + rel*curl(u)*curl(v)*dx('core') + (dHdB - rel)*curl(u)*curl(v)*dx('core')
        term1 = (1/mu0)*curl(u)*curl(v)*dx('outer') + rel*curl(u)*curl(v)*dx('core')
        jac=(dHdB - rel)*curl(u)*curl(v)*dx('core')
        a = BilinearForm(term1 + jac)
        a.Assemble()
 
        jacmat= BilinearForm(jac)
        jacmat.Assemble()

        m= BilinearForm(fes, symmetric=False)
        m += sigma*u*v*dx('core') + 1e-3*u*v*dx('outer') 
        m.Assemble() 

        mstar = m.mat.CreateMatrix()
        #print(f"m.mat.nze = {m.mat.nze}, a.mat.nze={a.mat.nze}, mstar.nze={mstar.nze}")
        #print(f"mstar.nze={mstar.nze}, len(mstar.AsVector())={len(mstar.AsVector())}")

        mstar.AsVector().data = m.mat.AsVector() + dt * a.mat.AsVector()
        invmstar = mstar.Inverse(freedofs=fes.FreeDofs())

        #sila=CF((0,0))
        #f = LinearForm(fes)
        #f += sila*v*dx
        #f.Assemble()

        res = m.mat * gfuPrev.vec - mstar *gfuD.vec + dt*jacmat.mat * oldgfu.vec
        gfu.vec.data = gfuD.vec + invmstar * res

        #NR pogreska:
        errfunc = (gfu - oldgfu).Norm()/(oldgfu.Norm()+1e-15) #/oldApot
        
        defon = mesh.Materials('core')
        error=Integrate(errfunc.Norm(), mesh, definedon=defon)
        #print('error =', error)
        errorlist.append(error)

        #oldgfu.vec.data= p*gfu.vec + (1-p)*old.vec  #damping
        oldgfu.vec.data= gfu.vec
        
        Babs= curl(gfu).Norm()
    
    print('errorlist', errorlist)

    DgfuDt.vec.data = (gfu.vec - gfuPrev.vec)
    gfuPrev.vec.data=gfu.vec
    Blist.append(curl(gfu)[0](mesh(0,0.495*brid)))
    Jlistoftuple.append(DgfuDt(mesh(0,0.495*brid)))
    #print(DgfuDt.Norm()(mesh(0,0.49*brid)))
    time_axis.append(time)
    #current power:
    currPow=Integrate(sigma*DgfuDt*DgfuDt/dt**2, mesh, 
                    definedon=mesh.Materials('core'))
    currBavg=Integrate(curl(gfu)/brid**2, mesh, definedon=mesh.Materials('core'))
    Plist.append(currPow)
    Bavglist.append(currBavg)
    
    if time>interval[0]*(tend-t0) and time<interval[1]*(tend-t0) and cnt%samplestep==0 :
        gfut.AddMultiDimComponent(gfu.vec)
        Jgfut.AddMultiDimComponent(DgfuDt.vec)
    cnt += 1; time = cnt * dt




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
pnt1=1
plt.figure()

plt.subplot(311)
plt.plot(time_axis[pnt1:],Blist[pnt1:],'b')
plt.grid(True)
#...
plt.subplot(312)
plt.plot(time_axis[pnt1:], Jlist[pnt1:],'r--')
plt.grid(True)

plt.subplot(313)
plt.plot(time_axis[pnt1:], Plist[pnt1:],'g')
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


