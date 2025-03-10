from ngsolve import *
from netgen.occ import *
import numpy as np
from extremaseek import FindDominantExtrema
from everett import Everett_atan,Everett_exp, Everett_exp2
from preisach import preisach_output

outer= Circle((0,0), 0.15).Face()
outer.edges.name = 'rub'
brid=0.16
core = MoveTo(-0.08,-0.08).Line(brid,0.0).Line(0,brid).Line(-brid,0).Close().Face()
core.edges.name="interface"
core.faces.maxh=0.2
outer = outer - core
core.faces.name="core"
core.faces.col = (1, 1, 0) 
outer.faces.name="outer"

geo = Glue([core, outer])
#Draw(OCCGeometry(geo));

mesh = Mesh(OCCGeometry(geo, dim=2).GenerateMesh(maxh=0.15, quad_dominated=False))
print('mesh.nv',mesh.nv)
print('mesh.nedge', mesh.nedge)
print('mesh.ne', mesh.ne)

#.............
fes = H1(mesh, order=0, dirichlet="rub", complex=False)
H, v = fes.TnT()

omega=2*pi*1
dt = 2*pi/omega /100
cnt = 0; t0=0.0 ; time = t0; tend=2.5
interval=(0.4,1.0); samplestep=2

mu0 = 1.257e-6
#rel = 200
sigma = 2e3
sig = mesh.MaterialCF({ "core" : 1 }, default=None)
Kf=1

#-------------
mu0=4*pi*1e-7

B_ref=[0.001, 0.04, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.83, 1000]
H_ref=[0.001, 3.5, 5.0, 7.6, 10.0, 12.07, 14.08, 16.0, 17.75, 19.42, 21.05, 22.64, 24.16, 25.76, 27.98, 29.66, 32.05, 35.33, 40.06, 47.17, 58.35, 77.19, 113.0, 196.0, 395.0, 800, 1000/mu0]

BHcurve = BSpline(2, [0]+list(H_ref), list(B_ref)) #ovo bi trebala biti instanca klase BSpline
diffBH= BHcurve.Differentiate() #HBcurve.Differentiate() metoda daje objekt klase BSpline
#-------------

gfu = GridFunction(fes)
oldgfu = GridFunction(fes)
t = Parameter(0.0)
#...
ic=CF(0) #initial condition
gfu.Set(ic)
oldgfu.Set(CF(0))
#..
bc = 30*sin(omega*t) #boundary condition
gfu.Set(bc, definedon=mesh.Boundaries('rub'))
#Draw(gfu, mesh, "gfu_bnd")

#inicijalizacija listi:::
time_axis=[]
Hlist=[]
Blist=[]
Jlistoftuple=[]
Plist=[]
Bavglist=[]

#inicijalizacija gridfunkcija:::
gfuD = GridFunction(gfu.space) #dirichlet
#DgfuDt= GridFunction(gfu.space) #dA/dt
gfuPrev=GridFunction(gfu.space)
gfuPrev.vec.data=gfu.vec

oldB =GridFunction(fes)
prevB =GridFunction(fes)
prevB.Set(ic*mu0)

gfut = GridFunction(gfu.space,multidim=0) #objekt za spremanje gfu u razlicitim trenucima; slovo t znaci time
gfut.AddMultiDimComponent(gfu.vec) #u prvom stupcu je inicijalno stanje gfu(t=0)

rot=CF( (0 , 1,  -1, 0), dims=(2,2) )

#HabsPrev= gfuPrev.Norm()
#Habs= oldgfu.Norm()
    
prevHmatrix= np.array([gfu.vec.FV().NumPy()])

Everett = Everett_exp
#............................
while time < tend - 0.5 * dt:  #Euler time-stepping petlja
#for iter in range(1,14):  #Euler time-stepping petlja

    t.Set(time)
    gfuD.Set(bc,definedon=mesh.Boundaries('rub')) #vrem. ovisan bc (t.Set(time))
    errorlist=[]
    
    #HabsPrev= gfuPrev.Norm()
    #Habs= oldgfu.Norm()
    
    for it in range(1, 8):  #NewtonRaphson petlja

        """ mu_fe = (BHcurve(Habs+1e-6))/(Habs+1e-6) #+ 1j*omega*sigma*d**2*1/12 #Babs+1e-6
        mu_fe_prev = (BHcurve(HabsPrev+1e-6))/(HabsPrev+1e-6) #+ 1j*omega*sigma*d**2*1/12 #Babs+1e-6
        dBdH = diffBH(Habs+1e-6) #+ 1j*omega*sigma*d**2*1/12 """

        """ mu_fe= 0.1*(Habs+1e-4)**(-0.5)
        mu_fe_prev= 0.1* (HabsPrev +1e-4)**(-0.5)
        dBdH= 0.1* 0.5 *(Habs +1e-4)**(-0.5) """
        
        """ prevB = 0.1* (HabsPrev +1e-4)**(0.5) * IfPos(gfuPrev,1,-1)
        oldB = 0.1* (Habs+1e-4)**(0.5) * IfPos(oldgfu,1,-1) """

        #print('rel=',rel(mesh(0,0.495*brid)), 'dHdB=',dHdB(mesh(0,0.495*brid)))

        #mu = mesh.MaterialCF({ "core" : mu_fe }, default=mu0)
        #mu_prev = mesh.MaterialCF({ "core" : mu_fe_prev }, default=mu0)
        #mud = mesh.MaterialCF({ "core" : dBdH }, default=mu0)

        #oldHin=oldgfu*sig 
        #prevHin=gfuPrev*sig

        #===================
        # h y s t e r e z a 
        Everett = Everett_exp
        Hvec = np.array([oldgfu.vec.FV().NumPy()])
        oldHmat=np.append(prevHmatrix, Hvec, axis=0)
        Hmax=100
        oldBlist=[]
        for el in range(gfu.vec.size):
            Helem=oldHmat[:,el]
            domIndx=FindDominantExtrema(Helem,Hmax)
            oldBlist.append(preisach_output(Helem,domIndx,Hmax))

        oldB.vec.FV().NumPy()[:] = np.array(oldBlist)
        mud=(oldB-prevB)/(oldgfu-gfuPrev+1e-4)
        #==================

        a = BilinearForm(fes, symmetric=False)
        #a += (1/mu0)*curl(u)*curl(v)*dx('outer') + rel*curl(u)*curl(v)*dx('core') + (dHdB - rel)*curl(u)*curl(v)*dx('core')
        
        term1 = dt*grad(H)*grad(v)*dx + sigma*mud*H*v*dx('core')
        a = BilinearForm(term1)
        a.Assemble()
 
        #jac = sigma*(mud - mu)*v*H*dx('core')
        jac = sigma* mud *v*H*dx('core')
        jacmat= BilinearForm(jac)
        jacmat.Assemble()

        
        #mstar = m.mat.CreateMatrix()
        #print(f"m.mat.nze = {m.mat.nze}, a.mat.nze={a.mat.nze}, mstar.nze={mstar.nze}")
        #print(f"mstar.nze={mstar.nze}, len(mstar.AsVector())={len(mstar.AsVector())}")

        inva=a.mat.Inverse(freedofs=fes.FreeDofs())


        f = LinearForm(fes)
        f += sigma*(prevB-oldB)*v*dx('core')#
        f.Assemble()

        #res = m.mat * gfuPrev.vec - a.mat *gfuD.vec + jacmat.mat * oldgfu.vec
        res = f.vec - a.mat *gfuD.vec + jacmat.mat * oldgfu.vec
        gfu.vec.data = gfuD.vec + inva * res

        #NR pogreska:
        errfunc = (gfu - oldgfu).Norm()/(oldgfu.Norm()+1e-15) #/oldApot
        
        defon = mesh.Materials('core')
        error=Integrate(errfunc.Norm(), mesh, definedon=defon)
        #print('error =', error)
        errorlist.append(error)

        #oldgfu.vec.data= p*gfu.vec + (1-p)*old.vec  #damping
        oldgfu.vec.data= gfu.vec
        
        #Habs= oldgfu.Norm()
    
    print('errorlist', errorlist)

    #DgfuDt.vec.data = (gfu.vec - gfuPrev.vec)
    gfuPrev.vec.data=gfu.vec
    #HabsPrev= gfuPrev.Norm()
    prevB.vec.data=oldB.vec

    prevHmatrix= np.append(prevHmatrix, [gfuPrev.vec.FV().NumPy()], axis=0)

    
    Hlist.append(gfu[0](mesh(0,0.495*brid)))
    #Blist.append((mu*gfu)[0](mesh(0,0.495*brid)))
    Jlistoftuple.append((rot*grad(gfu))(mesh(0,0.495*brid)))
    #print(DgfuDt.Norm()(mesh(0,0.49*brid)))
    time_axis.append(time)

    #current power:
    #grad(gfu)*grad(gfu)==rot*grad(gfu)*rot*grad(gfu)
    currPow=Integrate(sigma*grad(gfu)*grad(gfu), mesh, 
                    definedon=mesh.Materials('core'))
    currBavg=Integrate(oldB/brid**2, mesh, definedon=mesh.Materials('core'))
    Plist.append(currPow)
    Bavglist.append(currBavg)
    
    gfut.AddMultiDimComponent(gfu.vec)

    cnt += 1; time = cnt * dt


####....... P O S T P R O C E S I N G........
##

#Blist=[elem[0] for elem in Blistoftuple]
Jlist=[elem[0] for elem in Jlistoftuple] #elem[0] je Jx komponenta

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
plt.plot(time_axis[pnt1:],Hlist[pnt1:],'b')
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

H= gfu
Draw(H, mesh, f"H_{comp}")
J=rot*grad(gfu) *sig
Draw(J, mesh, f"J_{comp}")


