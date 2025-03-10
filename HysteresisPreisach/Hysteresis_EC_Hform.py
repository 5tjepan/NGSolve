from ngsolve import *
from netgen.occ import *
import numpy as np
#from extremaseek import FindDominantExtrema
#from everett import Everett_atan,Everett_exp, Everett_exp2
#from preisach import preisach_output

#-----------------------


def localmaxima(u):
    if len(u)<3: 
        lex=np.arange(len(u))
    else: #maksimum je ako je >= od prethodnog i iskljucivo veci od iduceg
        lex= np.nonzero((u[1:-1]>=u[0:-2]) & (u[1:-1]>u[2:]))[0] + 1 #np.nonzero() vraca tuple of arrays pa zato [0]
    
        if (u[0]>u[1]) : lex = np.concatenate(([0], lex))
        n=len(u)
        if (u[-1]>=u[-2]): lex = np.concatenate((lex, [n-1]))
    return lex


#------------------------

#pronalazi niz dominantnih ekstrema
def dominant_extrema(u, ekstremi):
    #pretpostavljamo da je len(ekstremi)>2, tj imamo barem tri lokalna ekstrema
    domEkstremi=[]  #inicijalizacija liste za dominantne ekstreme
    domEkstremi.append(ekstremi[-1])  #posljednji lokalni ekstrem je uvijek i dominantni ekstrem
    #if (u[ekstremi[-1]] == u_max): return domEkstremi

    globExtrInd=np.nonzero(np.abs(u[ekstremi])==np.max(np.abs(u[ekstremi])))[0][-1]
    ekstremi=ekstremi[globExtrInd:]
    
    if (len(ekstremi)==1): return domEkstremi

    flag = u[ekstremi[-2]] < u[domEkstremi[-1]]  #flag=True if domEkstrema[-1] je lokalni maksimum i obratno
    for i in range(len(ekstremi)-2,-1,-1):
        if (flag and (u[ekstremi[i]] > u[domEkstremi[-1]])):
            #domEkstremi.append(ekstremi[i] + np.argmin(u[ekstremi[i]:domEkstremi[-1]])) #daje krivi index kod konst funkcije, al tocan iznos ekstrema
            domindex=ekstremi[i] + np.nonzero(u[ekstremi[i]:domEkstremi[-1]]==np.min(u[ekstremi[i]:domEkstremi[-1]]))[0][-1]
            domEkstremi.append(domindex)
            flag=False
        elif (not flag and (u[ekstremi[i]] < u[domEkstremi[-1]])):
            #domEkstremi.append(ekstremi[i] + np.argmax(u[ekstremi[i]:domEkstremi[-1]]))
            domindex=ekstremi[i] + np.nonzero(u[ekstremi[i]:domEkstremi[-1]]==np.max(u[ekstremi[i]:domEkstremi[-1]]))[0][-1]
            domEkstremi.append(domindex)
            flag=True
    #...ova petlja ne moze appendat najstariji dominantni maksimum/minimum M0 pa moramo jos njega pronaci i provjeriti je li on globalni ekstrem
    #Cim smo u ovoj funkciji znamo da je len(ekstremi)>2, pa zbog prethodne petlje znamo da postoji \
    #...dominantni maksimum/minimum M0 stariji od u[domEkstremi[-1]] al ne znamo je li on ujedno i dominantni ekstrem jer svojstvo najstarijeg
    #dominantnog ekstrema je to da je on ujedno i globalni ekstrem!
    #Dakle, da bi M0 appendali u domEkstremi, mora vrijediti da je abs(M0)>abs(u[domEkstremi])
    #ono sto sigurno znamo je da postoji ekstrem E0 takav da je abs(E0) > abs(domEkstremi[-2]), ako postoji abs(domEkstremi[-2])
    #tj. ako je len(domEkstremi)>1
    #domEkstremi[-1] je trenutno vremenski najstariji dominantni ekstrem
    #ako je domEkstremi[-1] maksimum, moguce je da prebrise najstariji dominantni ekstrem (minimum) i obratno.
    #...to je moguce u slucaju da je np.abs(u[domEkstremi[-1]]) > od apsolutne vrijendosti najstarijeg dominantnog ekstrema

    domEkstremi.append(ekstremi[0])

    domEkstremi.reverse()
    return domEkstremi


#---------------------------

#sortiranje minimuma i maksimuma u listu ekstrema::
def localextrema(u):
    ekstremi=[]
    if len(u)<3: 
        ekstremi = localmaxima(u) #jer u ima samo dva clana
    #    print('ekstremi=',ekstremi)
        return ekstremi
    
    max_indeksi = localmaxima(u)
    min_indeksi = localmaxima(-u)
    #print('max_indeksi=', max_indeksi)
    #print('min_indeksi=', min_indeksi)
    if (min_indeksi[0]<=max_indeksi[0]):
        for i in range(len(max_indeksi)+len(min_indeksi)):
            if i%2==0:
                ekstremi.append(min_indeksi[i//2])
            else: 
                ekstremi.append(max_indeksi[i//2])
    else:
        for i in range(len(max_indeksi)+len(min_indeksi)):
            if i%2==0:
                ekstremi.append(max_indeksi[i//2])
            else: 
                ekstremi.append(min_indeksi[i//2])
    #print('ekstremi=',ekstremi)
    return ekstremi


#----------------------------

def FindDominantExtrema(u,u_max):

    ekstremi=localextrema(u)
    
    if abs(u[ekstremi[-1]])==abs(u_max):
        domindex=[ekstremi[-1]] #lista s jednim clanom
    elif len(ekstremi)==1:
        domindex = ekstremi
    else:
        domindex = dominant_extrema(u,ekstremi)
    #print('domindex=', domindex)

    return domindex


#((((((((((((((((((((((((((((((((((((()))))))))))))))))))))))))))))))))))


def Everett_exp(a,b):
    c=0.25
    y= 1/(c + np.exp(-0.03*a) + np.exp(+0.03*b) ) - \
          1/(c + np.exp(+0.03*a) + np.exp(-0.03*b) )
    return y

def Everett_exp2(a,b):
    y= 1/(0.25 + np.exp(-0.03*a) + np.exp(+0.03*b) + np.exp(-0.04*a) + np.exp(+0.04*b)) - \
          1/(0.25 + np.exp(+0.03*a) + np.exp(-0.03*b) +np.exp(+0.04*a) + np.exp(-0.04*b))
    return y

def Everett_atan(x,y):
    a=0.0196483
    b=2.95329554
    c=0.02211744
    d=1.04359946
    
    valid= x>=y

    alpha=valid*x
    beta=valid*y

    #E= (np.arctan(a*x) - np.arctan(a*y))**b + (np.arctan(c*x)**3 - np.arctan(c*y)**3)**d
    E= 0.1*((np.arctan(a*alpha) - np.arctan(a*beta))**b + (np.arctan(c*alpha)**3 - np.arctan(c*beta)**3)**d)
    return E


#[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]


#---------------
def preisach_output(u, domindex, u_max, Everett_fun):
    domex=u[domindex]

    Everett=Everett_fun
    if np.abs(domex[-1]) == u_max: 
        output= 0.5* np.sign(domex[-1]) * Everett(u_max,-u_max)
        return output

    output= np.sign(domex[0]) * 0.5*Everett(np.abs(domex[0]),-np.abs(domex[0]))
    #print('domex=',domex[0])
    #print('----->', output)
    if len(domex)>1: 
        flag=domex[1]>domex[0] #domex[0] je globalni ekstrem i zanima nas je li globalni min ili globalni max
        if len(domex)%2==0: domex=np.append(domex,domex[-1]) #appendanjem domex[-1] na domex ne mijenja se output jer je Everett(a,a)=0
        if flag:
            for i in range(1,len(domex),2):
                output+= Everett(domex[i],domex[i-1]) - Everett(domex[i], domex[i+1])
                #print('output2',output)
        else:
            for i in range(1,len(domex),2):
                output+= -Everett(domex[i-1],domex[i]) + Everett(domex[i+1], domex[i])
                #print('output3=',output)
    return output    
#-----------------
#
#
#
#
#
#|||||||||||||||||||||||||||
#||||  P  O  Z  I  V    ||||
#|||||||||||||||||||||||||||
#
#
#
#
#
outer= Circle((0,0), 0.15).Face()
outer.edges.name = 'rub'
d=0.0005
brid=d*5
core = MoveTo(-d/2,-brid/2).Line(d,0.0).Line(0,brid).Line(-d,0).Close().Face()
core.edges.name="rub"
#core.faces.maxh=0.005
outer = outer - core
core.faces.name="core"
core.faces.col = (1, 1, 0) 
outer.faces.name="outer"

geo = Glue([core, outer])
#Draw(OCCGeometry(geo));

mesh = Mesh(OCCGeometry(core, dim=2).GenerateMesh(maxh=d/12, quad_dominated=False))
print('mesh.nv',mesh.nv)
print('mesh.nedge', mesh.nedge)
print('mesh.ne', mesh.ne)

#.............
fes = H1(mesh, order=0, dirichlet="rub", complex=False)
H, v = fes.TnT()

omega=2*pi*50 
dt = 2*pi/omega /100
cnt = 0; 
t0=0.0
tend=1.22 *2*pi/omega
time = t0

mu0 = 1.257e-6
#rel = 200
sigma = mesh.MaterialCF({ "core" : 2e6 }, default=None)
sig = mesh.MaterialCF({ "core" : 1 }, default=None)
Kf=1
area= Integrate(CF(1), mesh)


#-------------
mu0=4*pi*1e-7

#B_ref=[0.001, 0.04, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.83, 1000]
#H_ref=[0.001, 3.5, 5.0, 7.6, 10.0, 12.07, 14.08, 16.0, 17.75, 19.42, 21.05, 22.64, 24.16, 25.76, 27.98, 29.66, 32.05, 35.33, 40.06, 47.17, 58.35, 77.19, 113.0, 196.0, 395.0, 800, 1000/mu0]

#BHcurve = BSpline(2, [0]+list(H_ref), list(B_ref)) #ovo bi trebala biti instanca klase BSpline
#diffBH= BHcurve.Differentiate() #HBcurve.Differentiate() metoda daje objekt klase BSpline
#-------------

gfu = GridFunction(fes)
oldgfu = GridFunction(fes)
t = Parameter(0.0)
#...
ic=CF(0) #initial condition
gfu.Set(ic)
oldgfu.Set(CF(0))
#..
bc = 80*sin(omega*t) #boundary condition
gfu.Set(bc, definedon=mesh.Boundaries('rub'))
#Draw(gfu, mesh, "gfu_bnd")

#inicijalizacija listi:::
time_axis=[]
Hlist=[]
Blist=[]
Jlistoftuple=[]
Plist=[]
Bavglist=[]
Hrublist=[]

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

Everett = Everett_atan
#............................
while time < tend - 0.5 * dt:  #Euler time-stepping 
#for iter in range(1,14):  #Euler time-stepping 

    t.Set(time)
    gfuD.Set(bc,definedon=mesh.Boundaries('rub')) #vrem. ovisan bc (t.Set(time))
    errorlist=[]
    
    #HabsPrev= gfuPrev.Norm()
    #Habs= oldgfu.Norm()
    
    for it in range(1, 10):  #NewtonRaphson

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

        #====================
        # h i s t e r e z a 
        Everett = Everett_atan
        Hvec = np.array([oldgfu.vec.FV().NumPy()])
        oldHmat=np.append(prevHmatrix, Hvec, axis=0)
        Hmax=100
        oldBlist=[]
        for el in range(gfu.vec.size):
            Helem=oldHmat[:,el]
            domIndx=FindDominantExtrema(Helem,Hmax)
            oldBlist.append(preisach_output(Helem,domIndx,Hmax, Everett_fun=Everett))

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

    
    Hlist.append(gfu[0](mesh(0,0.05*brid)))
    Blist.append(oldB[0](mesh(0,0.05*brid)))
    Jlistoftuple.append((rot*grad(gfu))(mesh(0,0.495*brid)))
    time_axis.append(time)
    Hrublist.append(bc(mesh(0,0)))

    #current power:
    #grad(gfu)*grad(gfu)==rot*grad(gfu)*rot*grad(gfu)
    currPow=Integrate(1/sigma*grad(gfu)*grad(gfu), mesh, 
                    definedon=mesh.Materials('core'))
    currBavg=Integrate(oldB/area, mesh, definedon=mesh.Materials('core'))
    Plist.append(currPow)
    Bavglist.append(currBavg)
    
    gfut.AddMultiDimComponent(gfu.vec)

    cnt += 1; time = cnt * dt


####....... P O S T P R O C E S I N G........
##

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

""" plt.subplot(411)
plt.plot(Hlist[pnt1:],Blist[pnt1:],'y')
plt.grid(True)
#... """

plt.subplot(511)
plt.plot(Hrublist[pnt1:],Bavglist[pnt1:],'y')  # BH-loop
plt.grid(True)
plt.xlabel("H_edge (T)")
plt.ylabel("B_avg (T)")
#...

plt.subplot(512)
plt.plot(time_axis[pnt1:],Hlist[pnt1:],'g')
plt.grid(True)
plt.ylabel("H_at_pnt1")
plt.xlabel("time")
#...

plt.subplot(513)
plt.plot(time_axis[pnt1:],Blist[pnt1:],'b')
plt.grid(True)
plt.ylabel("B_at_pnt1")
#...

plt.subplot(514)
plt.plot(time_axis[pnt1:], Jlist[pnt1:],'r--')
plt.grid(True)
plt.ylabel("J_at_pnt1")
#...

plt.subplot(515)
plt.plot(time_axis[pnt1:], Plist[pnt1:],'m')
plt.grid(True)
plt.ylabel("P (W/m)")
plt.show()

#...
comp=11
gfu.vec.data=gfut.vecs[comp]

H= gfu
Draw(H, mesh, f"H_{comp}")
J=rot*grad(gfu)
Draw(J, mesh, f"J_{comp}")

Draw(oldB,mesh, 'B')
#print(oldB.vec.data)

