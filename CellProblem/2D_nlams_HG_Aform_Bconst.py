
from ngsolve import *
from netgen.occ import *

#air= Circle((0,0), 0.008).Face()

size= 0.003 #0.001055 #0.000755 #0.00155 #0.003 #0.0026
n=9
d=0.001
ins=0.0001
fill_fact=n*d/(n*d+(n-1)*ins)
w=(n*d+(n-1)*ins)/2
b=w

air = MoveTo(-1.2*w, -1.2*b).Rectangle(2.4*w,2.4*b).Face()
air.edges.name = 'rub'

lamele = MoveTo(-w,-b).Rectangle(2*w,2*b).Face()
lamele.faces.name="lamele"
#lamele.faces.edges.name="lamele"
lamele.faces.maxh=size
air -=lamele
air.faces.name="air"
air.faces.col = (1, 1, 0)

geo = Glue([air, lamele])


mesh = Mesh(OCCGeometry(geo, dim=2).GenerateMesh(maxh=size, quad_dominated=True))

fes = HCurl(mesh, order=0, dirichlet="rub",  complex=True, nograds = False) #CMPLX
gfu = GridFunction(fes)
""" gfu.Load('my_sol.sol')
Draw(curl(gfu),mesh, 'gfu')
print(nvghch) """

print('mesh.nv',mesh.nv)
print('mesh.nedge', mesh.nedge)
print('mesh.ne', mesh.ne)

u, v = fes.TnT()

#bb=CF((0.1*y,-0.1*x))
B0=0.5
bb=CF((0*y,-B0*x/1.2**2))
gfu.Set(bb, definedon=mesh.Boundaries('rub'))

omega=314 * 1
mu0 = 1.257e-6
mu = 1000*mu0 *5



#ovo je za mu=5000mu0 i za f=50Hz ::: za N = 18 FEs
""" rel_cc= 176.9375561833921 + 56.528939111002565*1j
Pcoeff_cc = 8875.36
rel_cf= 176.8737766139486+45.42691850894148*1j
Pcoeff_cf = 7135.90265065109
rel = IfPos((y-2*b/3)*(y+2*b/3), rel_cf, rel_cc)
Pcoeff = IfPos((y-2*b/3)*(y+2*b/3), Pcoeff_cf, Pcoeff_cc)
#... """

#ovo je za mu=5000mu0 i za f=50Hz ::: stari za N = 9 FEs
""" rel_cc= 176.9375561833921 + 56.528939111002565*1j
Pcoeff_cc = 8875.36
rel_cf= 177.0647090736546 + 50.91894586044948*1j 
Pcoeff_cf = 7994.56
rel = IfPos((y-b/3)*(y+b/3), rel_cf, rel_cc)
Pcoeff = IfPos((y-b/3)*(y+b/3), Pcoeff_cf, Pcoeff_cc) """
#...

""" #ovo je za mu=5000mu0 i za f = 500 Hz
rel_cc= 373.941756964345+390.9513824417485*1j
Pcoeff_cc = 613840.768214015
rel_cf= 352.3232606912534+337.09896842871456*1j
Pcoeff_cf = 529283.6429409353
rel = IfPos((y-b/3)*(y+b/3), rel_cf, rel_cc)
Pcoeff = IfPos((y-b/3)*(y+b/3), Pcoeff_cf, Pcoeff_cc) """
#...


import numpy as np

Kcenter= np.loadtxt('K4_data.txt')
Kfrin= np.loadtxt('K1_data.txt')

freq = Kcenter[:, 0]
Re_nuc = Kcenter[:, 2]
Im_nuc = Kcenter[:, 3]
CPc = Kcenter[:, 4]

Re_nuf = Kfrin[:, 2]
Im_nuf = Kfrin[:, 3]
CPf = Kfrin[:, 4]

ind=2 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ovo je za mu=5000mu0 i za f = 500 Hz
rel_cc= Re_nuc[ind] + Im_nuc[ind]*1j
Pcoeff_cc = CPc[ind]
rel_cf= Re_nuf[ind] + Im_nuf[ind]*1j
Pcoeff_cf = CPf[ind]
rel = IfPos((y-b/3)*(y+b/3), rel_cf, rel_cc)
Pcoeff = IfPos((y-b/3)*(y+b/3), Pcoeff_cf, Pcoeff_cc)

print('f=',freq[ind],'Hz')


Draw(rel,mesh,'rel')

sigma = 2e6
sig = IfPos( (x+b)*(x-b), 0, 1)
Draw(sig, mesh, 'sig')

a = BilinearForm(fes)
a += 10*(1/mu0)*curl(u)*curl(v)*dx('air') + rel*curl(u)*curl(v)*dx('lamele')+ 1j*1e-4*u*v*dx 
a.Assemble()

#Jg=CF((0,8000*(x-2*w)/mu*sig)) 
#Jg=CF((0,50*(400*(x-w))/mu*sig)) 
Jg=CF((0,-50/mu*sig)) #*0
f = LinearForm(fes)
f += Jg*v*dx
f.Assemble()
Draw(Jg,mesh,'Jg')

r = f.vec - a.mat * gfu.vec
gfu.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs())*r

A=gfu
B = curl(gfu)

core = mesh.Materials('lamele')
#Pow=Pcoeff*(B*Conj(B)).Norm() #ili (B.Norm())**2 !!!
Pow=Pcoeff*(B.Norm())**2
Peddy= Integrate(Pow, mesh, definedon=core)
print('Peddy',Peddy)

area= Integrate(CF(1), mesh, definedon=core)
print('area= ',area)

tok= Integrate(B, mesh, definedon=core)
print('tok= ',tok)

print('Bavg= ',tok/area)


Draw(A, mesh, "A")
Draw (B, mesh, "B")
Draw (Pow, mesh, "Pow")


#gfu.Save('my_sol.txt')

""" import numpy as np
import matplotlib.pyplot as plt
import sys
sys.argv = ["dummy"]

X = np.linspace(-0.5*d, 0.5*d, num=21)
Y = np.zeros_like(X)+0.00065

plt.plot(X, (B).imag(mesh(X, Y)))
plt.xlabel('x')
plt.show() """


postfes = H1(mesh, order=1, definedon='lamele', complex=True) #dirichlet=Svi_dirich
Blin = postfes.TrialFunction()
beta = postfes.TestFunction()


posta = BilinearForm( ( Blin*beta)*dx).Assemble()
postf = LinearForm(B*beta*dx).Assemble()

Bset = GridFunction(postfes)
Bpost = GridFunction(postfes)
Bset.Set(B)
Bpost.vec.data = posta.mat.Inverse(postfes.FreeDofs()) * postf.vec

#c = Preconditioner(posta, type="direct", inverse="umfpack")
#solvers.BVP(posta,postf, gfpost, c)

Draw(Bpost, mesh, 'Bpost')
Draw(Bset, mesh, 'Bset')


#print('integral = ', Integrate(Bpost, mesh))

Powpost=Pcoeff*(Bpost.Norm())**2
Peddypost= Integrate(Powpost, mesh, definedon=core)
print('Peddyppost',Peddypost)

Powset=Pcoeff*(Bset.Norm())**2
Peddyset= Integrate(Powset, mesh, definedon=core)
print('Peddyset',Peddyset)

import matplotlib.pyplot as plt
import sys
sys.argv = ["dummy"]

X = np.linspace(-w, w, num=51)
Y = np.zeros_like(X)+ w*0.0

plt.plot(X, (Bpost).imag(mesh(X, Y)), color='green')
plt.plot(X, (Bset).imag(mesh(X, Y)), color='red')
plt.plot(X, (B).imag(mesh(X, Y)))
plt.xlabel('x')
plt.show()


