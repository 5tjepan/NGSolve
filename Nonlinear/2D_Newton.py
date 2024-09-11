
with open("M140.dat", "r") as file:
    # Citanje linija iz datoteke
    lines = file.readlines()

    # Inicijalizacija praznih lista za stupce
    stupac1 = []
    stupac2 = []

    # Iteriranje kroz svaku liniju u datoteci
    for line in lines:
        # Razdvajanje vrijednosti u svakoj liniji koristeći tabulator kao separator
        vrijednosti = line.strip().split("\t")

        # Dodavanje vrijednosti u odgovarajuce stupce
        stupac1.append(vrijednosti[0])
        stupac2.append(vrijednosti[1])

#pretvaranje liste stringova u listu floatova
Bref=list(map(float, stupac1))
Href=list(map(float, stupac2))


from ngsolve import *
from netgen.occ import *

import matplotlib.pyplot as plt
import sys
sys.argv = ["fun"]

outer= Circle((0,0), 0.02).Face()
outer.edges.name = 'rub'

inner = MoveTo(-0.008,-0.008).Line(0.016,0.0).Line(0,0.016).Line(-0.016,0).Close().Face()
inner.edges.name="interface"
inner.faces.maxh=0.002
outer = outer - inner

inner.faces.name="inner"
inner.faces.col = (1, 1, 0)  #colour
outer.faces.name="outer"

geo = Glue([inner, outer])
#Draw(OCCGeometry(geo));

mesh = Mesh(OCCGeometry(geo, dim=2).GenerateMesh(maxh=0.1, quad_dominated=False))
print('mesh.nv',mesh.nv)
print('mesh.nedge', mesh.nedge)
print('mesh.ne', mesh.ne)


fes = HCurl(mesh, order=0, dirichlet="rub", complex=True, nograds = False) #CMPLX
u, v = fes.TnT()

bb=CF((0.1*y,-0.1*x))

gfu = GridFunction(fes)
old = GridFunction(fes)


#gfu.Set(bb, VOL_or_BND = BND)
gfu.Set(bb, definedon=mesh.Boundaries('rub'))
#old.Set(bb, definedon=mesh.Boundaries('rub'))

omega=10
#rel = 1200
sigma = 2e2
sig = mesh.MaterialCF({ "inner" : 1 }, default=None)
#sigma=CoefficientFunction( (10, 0,  0, 10), dims=(2,2) )

mu0=4*pi*1e-7
#Href=[0,42,53,62,70,79,88,100,113,132,157,193,255,376,677,1624]
#Bref=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5]
#HBcurve = BSpline(2, [0]+list(Href), list(Bref)) #ovo bi trebala biti instanca klase BSpline
#diffHBcurve= HBcurve.Differentiate() #munonlin.Differentiate() metoda daje objekt klase BSpline

B=curl(old)
Babs= B.Norm()

for i in range(1,59):
    print(f"######iteration i={i}")
    
    print('Babs(0,0)=',Babs(mesh(0,0)))
    print('B(0,0) =',B(mesh(0,0)))
    
    gfu.Set(bb, definedon=mesh.Boundaries('rub'))

    #B+= 1e-6
    Babs += 1e-6
    rel= Babs*Babs
    dHdB=3*Babs*Babs
    #rel = HBcurve(Babs) 
    #dHdB = diffHBcurve(Babs)
    #print('rel=', rel(mesh(0.6,0.6)))
    #print('dHdB=', dHdB(mesh(0.6,0.6)))

    term1 = (1/mu0)*curl(u)*curl(v)*dx('outer') + rel*curl(u)*curl(v)*dx('inner') + \
        1*u*v*dx ('outer') +1j*omega*sigma*u*v*dx('inner')
    jac= (dHdB - rel)*(B*curl(u))*(B*curl(v))/(Babs*Babs)*dx('inner')

    a = BilinearForm(term1 + jac)
    a.Assemble()

    jacmat= BilinearForm(jac)
    jacmat.Assemble()

    dummy=CF((1e-17,1e-17))
    rhs1 = (dHdB - rel)*(B*curl(old))*(B*curl(v))/(Babs*Babs)*dx('inner')
    f = LinearForm(rhs1)
    f.Assemble()

    gfulist=[v for v in gfu.vec]
    oldlist=[vv for vv in old.vec]
    if i<=3: print('old.vec=', oldlist[0:3])
    if i<=3: print('gfu.vec=', gfulist[0:3])

    ##### SOLVER
    solvers.BVP(bf=a, lf=f, gf=gfu, pre='direct', maxsteps=200, print=True)
    #r = jacmat.mat * old.vec - a.mat * gfu.vec
    #gfu.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs())*r

    errfunc = (gfu - old)/old
    error=Integrate(errfunc.Norm(), mesh)
    print('error =', error)

    old.vec.data=gfu.vec.data
    
    B=curl(old)
    Babs= B.Norm()


#####POSTPROCESING
A=gfu
B = curl(gfu)
E = - 1j * omega * A
J = - 1j * omega * sigma * A *sig

Draw(A, mesh, "A")
Draw (B, mesh, "B")
Draw (J, mesh, "J")
