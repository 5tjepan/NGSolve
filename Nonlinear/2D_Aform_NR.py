


with open("M140.dat", "r") as file:
    # Citanje linija iz datoteke
    lines = file.readlines()

    # Inicijalzacija praznih lista za stupce
    stupac1 = []
    stupac2 = []

    # Iteriranje krooz svaku liniju u datoteci
    for line in lines:
        # Razdvajanje vrijednosti u svakoj liniji koristeći tabulator kao separator
        vrijednosti = line.strip().split("\t")

        # Dodavanje vrijednosti u odgovarajuce stupce
        stupac1.append(vrijednosti[0])
        stupac2.append(vrijednosti[1])

#pretvaranje liste stringova u listu floatova
#Bref=list(map(float, stupac1))
#Href=list(map(float, stupac2))


from ngsolve import *
from netgen.occ import *
import netgen.gui
import matplotlib.pyplot as plt
import sys
sys.argv = ["fun"]

outer= Circle((0,0), 0.2).Face()
outer.edges.name = 'rub'

inner = MoveTo(-0.08,-0.08).Line(0.16,0.0).Line(0,0.16).Line(-0.16,0).Close().Face()
#inner = MoveTo(-0.8,-0.8).Line(1.6,0).Line(0,1.6).Line(-1.6,0).Close().Face()
#inner = MoveTo(-0.8,-0.8).Line(1.6,-0.1).Line(0,1.7).Line(-1.6,-0.2).Close().Face()
#inner = MoveTo(0.0,-1.0).Line(1,1).Line(-1,1).Line(-1,-1).Close().Face()

inner.edges.name="interface"
inner.faces.maxh=0.01
outer = outer - inner

inner.faces.name="inner"
inner.faces.col = (1, 1, 0)  #colour
outer.faces.name="outer"

geo = Glue([inner, outer])
#Draw(OCCGeometry(geo));

mesh = Mesh(OCCGeometry(geo, dim=2).GenerateMesh(maxh=0.05, quad_dominated=False))
print('mesh.nv',mesh.nv)
print('mesh.nedge', mesh.nedge)
print('mesh.ne', mesh.ne)

########## VIZUALIZACIJA
""" mu0=4*pi*1e-7
H_KL_ref=[0,42,53,62,70,79,88,100,113,132,157,193,255,376,677,1624,1e9]
B_KL_ref=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1e9*mu0]
Hvis = Parameter(0) #klasa Parameter je dijete klase CF. Uvrstavanjem Hvis u pozivu BSpline() je poziv oblika BSpline(CF)...jer je __call__ overloadan

munonlin = BSpline(2, [0]+list(H_KL_ref), list(B_KL_ref)) #ovo bi trebala biti instanca klase BSpline
muder= munonlin.Differentiate() #munonlin.Differentiate() metoda daje objekt klase BSpline

munonlinvis = munonlin((1e-6+sqrt(Hvis*Hvis+1e-6)))#/(1e-6+sqrt(Hvis*Hvis+1e-6)) #prema __call__ poziv BSplinee(CF) daje CF.
mudervis= muder((1e-6+sqrt(Hvis*Hvis+1e-6))) #buduci da je muder objekt klase BSpline, prema __call__ poziv BSplinee(CF) daje CF

Bspl=[]
Bder=[]
Hvec=[i for i in range(2000)]
for k in Hvec:
    Hvis.Set(k)
    Bspl.append(munonlinvis(mesh())) #buduci da je munonlinvis CF funkcija, prema __call__ klase CF za dobivanje vrijednostii CF moze se koristiti input "mesh()"
    Bder.append(100*mudervis(mesh())) 
    

plt.xlabel('H')
plt.plot(H_KL_ref[:-1],B_KL_ref[:-1], linestyle='', marker='o')
plt.plot(Hvec, Bspl)
plt.plot(Hvec, Bder)
plt.show() """
##################



fes = HCurl(mesh, order=0, dirichlet="rub",  complex=True, nograds = False) #CMPLX
u, v = fes.TnT()

bb=CF((0.1*y,-0.1*x))
#bb=CF((10,10))

gfu = GridFunction(fes)
old = GridFunction(fes)
#dirich = GridFunction(fes)

#gfu.Set(bb, VOL_or_BND = BND)
gfu.Set(bb, definedon=mesh.Boundaries('rub'))
old.Set(bb, definedon=mesh.Boundaries('rub'))
#dirich.Set(bb, definedon=mesh.Boundaries('rub'))

omega=314
#rel = 1200
sigma = 2e2
sig = mesh.MaterialCF({ "inner" : 1 }, default=None)
#sigma=CoefficientFunction( (10, 0,  0, 10), dims=(2,2) )

mu0=4*pi*1e-7
Href=[0,42,53,62,70,79,88,100,113,132,157,193,255,376,677,1624,1e9]
Bref=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1e9*mu0]
nu=[Href[i]/Bref[i] if Bref[i]!=0 else 0 for i in range(len(Bref))]
#HBcurve = BSpline(4, [0,0,0]+list(Href), list(Bref)) #ovo bi trebala biti instanca klase BSpline
#diffHBcurve= HBcurve.Differentiate() #munonlin.Differentiate() metoda daje objekt klase BSpline

B=curl(old)
Babs= B.Norm() 
errorlist=[]
p=1

for i in range(1,81):
    print(f"####iteration i={i}")
    
    print('Babs=',Babs(mesh(0,0)))
    print('B =',B(mesh(0,0)))
    
    gfu.Set(bb, definedon=mesh.Boundaries('rub'))

    #Babs += 1e-7
    rel= Babs*100 + 1e-5
    dHdB= 2*Babs*100 + 1e-5
    #rel = HBcurve(Babs)  #Babs+1e-6
    #dHdB = diffHBcurve(Babs)

    #print('rel=', rel(mesh(0.6,0.6)))
    #print('dHdB=', dHdB(mesh(0.6,0.6)))

    term1 = (1/mu0)*curl(u)*curl(v)*dx('outer') + rel*curl(u)*curl(v)*dx('inner') + \
        1j*omega*sigma*u*v*dx('inner') + 1j*1e-4*u*v*dx('outer')
    #jac= (dHdB - rel)/Babs *(B*curl(v))/Babs * (Conj(B)*curl(u)) *dx('inner')
    jac= (dHdB - rel)*curl(u)*curl(v)*dx('inner')

    a = BilinearForm(term1+jac)
    a.Assemble()

    jacmat= BilinearForm(jac)
    jacmat.Assemble()

    dummy=CF((1e-17,1e-17))
    #rhs1 = dummy*v*dx
    #rhs1 = (dHdB - rel)*(Conj(B)*curl(old))*(B*curl(v))/(Babs*Babs)*dx('inner')
    rhs1 = (dHdB - rel)*curl(old)*curl(v) *dx('inner')
    f = LinearForm(rhs1)
    f.Assemble()

    gfulist=[v for v in gfu.vec]

    dirich=[]
    
    for i in range(len(fes.FreeDofs())):
        if fes.FreeDofs()[i]: dirich.append(gfu.vec[i])
    print('before dirich = ', dirich[:3]) 
    #if i<=5: print('before gfu.vec=', gfulist[:])

#    oldlist=[v for v in old.vec]
#    if i<=5: print('old.vec=', oldlist[:4])

    ##### SOLVER
    #solvers.BVP(bf=a, lf=f, gf=gfu, pre=None, maxsteps=2000, print=True)
    r = jacmat.mat * old.vec - a.mat * gfu.vec
    gfu.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs())*r

    dirich=[]
    for i in range(len(fes.FreeDofs())):
        if fes.FreeDofs()[i]: dirich.append(gfu.vec[i])
    print(' after dirich = ', dirich[:3]) 


    #errfunc = (gfu - old)/gfu
    errfunc = (gfu - old)/old
    defon = mesh.Materials('inner')
    error=Integrate(errfunc.Norm(), mesh, definedon=defon)
    print('error =', error)
    errorlist.append(error)

    old.vec.data= p*gfu.vec + (1-p)*old.vec
    
    B=curl(old)
    Babs= B.Norm()

print('errorlist', errorlist)
#####POSTPROCESING
A=gfu
E = - 1j * omega * A
J = - 1j * omega * sigma * A *sig
B = curl(gfu)

Pow=0.5*E*Conj(J) 
Peddy=Integrate(Pow, mesh, order=5)
print('Peddy',Peddy)

Draw(A, mesh, "A")
Draw (B, mesh, "B")
Draw (J, mesh, "J")


