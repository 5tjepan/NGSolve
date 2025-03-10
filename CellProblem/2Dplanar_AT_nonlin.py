
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
mu0=4*pi*1e-7

B_ref=[0.001, 0.04, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.74, 1.77]
H_ref=[0.001, 3.5, 5.0, 7.6, 10.0, 12.07, 14.08, 16.0, 17.75, 19.42, 21.05, 22.64, 24.16, 25.76, 27.98, 29.66, 32.05, 35.33, 40.06, 47.17, 58.35, 77.19, 107.0, 160.0, 270.0]

Bvis = Parameter(0) #klasa Parameter je dijete klase CF. Uvrstavanjem Hvis u pozivu BSpline() je poziv oblika BSpline(CF)...jer je __call__ overloadan

munonlin = BSpline(2, [0]+list(B_ref), list(H_ref)) #ovo bi trebala biti instanca klase BSpline
muder= munonlin.Differentiate() #munonlin.Differentiate() metoda daje objekt klase BSpline

munonlinvis = munonlin((1e-6+sqrt(Bvis*Bvis+1e-6)))#/(1e-6+sqrt(Hvis*Hvis+1e-6)) #prema __call__ poziv BSplinee(CF) daje CF.
mudervis= muder((1e-6+sqrt(Bvis*Bvis+1e-6))) #buduci da je muder objekt klase BSpline, prema __call__ poziv BSplinee(CF) daje CF

Hvec=[]
Hder=[]
Bvec=[i*1.75/100 for i in range(1,102)]
for k in Bvec:
    Bvis.Set(k)
    Hvec.append(munonlinvis(mesh())) #buduci da je munonlinvis CF funkcija, prema __call__ klase CF za dobivanje vrijednostii CF moze se koristiti input "mesh()"
    Hder.append(mudervis(mesh())) 
    

#plt.xlabel('B')
#plt.scatter(B_ref[:-1],H_ref[:-1])
#plt.plot(Bvec, Hvec)
#plt.plot(Bvec, Hder)
#plt.show()
##################

fsU = HCurl(mesh, order=0, dirichlet="rub",  complex=True, nograds = False) #CMPLX
fsV = H1(mesh, order=1, dirichlet='interface', definedon='inner', complex=True)
fes=fsU*fsV
mvp, csp = fes.TrialFunction()
alpha, theta = fes.TestFunction()

bb=CF((0.01*y,-0.01*x))

gfu = GridFunction(fes)
old = GridFunction(fes)

Apot, Tpot = gfu.components
oldApot, oldTpot = old.components


Apot.Set(bb, definedon=mesh.Boundaries('rub'))
oldApot.Set(bb, definedon=mesh.Boundaries('rub'))

omega=314
sigma = 2e2
rho=1/sigma
sig = mesh.MaterialCF({ "inner" : 1 }, default=None)

HBcurve = BSpline(2, [0]+list(B_ref), list(H_ref))
diffHBcurve= HBcurve.Differentiate() #munonlin.Differentiate() metoda daje objekt klase BSpline

B=curl(oldApot)
Babs= B.Norm() 
errorlist=[]
p=1

for i in range(1,20):
    print(f"####iteration i={i}")
    
    print('Babs=',Babs(mesh(0,0)))
    print('B =',B(mesh(0,0)))
    
    Apot.Set(bb, definedon=mesh.Boundaries('rub'))
    #Apot.Set(bb, definedon=mesh.Boundaries('rub|interface'))
    #oldApot.Set(bb, definedon=mesh.Boundaries('rub'))

    #Babs += 1e-7
    #rel= Babs*100 + 1e-5
    #dHdB= 2*Babs*100 + 1e-5

    #rel=munonlin(Babs+1e-6)
    #dHdB=muder(Babs+1e-6)

    rel= 300.001 * Babs
    dHdB=300.0 *2*Babs

    #rel = HBcurve(Babs+1e-6)  #Babs+1e-6
    #dHdB = diffHBcurve(Babs+1e-6)

    #print('rel=', rel(mesh(0.6,0.6)))
    #print('dHdB=', dHdB(mesh(0.6,0.6)))

    #term1 = (1/mu0)*curl(alpha)*curl(mvp)*dx('outer') + rel*curl(alpha)*curl(mvp)*dx('inner') + \
    # 1j*omega*sigma*alpha*mvp*dx('inner') + 1j*1e-4*alpha*mvp*dx('outer')
    
    term1 = rho * grad(csp)*grad(theta)*dx('inner') + 1j * omega*curl(mvp)*theta*dx('inner') + 1j*0.01*mvp*alpha*dx
    term2 = - 1j*omega/mu0*curl(mvp)*curl(alpha)*dx('outer') - 1j*omega* rel*curl(mvp)*curl(alpha)*dx('inner') + 1j*omega* csp*curl(alpha)*dx('inner') 
    jac= -1j*omega*(dHdB - rel)*curl(alpha)*curl(mvp)*dx('inner')

    a = BilinearForm(term1+term2 + jac)
    a.Assemble()

    jacmat= BilinearForm(jac)
    jacmat.Assemble()

    dummy=CF((1e-17,1e-17))
    #rhs1 = dummy*v*dx
    #rhs1 = (dHdB - rel)*(Conj(B)*curl(old))*(B*curl(v))/(Babs*Babs)*dx('inner')
    rhs1 = (dHdB - rel)*curl(oldApot)*curl(alpha) *dx('inner')
    f = LinearForm(rhs1)
    f.Assemble()


    dirich=[]
    
    for i in range(len(fes.FreeDofs())):
        if fes.FreeDofs()[i]: dirich.append(gfu.vec[i])
    print('before dirich = ', dirich[:3]) 
    #if i<=5: print('before gfu.vec=', gfulist[:])

    ##### SOLVER
    #solvers.BVP(bf=a, lf=f, gf=gfu, pre=None, maxsteps=2000, print=True)
    
    r = jacmat.mat * old.vec - a.mat * gfu.vec
    #r = f.vec - a.mat * gfu.vec
    gfu.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs())*r

    dirich=[]
    for i in range(len(fes.FreeDofs())):
        if fes.FreeDofs()[i]: dirich.append(gfu.vec[i])
    print(' after dirich = ', dirich[:3]) 

    #errfunc = (gfu - old)/old
    errfunc = (Apot - oldApot).Norm()/oldApot.Norm()
    defon = mesh.Materials('inner')
    error=Integrate(errfunc.Norm(), mesh, definedon=defon)
    print('error =', error)
    errorlist.append(error)

    old.vec.data= gfu.vec
    
    B=curl(oldApot)
    Babs= B.Norm()

print('errorlist', errorlist)
#####POSTPROCESING
potA=oldApot
rot=CF( (0 , 1,  -1, 0), dims=(2,2) )
J = rot*grad(Tpot)
B = curl(potA)

Pow=0.5*rho*J*Conj(J)
Peddy=Integrate(Pow, mesh, order=5)
print('Peddy',Peddy)

Draw(potA, mesh, "A")
Draw (B, mesh, "B")
Draw (J, mesh, "J")

#print('Bref=',Bref)
#print('Href=',Href)
