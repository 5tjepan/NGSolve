


""" with open("M140.dat", "r") as file:
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
#Href=list(map(float, stupac2)) """


from ngsolve import *
from netgen.csg import *
import netgen.gui
import matplotlib.pyplot as plt
import sys
sys.argv = ["fun"]

def MakeGeometry():
    geometry = CSGeometry()
    box = OrthoBrick(Pnt(-1,-1,-1),Pnt(1,1,1)).bc("outer")

    front= Plane (Pnt(-0.05,-0.05,-0.2), Vec(-1,0,0) ).bc("front").maxh(0.015)
    right= Plane (Pnt(-0.05,-0.05,-0.2), Vec(0,-1,0) ).bc("right")
    bot  = Plane (Pnt(-0.05,-0.05,-0.2), Vec(0,0,-1) ).bc("bot")
    back = Plane (Pnt(0.05,0.05,0.2), Vec(1, 0,0) ).bc("back").maxh(0.015)
    left = Plane (Pnt(0.05,0.05,0.2), Vec(0,1,0) ).bc("left")
    top  = Plane (Pnt(0.05,0.05,0.2), Vec(0,0, 1) ).bc("top")
    core = left * right * front * back * bot * top
    #core = OrthoBrick(Pnt(-0.05,-0.05,-0.2),Pnt(0.05,0.05,0.2))
    core.maxh(0.06)
    #front.maxh(0.01)
    core.mat("core")
    
    coil = (Cylinder(Pnt(0.0,0,-0.5), Pnt(0.0,0,0.5), 0.2) - \
            Cylinder(Pnt(0.0,0,-0.5), Pnt(0.0,0,0.5), 0.1)) * \
            OrthoBrick (Pnt(-1,-1,-0.1),Pnt(1,1,0.1))
    coil.maxh(0.14)
    coil.mat("coil")
    
    geometry.Add ((box-core-coil).mat("air"))
    geometry.Add (core)
    geometry.Add (coil)

    return geometry


ngmesh = MakeGeometry().GenerateMesh(maxh=0.4)
#ngmesh.Save("coil.vol")
mesh = Mesh(ngmesh)
mesh.Curve(5)  

Draw(mesh)
print(mesh.GetBoundaries())

#----------------------


fsU = HCurl(mesh, order=1, dirichlet="outer", complex=True, nograds = False)
fsV = H1(mesh, order=2, dirichlet='left|right|bot|top', definedon='core', complex=True)
fes=fsU*fsV
mvp, csp = fes.TrialFunction()
alpha, tau = fes.TestFunction()

gfu = GridFunction(fes)
old = GridFunction(fes)

Apot, Tpot = gfu.components
oldApot, oldTpot = old.components
#dirich = GridFunction(fes)

#gfu.Set(bb, definedon=mesh.Boundaries('rub'))
#old.Set(bb, definedon=mesh.Boundaries('rub'))

omega=314
#rel = 1200
rho=CF( (0.1 , 0, 0,   0, 0.0005, 0,  0, 0, 0.0005), dims=(3,3) )
sig = mesh.MaterialCF({ "core" : 1 }, default=None)
#sigma=CoefficientFunction( (10, 0,  0, 10), dims=(2,2) )

rot=CF( (0 , 0, 0,   0, 0, 1,  0, -1, 0), dims=(3,3) )

mu0=4*pi*1e-7
Href=[0,42,53,62,70,79,88,100,113,132,157,193,255,376,677,1624,1e9]
Bref=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1e9*mu0]
nu=[Href[i]/Bref[i] if Bref[i]!=0 else 0 for i in range(len(Bref))]
#HBcurve = BSpline(4, [0,0,0]+list(Href), list(Bref)) #ovo bi trebala biti instanca klase BSpline
#diffHBcurve= HBcurve.Differentiate() #munonlin.Differentiate() metoda daje objekt klase BSpline

B=curl(oldApot)
Babs= B.Norm() 
errorlist=[]
p=1.0


for i in range(1,7):
    print(f"####iteration i={i}")
    
    print('Babs=',Babs(mesh(0,0)))
    
    #Babs += 1e-7
    #relyz= 1200 
    #dHdByz= 1200
    relyz= Babs*500 + 1e-5
    dHdByz= 2*Babs*500 + 1e-5
    
    rel=CF( ( 32000, 0, 0,   0, relyz, 0,  0, 0, relyz), dims=(3,3) )
    dHdB=CF( ( 32000, 0, 0,   0, dHdByz, 0,  0, 0, dHdByz), dims=(3,3) )

    #rel = HBcurve(Babs)  #Babs+1e-6
    #dHdB = diffHBcurve(Babs)

    #print('rel=', rel(mesh(0.6,0.6)))
    #print('dHdB=', dHdB(mesh(0.6,0.6)))

    term1=(1/mu0)*curl(mvp)*curl(alpha)*dx('air|coil') + rel*curl(mvp)*curl(alpha)*dx('core') \
    - (rot*grad(csp))*alpha*dx('core') #+ 0.1*mvp*alpha*dx('air|coil')
    term2= -1j/omega*rho*(rot*grad(csp))*(rot*grad(tau))*dx('core') + (rot*grad(tau))*mvp*dx('core')
    term3=1*1e0*mvp*alpha*dx

    jac= (dHdB - rel)*curl(mvp)*curl(alpha)*dx('core')
  
    a = BilinearForm(term1+term2+term3+jac)
    a.Assemble()

    jacmat= BilinearForm(jac)
    jacmat.Assemble()

    #rhs1 = dummy*v*dx
    #rhs1 = (dHdB - rel)*(Conj(B)*curl(old))*(B*curl(v))/(Babs*Babs)*dx('inner')
    
    #rhs1 = (dHdB - rel)*curl(oldApot)*curl(alpha) *dx('core')
    
    f = LinearForm(fes)
    Ts_air= IfPos(z*z-0.01, (0,0,0), IfPos((x*x+y*y-0.01),CF((0,0,0.015)),(0,0,0)))
    Ts_coil= CF((0,0,0.5*(x*x+y*y)-0.005))
    f += Ts_coil *1e7* curl(alpha) * dx("coil") + Ts_air *1e7* curl(alpha) * dx("air|core") 
    #curr = Ts_coil *1e7*curl(alpha)*dx("coil") + Ts_air *1e7*curl(alpha)*dx("air|core")
    #f = LinearForm(curr) # + rhs1)
    f.Assemble()


    ##### SOLVER
    #solvers.BVP(bf=a, lf=f, gf=gfu, pre=None, maxsteps=2000, print=True)
    r = f.vec + jacmat.mat * oldApot.vec #old.vec <=ILI
    gfu.vec.data = a.mat.Inverse(freedofs=fes.FreeDofs())*r
    #sol.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec
    
    #errfunc = (gfu - old)/gfu
    errfunc = (Apot - oldApot)/oldApot
    defon = mesh.Materials('core')
    error=Integrate(errfunc.Norm(), mesh, definedon=defon)
    print('error =', error)
    errorlist.append(error)

    #old.vec.data= p*gfu.vec + (1-p)*old.vec
    old.vec.data= gfu.vec
    
    B=curl(oldApot)
    Babs= B.Norm()

print('errorlist', errorlist)
#####POSTPROCESING

Jpost= rot*grad(Tpot)
Bpost = curl(Apot)

Pow=0.5*rho*Jpost*Conj(Jpost) 
Peddy=Integrate(Pow, mesh, order=5, definedon=defon)
print('Peddy=',Peddy) 

Draw (Bpost, mesh, "B")
Draw (Jpost, mesh, "J")


print('Babs=',Babs(mesh(0.045,0,0)))
print('Jabs=',Jpost.Norm()(mesh(0.045,0)))

