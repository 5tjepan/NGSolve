from ngsolve import *
from netgen.csg import *
import netgen.gui
import matplotlib.pyplot as plt
import sys
sys.argv = ["fun"]


def MakeGeometry():
    tocx=0.0 #-0.05
    tocy=0.0 #-0.05
    dbox=0.6
    geometry = CSGeometry()
    box = OrthoBrick(Pnt(0,0,-dbox),Pnt(dbox,dbox,dbox)).bc("outer")

    front= Plane(Pnt(tocx,tocy,-0.1), Vec(-1,0,0) ).bc("front").maxh(0.002)
    right= Plane(Pnt(tocx,tocy,-0.1), Vec(0,-1,0) ).bc("right")
    bot  = Plane (Pnt(tocx,tocy,-0.1), Vec(0,0,-1) ).bc("bot")
    back = Plane (Pnt(0.005,0.02,0.1), Vec(1, 0,0) ).bc("back").maxh(0.002)
    left = Plane (Pnt(0.005,0.02,0.1), Vec(0,1,0) ).bc("left")
    top  = Plane (Pnt(0.005,0.02,0.1), Vec(0,0, 1) ).bc("top")
    core = left * right * front * back * bot * top
    #core = OrthoBrick(Pnt(-0.05,-0.05,-0.2),Pnt(0.05,0.05,0.2))
    core.maxh(0.003)
    #front.maxh(0.01)
    core.mat("core")
    
    dno=-0.09
    vrh=0.01
    centy=0.035
    rin=0.012
    rout=0.02
    coil = ((Cylinder(Pnt(0,centy,-0.2), Pnt(0,centy,0.2), rout) - \
            Cylinder(Pnt(0,centy,-0.2), Pnt(0,centy,0.2), rin)) * \
            OrthoBrick (Pnt(0,centy,dno),Pnt(rout,(centy+rout),vrh)))+ \
            OrthoBrick (Pnt(rin,0,dno),Pnt(rout,centy,vrh))
    coil.maxh(0.009)
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


fsU = HCurl(mesh, order=0, dirichlet="outer", complex=True, nograds = False)
fsV = H1(mesh, order=1, dirichlet='left|bot|top', definedon='core', complex=True) #bez right!
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
rho=CF( (0.1 , 0, 0,   0, 5e-7, 0,  0, 0, 5e-7), dims=(3,3) )
sig = mesh.MaterialCF({ "core" : 1 }, default=None)
#sigma=CoefficientFunction( (10, 0,  0, 10), dims=(2,2) )

rot=CF( (0 , 0, 0,   0, 0, 1,  0, -1, 0), dims=(3,3) )

mu0=4*pi*1e-7

B=curl(oldApot)
Babs= B.Norm() 
errorlist=[]
p=1.0


for i in range(1,4):
    print(f"####iteration i={i}")
    
    print('Babs=',Babs(mesh(0,0)))
    
    #Babs += 1e-7
    #relyz= 30
    #dHdByz= 30
    
    a=-144.54207430769281
    b= 304.6915669443657
    c=-209.76353623263986
    d=48.34162716175372
    #e = np.e
    e=2.718281828459
        
    exponent = a + b*Babs + c*Babs**2 + d*Babs**3
    numerator = e**exponent * ((b + 2*c*Babs + 3*d*Babs**2) * Babs - 1)
    
    fun = e**exponent/(Babs + 1e-5)
    dfun = numerator/(Babs**2 + 1e-8)

    fit1= 9.365857*Babs + 13.27562*Babs**0.5
    dfit1= 9.365857 + 13.27562*0.5*(Babs+1e-8)**(-0.5)

    relyz= IfPos(1.3-Babs, fit1, fun) + 1e-5
    dHdByz= IfPos(1.3-Babs, dfit1, dfun ) + 1e-5

    #relyz= 25.31501 + 0.000005078657*Babs**(eksp-1) + 1e-5   #Babs*500 + 1e-5
    #dHdByz=25.31501 + 0.000005078657*eksp*Babs**(eksp-1) + 1e-5      #2*Babs*500 + 1e-5
    
    rel=CF( ( 32000, 0, 0,   0, relyz, 0,  0, 0, relyz), dims=(3,3) )
    dHdB=CF( ( 32000, 0, 0,   0, dHdByz, 0,  0, 0, dHdByz), dims=(3,3) )

    term1=(1/mu0)*curl(mvp)*curl(alpha)*dx('air|coil') + rel*curl(mvp)*curl(alpha)*dx('core') \
    - (rot*grad(csp))*alpha*dx('core') #+ 0.1*mvp*alpha*dx('air|coil')
    term2= -1j/omega*rho*(rot*grad(csp))*(rot*grad(tau))*dx('core') + (rot*grad(tau))*mvp*dx('core')
    term3=1*1e0*mvp*alpha*dx

    jac= (dHdB - rel)*curl(mvp)*curl(alpha)*dx('core')
  
    a = BilinearForm(term1+term2+term3+jac)
    a.Assemble()

    jacmat= BilinearForm(jac)
    jacmat.Assemble()

    f = LinearForm(fes)
    I=5 #Amp
    zavoj=447
    dno=-0.09
    vrh=0.01
    centy=0.035
    rin=0.012
    rout=0.02
    R=(x**2 + (y-centy)**2)**0.5
    Js=1.414*I*zavoj/((rout-rin)*(vrh-dno))
    izvan=CF((0,0,Js*(rout-rin)))
    nula=CF((0,0,0))

    Ts_coil=IfPos(y-centy, CF((0,0,Js*(R-rin))),CF((0,0,Js*(x-rin))))
    Ts_air=IfPos(rin-R, nula, IfPos(rin-x, IfPos(y-centy,izvan,nula), izvan)) * IfPos((z-vrh)*(z-dno),0,1)
    f += Ts_coil *curl(alpha) * dx("coil") + Ts_air *curl(alpha) * dx("air|core") 
    
    f.Assemble()


    ##### SOLVER
    #solvers.BVP(bf=a, lf=f, gf=gfu, pre=None, maxsteps=2000, print=True)
    r = f.vec + jacmat.mat * oldApot.vec #old.vec <=ILI
    gfu.vec.data = a.mat.Inverse(freedofs=fes.FreeDofs())*r
    
    #errfunc = (gfu - old)/gfu
    errfunc = (Apot - oldApot).Norm()/(oldApot.Norm()+1e-15)
    #eps=CF((1e-15,1e-15,1e-15))
    #errfunc = (IfPos(Apot.Norm(),Apot,eps) - IfPos(oldApot.Norm(),oldApot,eps))/IfPos(oldApot.Norm(),oldApot,eps)

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

print('...')
Pow=0.5*rho*Jpost*Conj(Jpost) 
Peddy=Integrate(Pow, mesh, order=5, definedon=defon)
print('Peddy=',round(abs(Peddy),4),'W') 

volumen=Integrate(1,mesh, definedon=defon)
BdV=Integrate(Bpost.Norm(), mesh, order=5, definedon=defon)
print('Bavg=',round(BdV/volumen, 3),'T') 

Draw (Bpost, mesh, "B")
Draw (Jpost, mesh, "J")


