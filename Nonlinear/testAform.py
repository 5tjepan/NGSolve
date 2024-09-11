from ngsolve import *
from netgen.csg import *
import netgen.gui
import time
import matplotlib.pyplot as plt
import sys
sys.argv = ["fun"]


def MakeGeometry():
    tocx=0.0 #-0.05
    tocy=0.0 #-0.05
    dbox=0.6
    geometry = CSGeometry()
    box = OrthoBrick(Pnt(0,0,-dbox),Pnt(dbox,dbox,dbox)).bc("outer")

    front= Plane(Pnt(tocx,tocy,-0.1), Vec(-1,0,0) ).bc("front")#
    right= Plane(Pnt(tocx,tocy,-0.1), Vec(0,-1,0) ).bc("right")
    bot  = Plane (Pnt(tocx,tocy,-0.1), Vec(0,0,-1) ).bc("bot")
    back = Plane (Pnt(0.005,0.02,0.1), Vec(1, 0,0) ).bc("back").maxh(0.001)
    left = Plane (Pnt(0.005,0.02,0.1), Vec(0,1,0) ).bc("left")
    top  = Plane (Pnt(0.005,0.02,0.1), Vec(0,0, 1) ).bc("top")
    core = left * right * front * back * bot * top
    #core = OrthoBrick(Pnt(-0.05,-0.05,-0.2),Pnt(0.05,0.05,0.2))
    core.maxh(0.0015)
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
    coil.maxh(0.005)
    coil.mat("coil")
    
    geometry.Add ((box-core-coil).mat("air"))
    geometry.Add (core)
    geometry.Add (coil)

    return geometry


ngmesh = MakeGeometry().GenerateMesh(maxh=0.15)
#ngmesh.Save("coil.vol")
mesh = Mesh(ngmesh)
mesh.Curve(5)  

Draw(mesh)
print(mesh.GetBoundaries())


########## HBcurve --- BSpline
mu0=4*pi*1e-7

B_ref=[0.0, 0.04, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.74, 1.77, 1000]
H_ref=[0.0, 3.5, 5.0, 7.6, 10.0, 12.07, 14.08, 16.0, 17.75, 19.42, 21.05, 22.64, 24.16, 25.76, 27.98, 29.66, 32.05, 35.33, 40.06, 47.17, 58.35, 77.19, 107.0, 160.0, 270.0, 1000/mu0]

Bvis = Parameter(0) #klasa Parameter je dijete klase CF. Uvrstavanjem Hvis u pozivu BSpline() je poziv oblika BSpline(CF)...jer je __call__ overloadan

HBcurve = BSpline(2, [0]+list(B_ref), list(H_ref)) #ovo bi trebala biti instanca klase BSpline
diffHB= HBcurve.Differentiate() #HBcurve.Differentiate() metoda daje objekt klase BSpline

# v i s u a l i s a t i o n :::
HBcurve_vis = HBcurve((1e-6+sqrt(Bvis*Bvis+1e-6)))#/(1e-6+sqrt(Hvis*Hvis+1e-6)) #prema __call__ poziv BSplinee(CF) daje CF.
diffHB_vis= diffHB((1e-6+sqrt(Bvis*Bvis+1e-6))) #buduci da je diffHB objekt klase BSpline, prema __call__ poziv BSplinee(CF) daje CF

Hvec=[]
Hder=[]
Bvec=[i*1.78/100 for i in range(1,102)]
for k in Bvec:
    Bvis.Set(k)
    Hvec.append(HBcurve_vis(mesh())) #buduci da je HBcurve_vis CF funkcija, prema __call__ klase CF za dobivanje vrijednostii CF moze se koristiti input "mesh()"
    Hder.append(diffHB_vis(mesh())) 
    
#plt.xlabel('B')
#plt.scatter(B_ref[:-1],H_ref[:-1])
#plt.plot(Bvec, Hvec)
#plt.plot(Bvec, Hder)
#plt.show()
##################


def updateHCurlRegionOrder(fes, p, mat):
     for el in fes.Elements():
         if el.mat == mat:
             fes.SetOrder(NodeId(ELEMENT, el.nr), p)
             for f in el.faces:
                 fes.SetOrder(NodeId(FACE, f.nr), p)

             for ed in el.edges:
                 fes.SetOrder(NodeId(EDGE, ed.nr), p)

             for v in el.vertices:
                 fes.SetOrder(NodeId(EDGE, v.nr), p)

     fes.Update()


graddom = [True if mat == "core" else False for mat in mesh.GetMaterials()]
fes = HCurl(mesh, order=0, dirichlet="outer", complex=True)#, gradientdomains = graddom)

updateHCurlRegionOrder(fes, 1, "core")

print('fes.ndof=',fes.ndof)
print('...free =', sum(fes.FreeDofs()))

#fes = HCurl(mesh, order=1, dirichlet="outer", complex=True, nograds = False)
mvp = fes.TrialFunction()
alpha = fes.TestFunction()

Apot = GridFunction(fes)
oldApot = GridFunction(fes)


omega=314
#rel = 1200
d=0.00035 #m
Kf= 27*d/0.01 #0.945
kappa=2e6
rho=CF( (0.1 , 0, 0,   0, 5e-7, 0,  0, 0, 5e-7), dims=(3,3) )
sigma=CoefficientFunction( (0.01 , 0, 0,   0, kappa*Kf, 0,  0, 0, kappa*Kf), dims=(3,3) )
sig = mesh.MaterialCF({ "core" : 1 }, default=None)


B=curl(oldApot)
Babs= B.Norm() 
errorlist=[]
p=1.0


for i in range(1,2):
    print(f"####iteration i={i}")
    
    print('Babs=',Babs(mesh(0,0)))
    #print('B =',B(mesh(0,0)))
    
    #relyz= Babs*500 + 1e-5
    #dHdByz= 2*Babs*500 + 1e-5
    
    relyz= 160 + 1j*omega*kappa*d**2*1/12 #pazi, kappa_y doprinosi rel_z i obratno  
    dHdByz= 160 + 1j*omega*kappa*d**2*1/12

    #relyz = (HBcurve(Babs+1e-6))/(Babs+1e-6) + 1j*omega*kappa*d**2*1/12 #Babs+1e-6
    #dHdByz = diffHB(Babs+1e-6) + 1j*omega*kappa*d**2*1/12
    
    rel=CF( ( (1-Kf)/mu0, 0, 0,   0, relyz/Kf, 0,  0, 0, relyz/Kf), dims=(3,3) )
    dHdB=CF( ( (1-Kf)/mu0, 0, 0,   0, dHdByz/Kf, 0,  0, 0, dHdByz/Kf), dims=(3,3) )

    #rel = HBcurve(Babs)  #Babs+1e-6
    #dHdB = diffHBcurve(Babs)

    #print('rel=', rel(mesh(0.6,0.6)))
    #print('dHdB=', dHdB(mesh(0.6,0.6)))

    term1 = (1/mu0)*curl(mvp)*curl(alpha)*dx('air|coil') + rel*curl(mvp)*curl(alpha)*dx('core') + \
    1j*omega*sigma*mvp*alpha*dx('core')
    term2 = 0.01*1e0*mvp*alpha*dx('air|coil')
    
    jac= (dHdB - rel)*curl(mvp)*curl(alpha)*dx('core')
  
    a = BilinearForm(term1 + term2 + jac)
    a.Assemble()

    jacmat= BilinearForm(jac)
    jacmat.Assemble()

#######################   
    f = LinearForm(fes)
    I=4 #Amp
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
    start=time.time()
    r = f.vec + jacmat.mat * oldApot.vec #- a.mat * gfu.vec
    #Apot.vec.data = a.mat.Inverse(freedofs=fes.FreeDofs())*r
    
    r_bvp = LinearForm(fes).Assemble()
    r_bvp.vec.data += r
    solvers.BVP(bf=a, lf=r_bvp, gf=Apot, pre=None, maxsteps=2000, print=True, needsassembling=False)
    stop=time.time()
    #print('duration=',stop-start)
    #------------------

    #errfunc = (gfu - old)/gfu
    errfunc = (Apot - oldApot).Norm()/oldApot.Norm()
    defon = mesh.Materials('core')
    error=Integrate(errfunc, mesh, definedon=defon)
    print('error =', error)
    errorlist.append(error)

    #old.vec.data= p*gfu.vec + (1-p)*old.vec
    oldApot.vec.data= Apot.vec
    
    B=curl(oldApot)
    Babs= B.Norm()

print('errorlist', errorlist)
#####POSTPROCESING


Jpost = - 1j * omega*sig * sigma * Apot
Bpost = curl(Apot)
Btang= (Bpost[1].Norm()**2+Bpost[2].Norm()**2)**0.5

print('Aform:::')
""" Pow=0.5*rho*Jpost*Conj(Jpost) 
Peddy=Integrate(Pow, mesh, order=5, definedon=defon)
print('Peddy=',round(abs(Peddy),4),'W') """

print('...')

#p_e=3159.0*(Bpost.Norm())**2.074
p_e=3159.0*(Btang)**2.074
P_eps=Integrate(p_e,mesh, order=5, definedon=defon)
print('P_eps=',1e3*round(P_eps,4), 'mW')

#p_narrow= Kf * kappa/24 *(omega*d)**2 *(Bpost.Norm())**2 #vjerojatno treba ići kroz Kf
p_narrow= Kf * kappa/24 *(omega*d)**2 *(Btang)**2 #vjerojatno treba ići kroz Kf
P_xyz = Integrate(p_narrow, mesh,order=5, definedon=defon) 
print('P_xyz=',1e3*round(P_xyz,4), 'mW')

p_wide=0.5*rho*Jpost*Conj(Jpost) 
P_yz=Integrate(p_wide, mesh, order=5, definedon=defon)
print('P_yz=',1e3*round(abs(P_yz),4),'mW')

print('P_tot=',1e3*round(abs(P_eps+P_yz),4),'mW') 

volumen=Integrate(1,mesh, definedon=defon)
BdV=Integrate(Bpost.Norm(), mesh, order=5, definedon=defon)
print('Bavg=',round(BdV/volumen, 3),'T') 

Draw (Bpost, mesh, "B")
Draw (Jpost, mesh, "J")

print('...')
print('duration=', round(stop-start,1))


#print('Babs=',Babs(mesh(0.045,0,0)))
#print('Jabs=',Jpost.Norm()(mesh(0.045,0)))
