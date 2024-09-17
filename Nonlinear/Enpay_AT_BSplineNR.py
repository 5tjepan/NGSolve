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

    #box = OrthoBrick(Pnt(0,0,-dbox),Pnt(dbox,dbox,dbox)).bc("outer")
    front_outer= Plane(Pnt(tocx,tocy,-0.1), Vec(-1,0,0) ).bc("outfront")
    right_outer= Plane(Pnt(tocx,tocy,-0.1), Vec(0,-1,0) ).bc("outright")
    bot_outer  = Plane (Pnt(0,0,-dbox), Vec(0,0, -1) ).bc("outer")
    back_outer = Plane (Pnt(dbox,0,0), Vec(1, 0,0) ).bc("outer")
    left_outer = Plane (Pnt(0,dbox,0), Vec(0,1,0) ).bc("outer")
    top_outer  = Plane (Pnt(0,0,dbox), Vec(0,0, 1) ).bc("outer")
    box = (bot_outer * front_outer * right_outer * back_outer * left_outer * top_outer)

    front= Plane(Pnt(tocx,tocy,-0.1), Vec(-1,0,0) ).bc("front")#
    right= Plane(Pnt(tocx,tocy,-0.1), Vec(0,-1,0) ).bc("right")
    bot  = Plane (Pnt(tocx,tocy,-0.1), Vec(0,0,-1) ).bc("bot")
    back = Plane (Pnt(0.005,0.02,0.1), Vec(1, 0,0) ).bc("back").maxh(0.0021)
    left = Plane (Pnt(0.005,0.02,0.1), Vec(0,1,0) ).bc("left")
    top  = Plane (Pnt(0.005,0.02,0.1), Vec(0,0, 1) ).bc("top")
    core = left * right * front * back * bot * top
    #core = OrthoBrick(Pnt(-0.05,-0.05,-0.2),Pnt(0.05,0.05,0.1))
    core.maxh(0.0021) 
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
    #geometry.AddSurface(right,core) ##surface!

    return geometry


ngmesh = MakeGeometry().GenerateMesh(maxh=0.15)
mesh = Mesh(ngmesh)
mesh.Curve(5)  

#ngmesh.Save("coil.vol")

Draw(mesh)
print(mesh.GetBoundaries())
#----------------------


########## HBcurve --- BSpline
mu0=4*pi*1e-7

B_ref=[0.0, 0.04, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.83, 1000]
H_ref=[0.0, 3.5, 5.0, 7.6, 10.0, 12.07, 14.08, 16.0, 17.75, 19.42, 21.05, 22.64, 24.16, 25.76, 27.98, 29.66, 32.05, 35.33, 40.06, 47.17, 58.35, 77.19, 113.0, 196.0, 395.0, 800, 1000/mu0]
#H_ref=[0.0, 3.7, 5.5, 8.0, 10.1, 12.1, 14.08, 16.0, 17.75, 19.42, 21.05, 22.64, 24.16, 25.76, 27.98, 29.66, 32.05, 35.33, 40.06, 47.17, 58.35, 77.19, 113.0, 196.0, 395.0, 800, 1000/mu0]
#H_ref=[0.0, 2.48, 3.5, 5.4, 7.07, 8.61, 10.07, 11.44, 12.75, 14.03, 15.31, 16.62, 18.00, 19.50, 21.22, 22.23, 23.39, 24.75, 26.48, 28.80, 32.21, 37.69, 48.03, 71.89, 131.57, 240, 1000/mu0]

HBcurve = BSpline(2, [0]+list(B_ref), list(H_ref)) #ovo bi trebala biti instanca klase BSpline
diffHB= HBcurve.Differentiate() #HBcurve.Differentiate() metoda daje objekt klase BSpline

# v i s u a l i s a t i o n :::
Bvis = Parameter(0) #klasa Parameter je dijete klase CF. Uvrstavanjem Hvis u pozivu BSpline() je poziv oblika BSpline(CF)...jer je __call__ overloadan

HBcurve_vis = HBcurve((1e-6+sqrt(Bvis*Bvis+1e-6)))#/(1e-6+sqrt(Hvis*Hvis+1e-6)) #prema __call__ poziv BSplinee(CF) daje CF.
diffHB_vis= diffHB((1e-6+sqrt(Bvis*Bvis+1e-6))) #buduci da je diffHB objekt klase BSpline, prema __call__ poziv BSplinee(CF) daje CF

Hvec=[]
Hder=[]
Bvec=[i*1.0/100 for i in range(1,102)]
for k in Bvec:
    Bvis.Set(k)
    Hvec.append(HBcurve_vis(mesh())) #buduci da je HBcurve_vis CF funkcija, prema __call__ klase CF za dobivanje vrijednostii CF moze se koristiti input "mesh()"
    Hder.append(diffHB_vis(mesh())) 
    
#plt.xlabel('B')
#plt.scatter(B_ref[:-12],H_ref[:-12])
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
fsU = HCurl(mesh, order=0, dirichlet="outer|outfront|outright", complex=True, nograds = False)#, gradientdomains = graddom)

#updateHCurlRegionOrder(fsU, 0, "core")
#------------

#fsU = HCurl(mesh, order=0, dirichlet="outer", complex=True, nograds = False)
fsV = H1(mesh, order=1, dirichlet='left|bot|top', definedon='core', complex=True) #bez right!
fes=fsU*fsV
mvp, csp = fes.TrialFunction()
alpha, tau = fes.TestFunction()

print('fes.ndof=',fes.ndof)
print('...free =', sum(fes.FreeDofs()))

gfu = GridFunction(fes)
old = GridFunction(fes)

Apot, Tpot = gfu.components
oldApot, oldTpot = old.components

omega=314
d=0.00035 #m
Kf= 27*d/0.01 #0.945
kappa = 2e6
rho=CF( (0.1 , 0, 0,   0, 1/(Kf*kappa), 0,  0, 0, 1/(Kf*kappa)), dims=(3,3) )
sig = mesh.MaterialCF({ "core" : 1 }, default=None)

rot=CF( (0 , 0, 0,   0, 0, 1,  0, -1, 0), dims=(3,3) )


B=curl(oldApot)
Babs= B.Norm() 
errorlist=[]
p=1.0


for i in range(1,8):
    print(f"####iteration i={i}")
    
    print('Babs=',Babs(mesh(0,0)))
    
    #relyz= 30 + 1j*omega*kappa*d**2*1/12 #pazi, kappa_y doprinosi rel_z i obratno  
    #dHdByz= 30 + 1j*omega*kappa*d**2*1/12

    relyz = (HBcurve(Babs+1e-6))/(Babs+1e-6) + 1j*omega*kappa*d**2*1/12 #Babs+1e-6
    dHdByz = diffHB(Babs+1e-6) + 1j*omega*kappa*d**2*1/12
    
    rel=CF( ( (1-Kf)/mu0, 0, 0,   0, relyz/Kf, 0,  0, 0, relyz/Kf), dims=(3,3) )
    dHdB=CF( ( (1-Kf)/mu0, 0, 0,   0, dHdByz/Kf, 0,  0, 0, dHdByz/Kf), dims=(3,3) ) #?? should I divide with Kf

    term1=(1/mu0)*curl(mvp)*curl(alpha)*dx('air|coil') + rel*curl(mvp)*curl(alpha)*dx('core') \
    - (rot*grad(csp))*alpha*dx('core') #+ 0.1*mvp*alpha*dx('air|coil')
    term2= -1j/omega*rho*(rot*grad(csp))*(rot*grad(tau))*dx('core') + (rot*grad(tau))*mvp*dx('core')
    term3=1*1e0*mvp*alpha*dx

    jac= (dHdB - rel)*curl(mvp)*curl(alpha)*dx('core')
  
    a = BilinearForm(term1+term2+term3+jac)
    a.Assemble()

    jacmat= BilinearForm(jac)
    jacmat.Assemble()

    #:::::::::::::::::: source current
    #f = LinearForm(fes)
    I=0.75 #Amp
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
    Ts_total = Ts_coil *curl(alpha) * dx("coil") + Ts_air *curl(alpha) * dx("air|core") 
    #f += Ts_coil *curl(alpha) * dx("coil") + Ts_air *curl(alpha) * dx("air|core") 
    
    f=LinearForm(Ts_total)
    f.Assemble()
    #::::::::::::::::

    #pc = Preconditioner(a, type="multigrid")
    ##### SOLVER
    start=time.time()
    r = f.vec + jacmat.mat * oldApot.vec #old.vec <=ILI
    #gfu.vec.data = a.mat.Inverse(freedofs=fes.FreeDofs())*r
    
    r_bvp = LinearForm(fes).Assemble()
    r_bvp.vec.data += r
    solvers.BVP(bf=a, lf=r_bvp, gf=gfu, pre=None, maxsteps=2000, print=True, needsassembling=False)
    stop=time.time()
    #----------------------


    """ rhs_jac =  jacmat.mat * oldApot.vec #old.vec <=ILI
    r=LinearForm(Ts_total+rhs_jac)
    r.Assemble()
    solvers.BVP(bf=a, lf=r, gf=gfu, pre=None, maxsteps=2000, print=True) """
    
    #errfunc = (gfu - old).Norm()/(old.Norm()+1e-15)
    errfunc = (Apot - oldApot).Norm()/(oldApot.Norm()+1e-15)
    
    defon = mesh.Materials('core')
    error=Integrate(errfunc.Norm(), mesh) #, definedon=defon)
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

Btang= (Bpost[1].Norm()**2+Bpost[2].Norm()**2)**0.5

print('...')

""" #p_e=3159.0*(Bpost.Norm())**2.074
p_e = 3159.0*(Btang)**2.074
P_eps=Integrate(p_e,mesh, order=5, definedon=defon)
print('P_eps=',1e3*round(P_eps,4), 'mW') """

#p_narrow= Kf * kappa/24 *(omega*d)**2 *(Bpost.Norm())**2 #vjerojatno treba ići kroz Kf
p_narrow= Kf * kappa/24 *(omega*d)**2 *(Btang)**2 #vjerojatno treba ići kroz Kf
P_xyz = Integrate(p_narrow, mesh,order=5, definedon=defon) 
print('P_xyz=',1e3*round(P_xyz,4), 'mW')

p_wide=0.5*rho*Jpost*Conj(Jpost) 
P_yz=Integrate(p_wide, mesh, order=5, definedon=defon)
print('P_yz=', 1e3*round(abs(P_yz),4),'mW')

#print('P_tot-old=',1e3*round(abs(P_eps+P_yz),4),'mW')


B_eps=[0.0, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 2.5]
Pm3=[0, 612, 842, 1224, 1607, 2066, 2601, 3213, 3825, 4514, 5355, 5738, 6197, 6732, 7268, 7803, 8492, 9180, 10175, 11475,13158,35000]

EPScurve = BSpline(2, [0]+list(B_eps), list(Pm3)) #ovo bi trebala biti instanca klase BSpline
Peps_density = (EPScurve(Btang+1e-6))
Peps_total=Integrate(Peps_density, mesh, order=5, definedon=defon)
print('P_eps=',1e3*round(Peps_total,4), 'mW')


print('P_tot=',1e3*round(abs(Peps_total+P_yz),4),'mW') 

volumen=Integrate(1,mesh, definedon=defon)
BdV=Integrate(Bpost.Norm(), mesh, order=5, definedon=defon)
print('Bavg=',round(BdV/volumen, 3),'T') 

Draw (Bpost, mesh, "B")
Draw (Jpost, mesh, "J")


print('...')
print('duration=', round(stop-start,1))
