import netgen.gui
import ngsolve
from netgen.read_gmsh import ReadGmsh
from netgen.meshing import *
from ngsolve import *

#import matplotlib.pyplot as plt
#import sys
#sys.argv = ["fun"]

ngsglobals.msg_level = 5  #print informacije u log

mesh_g = ReadGmsh("bulkEnpay39k.msh")

mesh = ngsolve.Mesh(mesh_g)
Draw(mesh)

print(mesh.GetBoundaries())
#----------------------


#**************************************
######### HBcurve --- BSpline
mu0=4*pi*1e-7

B_ref=[0.0, 0.04, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.83, 1000]
H_ref=[0.0, 3.5, 5.0, 7.6, 10.0, 12.07, 14.08, 16.0, 17.75, 19.42, 21.05, 22.64, 24.16, 25.76, 27.98, 29.66, 32.05, 35.33, 40.06, 47.17, 58.35, 77.19, 113.0, 196.0, 395.0, 800, 1000/mu0]
#H_ref=[0.0, 3.7, 5.5, 8.0, 10.1, 12.1, 14.08, 16.0, 17.75, 19.42, 21.05, 22.64, 24.16, 25.76, 27.98, 29.66, 32.05, 35.33, 40.06, 47.17, 58.35, 77.19, 113.0, 196.0, 395.0, 800, 1000/mu0]
#H_rms=[0.0, 2.48, 3.5, 5.4, 7.07, 8.61, 10.07, 11.44, 12.75, 14.03, 15.31, 16.62, 18.00, 19.50, 21.22, 22.23, 23.39, 24.75, 26.48, 28.80, 32.21, 37.69, 48.03, 71.89, 131.57, 240, 1000/mu0]
#sljedeci H_ref je jednak Hrms*sqrt(2) iz tablice lima M140-35S:
#H_ref=[0.0, 3.5, 4.9, 7.6, 10.0, 12.2, 14.2, 16.2, 18.0, 19.8, 21.7, 23.5, 25.5, 27.6, 30.0, 31.4, 33.1, 35.0, 37.4, 40.7, 45.6, 53.3, 67.9, 101.7, 186.1, 339.4, 1414/mu0]

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


#graddom = [True if mat == "core" else False for mat in mesh.GetMaterials()]
fsU = HCurl(mesh, order=0, dirichlet="gamaB|front|gamaE|right", complex=True, nograds = False)#, gradientdomains = graddom)

#updateHCurlRegionOrder(fsU, 2, "core")
#------------

#fsU = HCurl(mesh, order=0, dirichlet="outer", complex=True, nograds = False)
fsV = H1(mesh, order=1, dirichlet='bottom|top|left', definedon='core', complex=True) #bez right!
fes=fsU*fsV
mvp, csp = fes.TrialFunction()
alpha, tau = fes.TestFunction()

print('fes.ndof=',fes.ndof)
print('...free =', sum(fes.FreeDofs()))

gfu = GridFunction(fes)
old = GridFunction(fes)

Apot, Tpot = gfu.components
oldApot, oldTpot = old.components

nu_fe=30.0  #420
omega=314
d=0.00035 #m
Kf= 27*d/0.01 #0.945
kappa = 2e6
rho=CF( (1 , 0, 0,   0, 1/(Kf*kappa), 0,  0, 0, 1/(Kf*kappa)), dims=(3,3) )
sig = mesh.MaterialCF({ "core" : 1 }, default=None)
#psi = mesh.MaterialCF({ "air_1" : 1 }, default=None)

rot=CF( (0 , 0, 0,   0, 0, 1,  0, -1, 0), dims=(3,3) )

B=curl(oldApot)
Babs= B.Norm() 
errorlist=[]
p=1.0


for i in range(1,2):
    print(f"####iteration i={i}")
    
    print('Babs=',Babs(mesh(0,0)))
    
    #relyz= 1/(5000*mu0) #1/(27000*mu0)  + 1j*omega*kappa*d**2*1/12 #pazi, kappa_y doprinosi rel_z i obratno  
    #dHdByz= 1/(500*mu0)

    #relyz = (HBcurve(Babs+1e-6))/(Babs+1e-6) + 1j*omega*kappa*d**2*1/12
    #dHdByz = diffHB(Babs+1e-6) + 1j*omega*kappa*d**2*1/12 
    
    relyz = nu_fe + 1j*omega*kappa*d**2*1/12 
    #dHdB=CF( ( (1-Kf)/mu0, 0, 0,   0, dHdByz/Kf, 0,  0, 0, dHdByz/Kf), dims=(3,3) ) #?? should I divide by Kf

    rel=CF( ( 1*(1-Kf)/mu0, 0, 0,   0, relyz/Kf, 0,  0, 0, relyz/Kf), dims=(3,3) )
    


    term1=(1/mu0)*curl(mvp)*curl(alpha)*dx('air|coil') + rel*curl(mvp)*curl(alpha)*dx('core') \
    - (rot*grad(csp))*alpha*dx('core') #+ 0.1*mvp*alpha*dx('air|coil')
    term2= -1j/omega*rho*(rot*grad(csp))*(rot*grad(tau))*dx('core') + (rot*grad(tau))*mvp*dx('core') #tau*CF((1,0,0))*curl(mvp)*dx('core') #
    rho_eff = 1*d**2 /12 /nu_fe/Kf
    term3=1j*1e2*mvp*alpha*dx + rho_eff *(rot*grad(csp))*(rot*grad(tau))*dx('core')

    #jac= (dHdB - rel)*curl(mvp)*curl(alpha)*dx('core')
  
    a = BilinearForm(term1+term2+term3) #+jac)
    a.Assemble()

    #jacmat= BilinearForm(jac)
    #jacmat.Assemble()

    #:::::::::::::::::: source current
    #f = LinearForm(fes)
    I=3.5 ##0.75*(1j+1)/1.414213562373 #Amp
    zavoj=447
    dno=-0.08
    vrh=0.02
    centy=0.035
    rin=0.012
    rout=0.02
    R=(x**2 + (y-centy)**2)**0.5
    Js=1.414*I*zavoj/((rout-rin)*(vrh-dno))
    izvan=CF((0,0,Js*(rout-rin)))
    nula=CF((0,0,0))

    Ts_coil=IfPos(y-centy, CF((0,0,Js*(R-rin))),CF((0,0,Js*(x-rin)))) * IfPos((z-vrh)*(z-dno),0,1)
    Ts_air=IfPos(rin-R, nula, IfPos(rin-x, IfPos(y-centy,izvan,nula), izvan)) * IfPos((z-vrh)*(z-dno),0,1)
    Ts = mesh.MaterialCF({ "coil" : Ts_coil, "core|air" : Ts_air })
    Ts_total = Ts *curl(alpha) * dx("coil") + Ts *curl(alpha) * dx("air|core") 
    #Ts_total = Ts_coil *curl(alpha) * dx("coil") + Ts_air *curl(alpha) * dx("air|core") 
    
    f=LinearForm(Ts_total)
    f.Assemble()
    #::::::::::::::::

    #pc = Preconditioner(a, type="multigrid")
    ##### SOLVER
    r = f.vec #+ jacmat.mat * oldApot.vec #old.vec <=ILI
    #gfu.vec.data = a.mat.Inverse(freedofs=fes.FreeDofs())*r
    
    r_bvp = LinearForm(fes).Assemble()
    r_bvp.vec.data += r
    solvers.BVP(bf=a, lf=r_bvp, gf=gfu, pre=None, maxsteps=200, print=True, inverse="pardiso", needsassembling=False)
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

#p_narrow= Kf * kappa/24 *(omega*d)**2 *(Bpost.Norm())**2 #
p_narrow= kappa/24 *(omega*d)**2 *(Btang)**2 /Kf  #!!! da, treba ici kroz Kf
P_xyz = Integrate(p_narrow, mesh,order=5, definedon=defon) 
print('P_xyz=',1e3*round(P_xyz,9), 'mW')

q_narrow=1* omega/24 * d**2/nu_fe * Jpost*Conj(Jpost) /Kf  #!!! da, treba ici kroz Kf
Q_xyz = Integrate(q_narrow, mesh,order=5, definedon=defon) 
print('Q_xyz=',1e3*round(abs(Q_xyz),8), 'mW')

p_wide=0.5*rho*Jpost*Conj(Jpost) 
P_yz=Integrate(p_wide, mesh, order=5, definedon=defon)
print('P_yz=', 1e3*round(abs(P_yz),9),'mW')

print('P_EC=',1e3*round(abs(P_xyz+P_yz),4),'mW')


B_eps=[0.0, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 2.5]
Pm3=[0, 612, 842, 1224, 1607, 2066, 2601, 3213, 3825, 4514, 5355, 5738, 6197, 6732, 7268, 7803, 8492, 9180, 10175, 11475,13158,35000]

EPScurve = BSpline(2, [0]+list(B_eps), list(Pm3)) #ovo bi trebala biti instanca klase BSpline
Peps_density = (EPScurve(Btang+1e-6))
Peps_total=Integrate(Peps_density, mesh, order=5, definedon=defon)
#print('P_eps=',1e3*round(Peps_total,4), 'mW')

#===============
volumen=Integrate(1,mesh, definedon=defon)
#BdV=Integrate(Bpost.Norm(), mesh, order=5, definedon=defon)
BzdV=Integrate(Bpost[2], mesh, order=5, definedon=defon)
print('Bavg=',round(abs(BzdV)/volumen, 4),'T') 
#print('Bavg=',(BzdV/volumen),'T') 

BxdV=Integrate(Bpost[0].Norm(), mesh, order=5, definedon=defon)
print('Bx_avg=',round(BxdV/volumen, 8),'T') 

#Btang= (Bpost[1].Norm()**2+Bpost[2].Norm()**2)
#BdV=Integrate(Btang, mesh, order=5, definedon=defon)
#print('Btang**2=',round(BdV/volumen, 6),'T') 

nu_CF=CF( ( 1*(1-Kf)/mu0, 0, 0,   0, nu_fe/Kf, 0,  0, 0, nu_fe/Kf), dims=(3,3) )
BBdV=Integrate( omega/2*Bpost*(nu_CF*Conj(Bpost)), mesh, order=5, definedon=mesh.Materials('core'))
print('Q =',1e3*round(BBdV.real, 5),'mVAr') 
#relBBdV=Integrate( omega/2*Bpost*(rel*Conj(Bpost)), mesh, order=5, definedon=mesh.Materials('core'))
relBBdV=Integrate( 1j*omega/2*Bpost*Conj(rel*Bpost), mesh, order=5, definedon=mesh.Materials('core'))
print('S =',1e3*(P_yz + relBBdV),'mVA') 

#print('P_tot=',1e3*round(abs(Peps_total+P_yz),4),'mW') 
Draw (Bpost, mesh, "B",draw_surf=False)
Draw (Jpost, mesh, "J",draw_surf=False)
Draw (Ts, mesh, "Ts")


#mesh.ngmesh.SetMaterial(2,'two')
#or even better:
#domain=CF([i+1 for i in range(len(mesh.GetMaterials()))])
#...
#Integrate(1,mesh,region_wise=True)
