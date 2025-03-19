import netgen.gui
import ngsolve
from netgen.read_gmsh import ReadGmsh
from netgen.meshing import *
from ngsolve import *

import matplotlib.pyplot as plt
import sys
sys.argv = ["fun"]

ngsglobals.msg_level = 5

#Naziv 8BCC oznacava OsminuBoxCoreCoil
# Učitavanje Gmsh mreže
mesh_g = ReadGmsh("BoxCoreCoil_eighth1.msh")

#mesh = Mesh(mesh)
#mesh.ngmesh.Save("disk.vol")
#Draw(mesh)

mesh = ngsolve.Mesh(mesh_g)
Draw(mesh)

print(mesh.GetBoundaries())
#----------------------


########## HBcurve --- BSpline
mu0=4*pi*1e-7

B_ref=[0.0, 0.04, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.83, 1000]
H_ref=[0.0, 3.5, 4.9, 7.6, 10.0, 12.2, 14.2, 16.2, 18.0, 19.8, 21.7, 23.5, 25.5, 27.6, 30.0, 31.4, 33.1, 35.0, 37.4, 40.7, 45.6, 53.3, 67.9, 101.7, 186.1, 339.4, 1414/mu0]
B_ref2=[elem*1.25 for elem in B_ref]
H_ref2=[elem*0.9 for elem in H_ref]

HBcurve = BSpline(2, [0]+list(B_ref), list(H_ref)) #ovo bi trebala biti instanca klase BSpline
diffHB= HBcurve.Differentiate() #HBcurve.Differentiate() metoda daje objekt klase BSpline

HBcurve2 = BSpline(2, [0]+list(B_ref2), list(H_ref2)) #ovo bi trebala biti instanca klase BSpline
diffHB2= HBcurve.Differentiate() #HBcurve.Differentiate() metoda daje objekt klase BSpline

#------------------------------------------
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
#-------------------------------------------

#graddom = [True if mat == "core" else False for mat in mesh.GetMaterials()]
fsU = HCurl(mesh, order=0, dirichlet="gamaE|right|gamaB|front", complex=False, nograds = False)#, gradientdomains = graddom)
#updateHCurlRegionOrder(fsU, 0, "core")
#------------
#fsU = HCurl(mesh, order=0, dirichlet="outer", complex=False, nograds = False)
fsV = H1(mesh, order=1, definedon='core', dirichlet='right|gamaE',complex=False) #bez right!
fes=fsU*fsV
mvp, esp = fes.TrialFunction()
alpha, phi = fes.TestFunction()

print('fes.ndof=',fes.ndof)
print('...free =', sum(fes.FreeDofs()))

#----------------------------------
omega=2*pi*50 #50Hz
dt = 2*pi/omega /51
cnt = 0; t0=0.0 ; time = t0; tend=0.35*(2*pi/omega)
interval=(0.4,1.0); samplestep=2

d=0.00035 #m
Kf= 27*d/0.01 #0.945
#sigma = 2e6
sigma=CF( (1 , 0, 0,   0, 2e6, 0,  0, 0, 2e6), dims=(3,3) )
#rho=CF( (0.1 , 0, 0,   0, 1/(Kf*kappa), 0,  0, 0, 1/(Kf*kappa)), dims=(3,3) )
sig = mesh.MaterialCF({ "core" : 1 }, default=None)

#rot=CF( (0 , 0, 0,   0, 0, 1,  0, -1, 0), dims=(3,3) )
#----------------------------------

gfu = GridFunction(fes)
old = GridFunction(fes)
prev= GridFunction(gfu.space)

Apot, Vpot = gfu.components
oldApot, oldVpot = old.components
prevApot, prevVpot = prev.components

t = Parameter(0.0)
#...
ic=CF((0,0,0)) #initial condition
Apot.Set(ic)
oldApot.Set(ic)
#prev.vec.data=gfu.vec
#..

#inicijalizacija listi:::
time_axis=[]
Hlist=[]
Blist=[]
Jlistoftuple=[]
Plist=[]
Bavglist=[]

gfut = GridFunction(gfu.space,multidim=0) #objekt za spremanje gfu u razlicitim trenutcima t
Jgfut = GridFunction(gfu.space,multidim=0) #objekt za spremanje Jgfu u razlicitim trenutcima t
gfut.AddMultiDimComponent(gfu.vec) #u prvom stupcu je inicijalno stanje gfu(t=0)
Jgfut.AddMultiDimComponent(gfu.vec)

B=curl(oldApot)
Babs= B.Norm() 

volumen=Integrate(CF(1), mesh, definedon=mesh.Materials('core'))

Ts_of_t = sin(omega*t) #vremesnka ovisnost pobude Ts

while time < tend-0.5*dt:
    print('*time*',time)
    t.Set(time)
#    gfuD.Set(bc,definedon=mesh.Boundaries('interface')) #vrem. ovisan bc (t.Set(time))
    errorlist=[]
    prev.vec.data=gfu.vec

    for i in range(1,4):
        print(f"####iteration i={i}")
        
        print('Babs=',Babs(mesh(0.0025,0.0025,0.0025)))
        
        #relyz = (HBcurve(Babs+1e-6))/(Babs+1e-6) + 1j*omega*kappa*d**2*1/12 *(Babs+1e-6)**0.5 #PAZI
        #dHdByz = diffHB(Babs+1e-6) + 1j*omega*kappa*d**2*1/12 *1.5*(Babs+1e-6)**0.2 
        
        #relyz = (HBcurve(Babs+1e-6))/(Babs+1e-6) #+ 1j*(HBcurve2(Babs+1e-6))/(Babs+1e-6) #PAZI
        #dHdByz = diffHB(Babs+1e-6) #+ 1j*diffHB2(Babs+1e-6)

        relyz = 50.0
        dHdByz=50.0001

        """ prevH = HBcurve(Babs+1e-6)
        prevH = 0.1* (HabsPrev +1e-4)**(0.5) * IfPos(curl(prevApot),1,-1)
        oldH = 0.1* (Habs+1e-4)**(0.5) * IfPos(curl(oldApot),1,-1)
        #mud = mesh.MaterialCF({ "core" : dBdH }, default=mu0)  # OSNOVNA DEFINICIJA za mud
        nud=(oldH-prevH)/(curl(oldApot)-curl(prevApot)+1e-4)  # ALTERNATIVNA DEFINICIJA za mud koja se oslanja na vremensku derivaciju """

        rel=CF( ( (1-Kf)/mu0, 0, 0,   0, relyz/Kf, 0,  0, 0, relyz/Kf), dims=(3,3) )
        dHdB=CF( ( (1-Kf)/mu0, 0, 0,   0, dHdByz/Kf, 0,  0, 0, dHdByz/Kf), dims=(3,3) ) #?? should I divide by Kf

        term1= (1/mu0)*curl(mvp)*curl(alpha)*dx('air|coil') + dt*rel*curl(mvp)*curl(alpha)*dx('core') + \
            sigma*mvp*alpha*dx('core') + sigma*grad(esp)*alpha*dx('core')
        term2= sigma*mvp*grad(phi)*dx('core') + sigma*grad(esp)*grad(phi)*dx('core')
        reg_fact=10 #sto manji elementi to veca regularizacijska sigma
        term3= reg_fact*10.0*mvp*alpha*dx('air|coil') + reg_fact*10.0*(esp)*(phi)*dx('core')

        
        jac= dt* (dHdB - rel)*curl(mvp)*curl(alpha)*dx('core')
    
        a = BilinearForm(term1+term2+term3+jac)

        pc = Preconditioner(a, type="jacobi")
        a.Assemble()

        jacmat= BilinearForm(jac)
        jacmat.Assemble()

        prevterm1= sigma*mvp*alpha*dx('core') + sigma*grad(esp)*alpha*dx('core')
        prevterm2= sigma*mvp*grad(phi)*dx('core') + sigma*grad(esp)*grad(phi)*dx('core')
        prevmat=BilinearForm(prevterm1+prevterm2)
        prevmat.Assemble()



        #>>> test=H1(mesh, dirichlet='topleft')
        #>>> gft=GridFunction(test)
        #>>> gft.Set(1,BND)
        #>>> Draw(gft)




        #:::::::::::::::::: source current::::::::::::::::::
        e = 0.0026 #debljina jezgre
        #I=3 #Amp 0.75
        #zavoj=1
        dno=-0.00001
        vrh=3*e
        rin=3*e
        rout=4*e
        R=(x**2 +  y**2)**0.5
        #Js=1.414*I*zavoj/((rout-rin)*(vrh-dno))
        Js=4*1e6
        izvan=CF((0,0,Js*(rout-rin)))
        nula=CF((0,0,0))

        Ts_coil= CF( (0,0,Js*(R-rin)) )
        Ts =IfPos(rin-R, nula, IfPos(rout-R, Ts_coil, izvan)) * IfPos((z-vrh)*(z-dno),0,1)
        Ts_term = Ts_of_t* Ts *curl(alpha) * dx("air|coil",bonus_intorder=5) #+ Ts_of_t*Ts *curl(alpha) * dx("coil") 
        #f += Ts_coil *curl(alpha) * dx("coil") + Ts_air *curl(alpha) * dx("air|core") 
        
        f=LinearForm(Ts_term)
        f.Assemble()
        #::::::::::::::::
        #::::::::::::::::


        #inva=a.mat.Inverse(freedofs=fes.FreeDofs())
        #res = f.vec - a.mat *gfuD.vec + jacmat.mat * oldgfu.vec
        #gfu.vec.data = gfuD.vec + inva * res

        #pc = Preconditioner(a, type="multigrid")
        ##### SOLVER
        #start=time.time()
        r = f.vec + jacmat.mat * oldApot.vec + prevmat.mat * prev.vec #old.vec <=ILI
        #r = f.vec + jacmat.mat * oldApot.vec + prevmat.mat * prevApot.vec #old.vec <=ILI
        #gfu.vec.data = a.mat.Inverse(freedofs=fes.FreeDofs())*r
        
        r_bvp = LinearForm(fes).Assemble()
        r_bvp.vec.data += r
        solvers.BVP(bf=a, lf=r_bvp, gf=gfu, pre=pc, maxsteps=2000, needsassembling=False)
        #stop=time.time()
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
    
    #Hlist.append(gfu[0](mesh(0,0.495*brid)))
    #Blist.append(curl(Apot)[2](mesh(0.0025,0.0025,0.0025)))
    #Jlistoftuple.append(sig*sigma/dt*(prevApot-Apot+ grad(prevVpot)-grad(Vpot))(mesh(0.0025,0.0025,0.0025)))
    time_axis.append(time)

    #current power:
    #grad(gfu)*grad(gfu)==rot*grad(gfu)*rot*grad(gfu)
    currPow=Integrate(sigma/dt**2 *(prevApot-Apot+ grad(prevVpot)-grad(Vpot))* \
                      (prevApot-Apot+ grad(prevVpot)-grad(Vpot)), mesh, definedon=mesh.Materials('core'))
    currBavg=Integrate(Babs/volumen, mesh, definedon=mesh.Materials('core'))
    Plist.append(currPow)
    Bavglist.append(currBavg)
    
    #if time>interval[0]*(tend-t0) and time<interval[1]*(tend-t0) and cnt%samplestep==0 :
    #    gfut.AddMultiDimComponent(gfu.vec)

    #prev.vec.data=gfu.vec

    cnt += 1; time = cnt * dt

#print('errorlist', errorlist)
#####POSTPROCESING

Jpost= sig*sigma/dt*(prevApot-Apot+ grad(prevVpot)-grad(Vpot))
Bpost = curl(Apot)

print('...')


#volumen=Integrate(1,mesh, definedon=defon)

BdV=Integrate(Bpost.Norm(), mesh, order=5, definedon=defon)
print('Bavg=',round(BdV/volumen, 3),'T') 

Draw (Bpost, mesh, "B",draw_surf=False)
Draw (Jpost, mesh, "J",draw_surf=False)
Draw (Ts, mesh, "Ts")

print('...')
print('Bavglist=',Bavglist)
print('Plist=',Plist)
print('time_axis=',time_axis)

