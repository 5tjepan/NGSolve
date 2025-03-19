#from netgen.occ import *
import ngsolve
from netgen.read_gmsh import ReadGmsh
from netgen.meshing import *
from ngsolve import *

# Učitavanje Gmsh mreže
mesh_g = ReadGmsh("gmsh_BoxCoilCore_full/box_coil_box61k.msh")

#mesh = Mesh(mesh)
#mesh.ngmesh.Save("disk.vol")
#Draw(mesh)

mesh = ngsolve.Mesh(mesh_g)

Draw(mesh)

print(mesh.GetMaterials())
print(mesh.GetBoundaries())


fsU = HCurl(mesh, order=0, dirichlet="gamaB", complex=True, nograds = False)
fsV = H1(mesh, order=1, dirichlet='tblr', definedon='core', complex=True) #bez right!
fes=fsU*fsV
mvp, csp = fes.TrialFunction()
alpha, tau = fes.TestFunction()

print('fes.ndof=',fes.ndof)
print('...free =', sum(fes.FreeDofs()))

gfu = GridFunction(fes)
old = GridFunction(fes)

Apot, Tpot = gfu.components
oldApot, oldTpot = old.components

mu0=4*pi*1e-7
omega=314
d=0.00035 #m
Kf= 27*d/0.01 #0.945
kappa = 2e5
rho=CF( (1 , 0, 0,   0, 1/(Kf*kappa), 0,  0, 0, 1/(Kf*kappa)), dims=(3,3) )
sig = mesh.MaterialCF({ "core" : 1 }, default=None)

rot=CF( (0 , 0, 0,   0, 0, 1,  0, -1, 0), dims=(3,3) )

B=curl(oldApot)
Babs= B.Norm() 
errorlist=[]
p=1.0


for i in range(1,2):
    print(f"####iteration i={i}")
    
    #print('Babs=',Babs(mesh(0,0)))
    
    relyz= 1/(500*mu0) #+ 1j*omega*kappa*d**2*1/12 #pazi, kappa_y doprinosi rel_z i obratno  
    #dHdByz= 1/(500*mu0) #+ 1j*omega*kappa*d**2*1/12

    #relyz = (HBcurve(Babs+1e-6))/(Babs+1e-6) + 1j*omega*kappa*d**2*1/12 *(Babs+1e-6)**0.5 #PAZI
    #dHdByz = diffHB(Babs+1e-6) + 1j*omega*kappa*d**2*1/12 *1.5*(Babs+1e-6)**0.2 
    
    #relyz = (HBcurve(Babs+1e-6))/(Babs+1e-6) #+ 1j*(HBcurve2(Babs+1e-6))/(Babs+1e-6) #PAZI
    #dHdByz = diffHB(Babs+1e-6) #+ 1j*diffHB2(Babs+1e-6) 

    rel=CF( ( (1-Kf)/mu0, 0, 0,   0, relyz/Kf, 0,  0, 0, relyz/Kf), dims=(3,3) )
    #dHdB=CF( ( (1-Kf)/mu0, 0, 0,   0, dHdByz/Kf, 0,  0, 0, dHdByz/Kf), dims=(3,3) ) #?? should I divide by Kf

    term1=(1/mu0)*curl(mvp)*curl(alpha)*dx('air|coil') + rel*curl(mvp)*curl(alpha)*dx('core') \
    - (rot*grad(csp))*alpha*dx('core') #+ 0.1*mvp*alpha*dx('air|coil')
    term2= 1j/omega*rho*(rot*grad(csp))*(rot*grad(tau))*dx('core') - (rot*grad(tau))*mvp*dx('core')
    term3=1*1e0*mvp*alpha*dx

    #jac= (dHdB - rel)*curl(mvp)*curl(alpha)*dx('core')
  
    a = BilinearForm(term1+term2+term3) #+jac)
    a.Assemble()

    #jacmat= BilinearForm(jac)
    #jacmat.Assemble()

    #:::::::::::::::::: source current
    #f = LinearForm(fes)
    I=1.2*(1j+1)/1.414213562373 #Amp 0.75
    zavoj=447
    dno=-0.1
    vrh=0.1
    rin=0.12
    rout=0.2
    R=(x**2 +  y**2)**0.5
    Js=1.414*I*zavoj/((rout-rin)*(vrh-dno))
    izvan=CF((0,0,Js*(rout-rin)))
    nula=CF((0,0,0))

    Ts_coil= CF( (0,0,Js*(R-rin)) )
    Ts_air=IfPos(rin-R, nula, IfPos(rout-R, Ts_coil, izvan)) * IfPos((z-vrh)*(z-dno),0,1)
    Ts_total = Ts_coil *curl(alpha) * dx("coil") + Ts_air *curl(alpha) * dx("air|core") 
    #f += Ts_coil *curl(alpha) * dx("coil") + Ts_air *curl(alpha) * dx("air|core") 
    
    f=LinearForm(Ts_total)
    f.Assemble()
    #::::::::::::::::

    #pc = Preconditioner(a, type="multigrid")
    ##### SOLVER
    r = f.vec #+ jacmat.mat * oldApot.vec
    #gfu.vec.data = a.mat.Inverse(freedofs=fes.FreeDofs())*r
    
    r_bvp = LinearForm(fes).Assemble()
    r_bvp.vec.data += r
    solvers.BVP(bf=a, lf=r_bvp, gf=gfu, pre=None, maxsteps=2000, print=True, needsassembling=False)
    #----------------------


    """ rhs_jac =  jacmat.mat * oldApot.vec #old.vec <=ILI
    r=LinearForm(Ts_total+rhs_jac)
    r.Assemble()
    solvers.BVP(bf=a, lf=r, gf=gfu, pre=None, maxsteps=2000, print=True) """
    
    #errfunc = (gfu - old).Norm()/(old.Norm()+1e-15)
    #errfunc = (Apot - oldApot).Norm()/(oldApot.Norm()+1e-15)
    
    defon = mesh.Materials('core')
    #error=Integrate(errfunc.Norm(), mesh) #, definedon=defon)
    #print('error =', error)
    #errorlist.append(error)

    #old.vec.data= p*gfu.vec + (1-p)*old.vec
    #old.vec.data= gfu.vec
    
    B=curl(Apot)
    Babs= B.Norm()

#print('errorlist', errorlist)
#####POSTPROCESING

Jpost= rot*grad(Tpot)
Bpost = curl(Apot)

Draw (Bpost, mesh, "B", draw_surf=False)
Draw (Jpost, mesh, "J", draw_surf=False)
#Draw (Ts, mesh, "Ts")

p_wide=0.5*rho*Jpost*Conj(Jpost) 
P_yz=Integrate(p_wide, mesh, order=5, definedon=defon)
print('P_yz=', 1e3*round(abs(P_yz),5),'mW')


volumen=Integrate(1,mesh, definedon=defon)
BdV=Integrate(Bpost.Norm(), mesh, order=5, definedon=defon)
print('Bavg=',round(BdV/volumen, 4),'T') 

BxdV=Integrate(Bpost[0].Norm(), mesh, order=5, definedon=defon)
print('Bx_avg=',round(BxdV/volumen, 8),'T') 

Btang= (Bpost[1].Norm()**2+Bpost[2].Norm()**2)
BdV=Integrate(Btang, mesh, order=5, definedon=defon)
print('Btang**2=',round(BdV/volumen, 6),'T') 

