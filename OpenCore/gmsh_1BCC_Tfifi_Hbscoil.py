#from netgen.occ import *
import ngsolve
from netgen.read_gmsh import ReadGmsh
from netgen.meshing import *
from ngsolve import *

# Učitavanje Gmsh mreže
mesh_g = ReadGmsh("box_coil_box.msh")

#mesh = Mesh(mesh)
#mesh.ngmesh.Save("disk.vol")
#Draw(mesh)

mesh = ngsolve.Mesh(mesh_g)

Draw(mesh)

print(mesh.GetMaterials())
print(mesh.GetBoundaries())


#========================

fesH0 = HCurl(mesh, order=0, dirichlet="gamaB", complex=True, nograds = False)
mvp = fesH0.TrialFunction()
alpha = fesH0.TestFunction()

gfsource = GridFunction(fesH0)

#:::::::::::::::::: source current
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
#Ts_air=IfPos(rin-R, nula, IfPos(rout-R, Ts_coil, izvan)) * IfPos((z-vrh)*(z-dno),0,1)
Ts=IfPos(rin-R, nula, IfPos(rout-R, Ts_coil, izvan)) * IfPos((z-vrh)*(z-dno),0,1)

part1= curl(mvp)*curl(alpha)*dx 
stiff = BilinearForm(part1)
stiff.Assemble()

fs=LinearForm(fesH0)
fs +=Ts *curl(alpha)*dx
fs.Assemble()
#::::::::::::::::

solvers.BVP(bf=stiff, lf=fs, gf=gfsource, pre=None, maxsteps=2000, print=True, needsassembling=False)
#----------------------

Hbs=curl(gfsource)

#Draw(Hbs,mesh,'Hbs')
#||||||||||||||||||||||||


fsU = HCurl(mesh, order=0, dirichlet="tblr|front_back", definedon='core', complex=True, nograds = False)
fsV = H1(mesh, order=1, complex=True) #bez right!
fes=fsU*fsV
cvp, msp = fes.TrialFunction()
tau, psi = fes.TestFunction()

print('fes.ndof=',fes.ndof)
print('...free =', sum(fes.FreeDofs()))

gfu = GridFunction(fes)
#old = GridFunction(fes)

Tpot, Fpot = gfu.components
#oldTpot, oldFpot = old.components

mu0=4*pi*1e-7
omega=314
d=0.00035 #m
Kf= 27*d/0.01 #0.945
kappa = 2e3
rho=CF( (0.1 , 0, 0,   0, 1/(Kf*kappa), 0,  0, 0, 1/(Kf*kappa)), dims=(3,3) )
sig = mesh.MaterialCF({ "core" : 1 }, default=None)

rot=CF( (0 , 0, 0,   0, 0, 1,  0, -1, 0), dims=(3,3) )

#B=curl(oldApot)
#Babs= B.Norm() 
errorlist=[]
p=1.0


for i in range(1,2):
    print(f"####iteration i={i}")
    
    #print('Babs=',Babs(mesh(0,0)))
    
    #relyz= 500 + 1j*omega*kappa*d**2*1/12 #pazi, kappa_y doprinosi rel_z i obratno  
    #dHdByz= 500 + 1j*omega*kappa*d**2*1/12

    #relyz = (HBcurve(Babs+1e-6))/(Babs+1e-6) + 1j*omega*kappa*d**2*1/12 *(Babs+1e-6)**0.5 #PAZI
    #dHdByz = diffHB(Babs+1e-6) + 1j*omega*kappa*d**2*1/12 *1.5*(Babs+1e-6)**0.2 
    
    #relyz = (HBcurve(Babs+1e-6))/(Babs+1e-6) #+ 1j*(HBcurve2(Babs+1e-6))/(Babs+1e-6) #PAZI
    #dHdByz = diffHB(Babs+1e-6) #+ 1j*diffHB2(Babs+1e-6) 
    muyz=500*mu0
    mucore=CF( ( mu0/(1-Kf), 0, 0,  0, muyz*Kf, 0,  0, 0, muyz*Kf), dims=(3,3) )
    muair=CF( ( mu0, 0, 0,  0, mu0, 0,  0, 0, mu0), dims=(3,3) )
    
    mu= mesh.MaterialCF({ "core" : mucore }, default=muair)
    #dHdB=CF( ( (1-Kf)/mu0, 0, 0,  0, dHdByz/Kf, 0,  0, 0, dHdByz/Kf), dims=(3,3) ) #?? should I divide by Kf

    term1= rho*curl(cvp)*curl(tau)*dx('core') + 1j*omega*mucore*cvp*tau*dx('core') \
    - 1j*omega*mucore*grad(msp)*tau*dx('core') #+ 0.1*mvp*alpha*dx('air|coil')
    term2= -1j*omega*mucore*cvp*grad(psi)*dx('core') + 1j*omega*mucore*grad(msp)*grad(psi)*dx('core') \
    + 1j*omega*muair*grad(msp)*grad(psi)*dx('air|coil')
    #term3=1*1e0*mvp*alpha*dx

    #jac= (dHdB - rel)*curl(mvp)*curl(alpha)*dx('core')
    a = BilinearForm(term1+term2) #+jac)
    a.Assemble()

    #jacmat= BilinearForm(jac)
    #jacmat.Assemble()


    #Ts_total = -1j*omega*mu*Hbs *tau*dx('core') + 1j*omega*mu*Hbs*grad(psi)*dx 
    
    #Ts_total = -1j*omega*mu*Ts_coil *tau * dx("coil") + Ts_air *curl(alpha) * dx("air|core") 
    #f += Ts_coil *curl(alpha) * dx("coil") + Ts_air *curl(alpha) * dx("air|core") 
    
    f=LinearForm(fes)
    f += -1j*omega*mucore*Hbs *tau*dx('core') + 1j*omega*mucore*Hbs*grad(psi)*dx('core') \
    -1j*omega*muair*Hbs *tau*dx('air|coil') + 1j*omega*muair*Hbs*grad(psi)*dx('air|coil') 
    
    f.Assemble()
    #::::::::::::::::

    #pc = Preconditioner(a, type="multigrid")
    ##### SOLVER
    #r = f.vec #+ jacmat.mat * oldApot.vec
    ###gfu.vec.data = a.mat.Inverse(freedofs=fes.FreeDofs())*r
    
    #r_bvp = LinearForm(fes).Assemble()
    #r_bvp.vec.data += r
    solvers.BVP(bf=a, lf=f, gf=gfu, pre=None, maxsteps=2000, print=True, needsassembling=False)
    #----------------------


    #errfunc = (gfu - old).Norm()/(old.Norm()+1e-15)
    #errfunc = (Apot - oldApot).Norm()/(oldApot.Norm()+1e-15)
    
    defon = mesh.Materials('core')
    #error=Integrate(errfunc.Norm(), mesh) #, definedon=defon)
    #print('error =', error)
    #errorlist.append(error)

    #old.vec.data= p*gfu.vec + (1-p)*old.vec
    #old.vec.data= gfu.vec
    
    #B=mu*(Tpot - grad(Fpot) + Hbs)
    #Babs= B.Norm()

#print('errorlist', errorlist)
#####POSTPROCESING

Jpost=curl(Tpot)

sig_air = mesh.MaterialCF({ "air|coil" : 1 }, default=None)
Bcore =mucore*(Tpot - grad(Fpot) + Hbs) *sig
Bair =muair*(Tpot - grad(Fpot) + Hbs) * sig_air


Draw (Bcore, mesh, "B")
Draw (Bair, mesh, "Bair")
Draw (Jpost, mesh, "J")

Draw(Hbs,mesh,'Hbs')
Draw (Fpot, mesh, "Fpot")
Draw (grad(Fpot), mesh, "grad(Fpot)")



p_wide=0.5*rho*Jpost*Conj(Jpost) 
P_yz=Integrate(p_wide, mesh, order=5, definedon=defon)
print('P_yz=', 1e3*round(abs(P_yz),4),'mW')
