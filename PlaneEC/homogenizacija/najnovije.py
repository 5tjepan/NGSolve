from ngsolve import *
from netgen.occ import *

import netgen.gui

flag=1

d=0.001*40
h=d

size= d/6 #0.0118 #0.0059 #0.00145 #0.00046  #0.00156 #0.0008

core1 = MoveTo(-d/2, -h/2).Rectangle(d,h).Face()
core1.edges.Min(X).name = "l1"
core1.edges.Min(Y).name = "b1"
core1.edges.Max(Y).name = "t1"
core1.edges.Max(X).name = "r1"

core1.faces.name="core1"


#####################################

mesh = Mesh(OCCGeometry(core1, dim=2).GenerateMesh(maxh=size, quad_dominated=False))

fsU = HCurl(mesh, order=0, dirichlet="l1|t1|b1|r1", complex=True, nograds = False)
fsV = H1(mesh, order=1, dirichlet="l1|r1", complex=True)
fes=fsU*fsV
mvp, csp = fes.TrialFunction()
alpha, theta = fes.TestFunction()

print('mesh.nv',mesh.nv)
print('mesh.nedge', mesh.nedge)
print('mesh.ne', mesh.ne)

sol = GridFunction(fes)
Apot, Tpot = sol.components


#sig = mesh.MaterialCF({ "core|Km" : 1 }, default=0.000000001)
#B0=6e2*x+0.4
#B0=6e5*(x**2-d**2/12)+0.4
B0 = 0.4
#Draw(B0,mesh,'B0')

fesl2=L2(mesh, order=0,complex=True)
fesl2cmplx=L2(mesh, order=0,dim=2)

PowLoss=[]
PLref=[]
#PLref=[85.26727719592759, 308.7850278047253, 615.3268050938642, 972.9249103398167, 1371.349833143259, 1808.412872976427, 2283.163375472533, 2794.1953925824864, 3339.664030250664, 3917.5932891902735, 4526.102318750874, 5163.51284362164, 5828.375392190643, 6519.453188119096, 7235.689482097641, 7976.172420612974]

for ind in range(1,17):
    Apot.Set(CF((0,0)))
    Tpot.Set(CF(0))
    omega=314.15*ind/4
    mu0 = 1.257e-6

    rel = 20000 
    rho= 5e-7 

    kh= d/4 #d/2/sqrt(5) #size # d/4

    rot=CF( (0 , 1,  -1, 0), dims=(2,2) )
    diry=CF((0,1))

    #::::::::::::::::::::::::::::::::::::::
    
    #gavg=GridFunction(VectorValued(fesl2))
    gavg=GridFunction(fesl2cmplx)
    favg=GridFunction(fesl2)
    #favgim=GridFunction(fesl2)
    #favgre=GridFunction(fesl2)
    

    arealist=Integrate(1,mesh,element_wise=True)
    intlist=Integrate(x, mesh, element_wise=True)

    #for i in range(len(intlist)):   xavg.vec[i]=intlist[i]/arealist[i]
    
    
    ax=CF((1,0))
    ay=CF((0,1))
    
    zaf= -0.5*(x**2+y**2)
    zag= 0.5*(-y*ax + x*ay)
    
    
    for iter in range(2):
        gavg.Set(zag)
        
        if iter==0:
            #favgre.Set(zaf/sqrt(2))
            #favgim.Set(zaf/sqrt(2))
            favg.Set(zaf)
            aj=gavg/Norm(gavg)
        else:
            #zaf= y*(J*ax)/Norm(J) - x*(J*ay)/Norm(J)
            #favgre.Set(zaf.real)
            #favgim.Set(zaf.imag)
            
            #zaf= y*Norm(J*ax)/Norm(J) - x*Norm(J*ay)/Norm(J)
            zaf= y*(J*ax)/Norm(J) - x*(J*ay)/Norm(J)
            favg.Set(zaf)
            aj=J/Norm(J) #CMPLX vector!


        #Draw(favgre,mesh,"favgre")
        #Draw(favgim,mesh,"favgim")
        Draw(gavg,mesh,"gavg")
        Draw(zaf,mesh,"zaf")
            
        gfun= (zag-1.0*gavg) *(-1j)*omega/rho #
        ffun= (zaf - favg)/rel #  #vjerojatno NIJE minus! ??????
        #ffun= (zaf - favgre - favgim)/rel #  #vjerojatno NIJE minus! ??????

        Draw(gfun, mesh, "gfun")
        Draw(ffun, mesh, "ffun")
        #::::::::::::::::::::::::::::::::::::::


        #term1 = rho*grad(csp)*grad(theta)*dx + 1j*omega*curl(mvp)*theta*dx
        term1 = rho*grad(csp)*grad(theta)*dx + 1j*omega*mvp*(rot*grad(theta))*dx
        #term2 = - 1j*omega* rel*curl(mvp)*curl(alpha)*dx + 1j*omega*csp*curl(alpha)*dx
        term2 = - 1j*omega*rel*curl(mvp)*curl(alpha)*dx + 1j*omega*alpha*(rot*grad(csp))*dx#(bonus_intorder=2)

        #term1a= 1j*omega*theta*ffun* rot*grad(csp)*CF((0,1))*dx(bonus_intorder=2)
        term1a= 1j*omega*theta*ffun* Norm(rot*grad(csp))*dx(bonus_intorder=2)
        #term1a= 1j*omega*theta*ffun* (rot*grad(csp))*aj*dx(bonus_intorder=2)
        term2a= 1j*omega*gfun*alpha*curl(mvp)*dx(bonus_intorder=2)  

        #term1b= 1j*(rho*(gfun*gfun))*curl(alpha)*curl(mvp)*dx(bonus_intorder=2)  
        #term2b= -(1j)*omega*ffun*ffun*(rot*(grad(theta)))*(rot*(grad(csp)))*dx(bonus_intorder=2) 


        a = BilinearForm(1*term1 + 1*term2 + flag*term1a + flag*term2a)
        #a = BilinearForm(1*term1 - 1*term2 - 1*term1a + 1*term2a + 0*term1b + 0*term2b)
        #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


        #.....................
        term1hg= 1j*omega*(1)/(12*rel)*(kh**2) * grad(csp)*grad(theta)*dx
        term2hg= 1/12*(kh*omega)**2 /rho *curl(mvp)*curl(alpha)*dx
        #a = BilinearForm(term1 + term2 + term1hg + term2hg)
        #.....................

        a.Assemble()

        force = -1j * omega*B0*theta*dx
        #force = -1j * omega*B0*theta*dx - rho*gfun*gfun*curl(alpha)*B0*dx - 2*1j*omega*gfun*alpha*B0*dx
        #force = -1j * omega*B0*theta*dx - rho*gfun*gfun*curl(alpha)*B0*dx
        #force = -1j * omega*B0*theta*dx - 1*1j*omega*gfun*alpha*B0*dx
        f=LinearForm(force)
        f.Assemble()

        #print(a.mat)
        r = f.vec - a.mat * sol.vec
        sol.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs()) * r

        #...................
        delta=1/sqrt(omega/(rho*rel*2))
        gama = (1+1j)/delta
        Bfun = B0*gama*d/2 *(cosh(gama*x)/sinh(gama*d/2))
        Jfun = -1j*omega*d/(2*rho)* B0 *(sinh(gama*x)/sinh(gama*d/2))

        #Draw(Bfun, mesh, 'Bfun')

        Pfun=( (h*d)*(B0*omega)**2 *d*delta/(8*rho) * ( (sinh(d/delta)-sin(d/delta))/(cosh(d/delta)-cos(d/delta))) )
        #print('Pfun=',Pfun)
        
        # P  O  S  T

        J = rot*grad(Tpot)
        B=curl(Apot)


    PLref.append(Integrate(Pfun/(h*d),mesh,order=5).real)
    print(f'PLref[{ind}] ={PLref[ind-1]}')
    #Power=0.5*rho*(J + gfun*(B0+B)) * Conj(J + gfun*(B0+B))

    Power=0.5*rho*J*Conj(J) + flag*0.5*rho*(gfun*(B0+B))*Conj((B0+B)*gfun)
    #Power=0.5*rho*J*Conj(J) + 0.5* 1/12*(kh*omega)**2 /rho * (B0+B)*Conj(B0+B)

    PowLoss.append(Integrate(Power,mesh,order=5).real)
    print(f"PowLoss[{ind}] = {PowLoss[ind-1]}")
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #print(Integrate(1,mesh,element_wise = True))
    print('%ER =', round((PowLoss[ind-1] - PLref[ind-1])/PLref[ind-1]*100,1),'%')



#-------------------------------------
f_Km = IfPos( (x-d/2)*(x), 0, IfPos( (y-4*h/12)*(y-h/2), 0, 1))
f_core = mesh.MaterialCF({ "core" : 1 }, default=0.0)
reg_Km =f_Km
#pozivom metode mesh.Materials() dobiva se objekt klase ngsolve.comp.Region
defon = mesh.Materials('Km')
#int_Km = Integrate(CF(1)*reg_Km, mesh, definedon=defon)
int_Km = Integrate(CF(1)*reg_Km, mesh)

B_bar= Integrate((B+B0)*reg_Km, mesh)/int_Km
J_bar= Integrate(J*CF((0,1))*reg_Km, mesh)/int_Km

print('B_bar = ',B_bar)
print('J_bar = ', J_bar)


Draw(Tpot, mesh, "T")
Draw(Apot, mesh, "A")
Draw (B, mesh, "B")
Draw (J, mesh, "J")
#Draw(grad(Tpot),mesh,"gradT")


Jfuntilda=Jfun-J_bar#*CF((0,1))
Jtilda=J-J_bar*CF((0,1))
ggint= Integrate(Jfuntilda*Jfuntilda*reg_Km, mesh)/int_Km/B0**2
#print('ggint=',ggint)
#print('rho*ggint=',rho*ggint)

Pow_ggint=Integrate(Jfuntilda*Conj(Jfuntilda)*reg_Km, mesh)/int_Km/B0**2
#print('Pow_ggint=',Pow_ggint)


""" import numpy as np
import matplotlib.pyplot as plt
import sys
sys.argv = ["dummy"]

X = np.linspace(-0.5*d, 0.5*d, num=100)
Y = np.zeros_like(X)+0.00065

plt.plot(X, gfun.imag(mesh(X, Y)))
plt.plot(X, gfun.imag(mesh(X, Y)))
plt.xlabel('x')
#plt.show() """




