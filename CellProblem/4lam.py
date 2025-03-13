# 2D AT formulacija 
# 3lamele 
# nonlinear 
# pobuda su B0 i J0


from ngsolve import *
from netgen.occ import *
import netgen.gui
import matplotlib.pyplot as plt
import sys
sys.argv = ["fun"]

outer= Circle((0,0), 0.2).Face()
#outer.edges.name = 'rub'

#inner = MoveTo(-0.08,-0.08).Line(0.16,0.0).Line(0,0.16).Line(-0.16,0).Close().Face()
#inner = MoveTo(-0.8,-0.8).Line(1.6,0).Line(0,1.6).Line(-1.6,0).Close().Face()


width=0.002
d=width/4
h=4*d
izOVERd=1/25
fill_fact=4/(4+3*izOVERd)
lam1 = MoveTo(0,-h).Line(d,0.0).Line(0,h).Line(-d,0).Close().Face()
lam1.edges.Min(X).name = 'l1'
lam1.edges.Max(X).name = 'r1'
lam1.edges.Min(Y).name = 'b1'
lam1.edges.Max(Y).name = 't1'
#lam1.edges.name="interface1"
#lam1.faces.maxh=0.01
lam1.faces.name="lam1"
lam1.faces.col = (1, 1, 0)  #colour

lam2 = MoveTo(d,-h).Line(d*izOVERd,0.0).Line(0,h).Line(-d*izOVERd,0).Close().Face()
#lam2.edges.Min(X).name = 'l2'
lam2.edges.Max(X).name = 'r2'
lam2.edges.Min(Y).name = 'b2'
lam2.edges.Max(Y).name = 't2'
#lam2.faces.maxh=0.01
lam2.faces.name="lam2"
lam2.faces.col = (1, 0.5, 0)  #colour

lam3 = MoveTo((1+izOVERd)*d,-h).Line(d,0.0).Line(0,h).Line(-d,0).Close().Face()
#lam3.edges.Min(X).name = 'l3'
lam3.edges.Max(X).name = 'r3'
lam3.edges.Min(Y).name = 'b3'
lam3.edges.Max(Y).name = 't3'
#lam3.faces.maxh=0.01
lam3.faces.name="lam3"
lam3.faces.col = (1, 0.1, 0)  #colour


lam4 = MoveTo((2+izOVERd)*d,-h).Line(d*izOVERd,0.0).Line(0,h).Line(-d*izOVERd,0).Close().Face()
#lam4.edges.Min(X).name = 'l4'
lam4.edges.Max(X).name = 'r4'
lam4.edges.Min(Y).name = 'b4'
lam4.edges.Max(Y).name = 't4'
#lam4.faces.maxh=0.01
lam4.faces.name="lam4"
lam4.faces.col = (1, 0.5, 0)  #colour

lam5 = MoveTo((2+2*izOVERd)*d,-h).Line(d,0.0).Line(0,h).Line(-d,0).Close().Face()
#lam5.edges.Min(X).name = 'l5'
lam5.edges.Max(X).name = 'r5'
lam5.edges.Min(Y).name = 'b5'
lam5.edges.Max(Y).name = 't5'
#lam5.faces.maxh=0.01
lam5.faces.name="lam5"
lam5.faces.col = (1, 0.1, 0)  #colour


lam6 = MoveTo((3+2*izOVERd)*d,-h).Line(d*izOVERd,0.0).Line(0,h).Line(-d*izOVERd,0).Close().Face()
#lam6.edges.Min(X).name = 'l6'
lam6.edges.Max(X).name = 'r6'
lam6.edges.Min(Y).name = 'b6'
lam6.edges.Max(Y).name = 't6'
#lam6.faces.maxh=0.01
lam6.faces.name="lam6"
lam6.faces.col = (1, 0.5, 0)  #colour

lam7 = MoveTo((3+3*izOVERd)*d,-h).Line(d,0.0).Line(0,h).Line(-d,0).Close().Face()
#lam7.edges.Min(X).name = 'l7'
lam7.edges.Max(X).name = 'r7'
lam7.edges.Min(Y).name = 'b7'
lam7.edges.Max(Y).name = 't7'
#lam7.faces.maxh=0.01
lam7.faces.name="lam7"
lam7.faces.col = (1, 0.1, 0)  #colour


rub='l1|t1|t2|t3|t4|t5|t6|t7|r7|b7|b6|b5|b4|b3|b2|b1'
svi='l1|t1|t2|t3|t4|t5|t6|t7|r7|b7|b6|b5|b4|b3|b2|b1|r1|r2|r3|r4|r5|r6'
A_bnd= rub # A_bnd='' <- bez rubnih uvjeta na A
top='l1|r1|r2|r3|r4|r5|r6|r7|t1|t2|t3|t4|t5|t6|t7'
inf='l1|r1|r2|r3|r4|r5|r6|r7'
T_bnd=top

geo = Glue([lam1,lam2,lam3,lam4,lam5,lam6,lam7])
lams='lam1|lam3|lam5|lam7'
izol='lam2|lam4|lam6'
all='lam1|lam3|lam5|lam7|lam2|lam4|lam6'

mesh = Mesh(OCCGeometry(geo, dim=2).GenerateMesh(maxh=0.000028, quad_dominated=False))
print('mesh.nv',mesh.nv)
print('mesh.nedge', mesh.nedge)
print('mesh.ne', mesh.ne)
print(mesh.GetBoundaries())

########## VIZUALIZACIJA
mu0=4*pi*1e-7

B_ref=[0.0001, 0.04, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.83, 2, 2.5, 10000]
H_ref0=[0.01, 3.5, 5.0, 7.6, 10.0, 12.07, 14.08, 16.0, 17.75, 19.42, 21.05, 22.64, 24.16, 25.76, 27.98, 29.66, 32.05, 35.33, 40.06, 47.17, 58.35, 77.19, 113.0, 196.0, 395.0, 800, 8000, 68000, 10000/mu0]
H_ref=[elem*1.8 for elem in H_ref0]

HBcurve = BSpline(3, [0,0]+list(B_ref), list(H_ref)) #ovo bi trebala biti instanca klase BSpline
diffHB= HBcurve.Differentiate() #HBcurve.Differentiate() metoda daje objekt klase BSpline

#munonlin = BSpline(2, [0]+list(B_ref), list(H_ref)) #ovo bi trebala biti instanca klase BSpline
#muder= munonlin.Differentiate() #munonlin.Differentiate() metoda daje objekt klase BSpline
##################

fsU = HCurl(mesh, order=0, dirichlet=A_bnd,  complex=True, nograds = False) #CMPLX
fsV = H1(mesh, order=1, dirichlet=T_bnd, definedon=lams, complex=True)
fes=fsU*fsV
mvp, csp = fes.TrialFunction()
alpha, theta = fes.TestFunction()

#bb=20*CF((0.01*y,-0.01*x))

gfu = GridFunction(fes)
old = GridFunction(fes)

Apot, Tpot = gfu.components
oldApot, oldTpot = old.components

Apot.Set(CF((0,0)), definedon=mesh.Boundaries(rub))
oldApot.Set(CF((0,0)), definedon=mesh.Boundaries(rub))

#Apot.Set(bb, definedon=mesh.Boundaries(rub))
#oldApot.Set(bb, definedon=mesh.Boundaries(rub))

omega=314*3
sigma = 2e6
rho=1/sigma
sig = mesh.MaterialCF({ lams : 1 }, default=None)
feromag = mesh.MaterialCF({ izol : 1/1000 }, default=1)

B0= 1.0/fill_fact*feromag #+ 50*(x - width/2)
J0= 0.01*CF((0,1)) #9e6*(x-width/2)*CF((0,1))

B=curl(oldApot)
#B=B0
Babs= B.Norm() 
errorlist=[]
p=1

nonlin= True

for i in range(1,7):
    print(f"####iteration i={i}")
    
    print('Babs=',Babs(mesh(0,0)))
    print('B =',B(mesh(0,0)))
    
    #Apot.Set(CF((0,0)), definedon=mesh.Boundaries(rub))
    #Tpot.Set(CF(0), definedon=mesh.Boundaries(T_bnd))
    
    #Apot.Set(bb, definedon=mesh.Boundaries(rub))
    
    #Apot.Set(bb, definedon=mesh.Boundaries('rub|interface1'))
    #oldApot.Set(bb, definedon=mesh.Boundaries('rub'))

    if nonlin:
        rel = (HBcurve(Babs+1e-3))/(Babs+1e-3)
        dHdB = diffHB(Babs+1e-3) 
        #rel= Babs*100 + 1e-5
        #dHdB= 2*Babs*100 + 1e-5
        #rel= 100.0001 * Babs**0.5 + 1e-5  #63.257
        #dHdB= 100.0 *1.5*Babs**0.5 + 1e-5  #63.258
    else:
        rel= 64.48 #83.45
        dHdB= 64.481 #83.451
    
    #print('rel=', rel(mesh(0.6,0.6)))
    #print('dHdB=', dHdB(mesh(0.6,0.6)))

    #term1 = (1/mu0)*curl(alpha)*curl(mvp)*dx('outer') + rel*curl(alpha)*curl(mvp)*dx('lam1') + \
    # 1j*omega*sigma*alpha*mvp*dx('lam1') + 1j*1e-4*alpha*mvp*dx('outer')
    
    term1 = rho * grad(csp)*grad(theta)*dx(lams) + 1j*omega*curl(mvp)*theta*dx(lams) + 1j*1*mvp*alpha*dx
    term2 = -rel*curl(mvp)*curl(alpha)*dx(lams) + csp*curl(alpha)*dx(lams) - 30000*curl(mvp)*curl(alpha)*dx(izol)
    jac= -(dHdB - rel)*curl(alpha)*curl(mvp)*dx(lams)
    #term2 = - 1j*omega* rel*curl(mvp)*curl(alpha)*dx(lams) + 1j*omega* csp*curl(alpha)*dx(lams) - 1j*omega*30000*curl(mvp)*curl(alpha)*dx(izol)
    #jac= -1j*omega*(dHdB - rel)*curl(alpha)*curl(mvp)*dx(lams)


    a = BilinearForm(term1+term2 + jac)
    a.Assemble()

    jacmat= BilinearForm(jac)
    jacmat.Assemble()

    force = -1j * omega*B0*theta*dx #(lams) -1j * omega*J0*alpha*dx(lams)
    f=LinearForm(force)
    f.Assemble()

    """     dummy=CF((1e-17,1e-17))
    #rhs1 = dummy*v*dx
    #rhs1 = (dHdB - rel)*(Conj(B)*curl(old))*(B*curl(v))/(Babs*Babs)*dx('lam1')
    rhs1 = (dHdB - rel)*curl(oldApot)*curl(alpha) *dx(lams)
    f = LinearForm(rhs1)
    f.Assemble() """

    dirich=[]
    
    for i in range(len(fes.FreeDofs())):
        if fes.FreeDofs()[i]: dirich.append(gfu.vec[i])
    print('before dirich = ', dirich[:3]) 
    #if i<=5: print('before gfu.vec=', gfulist[:])

    ##### SOLVER
    #solvers.BVP(bf=a, lf=f, gf=gfu, pre=None, maxsteps=2000, print=True)
    
    r = f.vec + jacmat.mat * old.vec #- a.mat * gfu.vec
    #r = f.vec - a.mat * gfu.vec
    gfu.vec.data = a.mat.Inverse(freedofs=fes.FreeDofs())*r  #PAZI koristim '=' umjesto '+='

    dirich=[]
    for i in range(len(fes.FreeDofs())):
        if fes.FreeDofs()[i]: dirich.append(gfu.vec[i])
    print(' after dirich = ', dirich[:3]) 

    errfunc = (Apot - oldApot).Norm()/oldApot.Norm()
    defon = mesh.Materials(lams)
    error=Integrate(errfunc.Norm(), mesh, definedon=defon)
    print('error =', error)
    errorlist.append(error)

    old.vec.data= gfu.vec
    
    B=curl(oldApot)+B0
    
    Babs= B.Norm()

print('errorlist', errorlist)

#>>>>>>>>>>>>>>>>>>
#####POSTPROCESING

potA=oldApot
rot=CF( (0 , 1,  -1, 0), dims=(2,2) )
J = rot*grad(oldTpot) +J0
B = (curl(potA) + B0)*sig

Draw(potA, mesh, "A")
Draw (B, mesh, "B")
Draw (J, mesh, "J")
#B0_fig= 1.2/fill_fact + 50*(x - width/2)
#Draw (B0_fig, mesh, "B0_fig")
#Draw (rel, mesh, "rel")

Pow=0.5*rho*J*Conj(J)
print('PowLoss =', Integrate(Pow,mesh,order=5).real)

#defon = mesh.Materials(lams)
defon = mesh.Materials(all)
Peddy=Integrate(Pow, mesh, order=5)
print('Peddy',Peddy)


int_Km = Integrate(CF(1), mesh, definedon=defon)

mu_int= Integrate(1/(rel+1e-6), mesh, definedon=defon)/int_Km
print('avg of mu_rel=',mu_int/mu0)

b_bar= Integrate(B, mesh, definedon=defon)/int_Km
print('b_bar=',b_bar)
b_sqr=B*B
#int_bsqr = Integrate(b_sqr, mesh, definedon=defon)
int_bsqr = Integrate(b_sqr*rel, mesh, definedon=defon)
FF_avg=(1/int_Km * int_bsqr/(b_bar**2))
print('FF=1/Km * int_bsqr/bbar^2= rel* =',FF_avg)

#j_sqr=(J*Conj(J)).Norm()
j_sqr=J*J
int_jsqr = Integrate(j_sqr, mesh, definedon=defon)
GG_avg =(1/int_Km * int_jsqr/(b_bar**2))
print('GG=1/Km * int_jsqr/jbar^2= rho* =',GG_avg)


#nuzz=1/mu*FF_avg - 1j/omega/sigma *GG_avg
nuzz=FF_avg - 1j/omega/sigma *GG_avg
print('nuzz',nuzz)


jabs_sqr=(J.Norm())**2
int_jabs_sqr = Integrate(jabs_sqr, mesh, definedon=defon)
GGabs_avg =(1/int_Km * int_jabs_sqr/abs(b_bar)**2)

Pcoeff=0.5/sigma * GGabs_avg #(1/int_Km * abs(int_jsqr)/((abs(b_bar))**2))
print('Pcoeff =', Pcoeff)
print('P_Km = Pcoeff*b_bar^2 * int_Km=', Pcoeff*abs(b_bar)**2 * int_Km)

area= Integrate(feromag, mesh)
print('area',area)
b_bar_fe=Integrate(B, mesh)/area
print('Bavg=',b_bar_fe)

if nonlin:
    mu=abs(b_bar_fe)/HBcurve(CF(abs(b_bar_fe)))(mesh(0,0)) #mesh(0,0) bi trebao postojati
    #mu = 1/(100.0 * b_bar_fe.real**0.5 + 1e-5)
    #mu = 1/(100.0 * b_bar_fe.real + 1e-5)
else:
    mu=1/rel #mesh(0,0) bi trebao postojati
print('mu_rel=',mu/mu0)

delta = sqrt(2/(omega*sigma*mu))

F_B = sigma*d*1j*omega/(8*(sinh((1+1j)*d/2/delta))**2) * (delta*(1-1j)/2*sinh((1+1j)*d/delta)+ d)
F_J = sigma*d*1j*omega/(8*(sinh((1+1j)*d/2/delta))**2) * (delta*(1-1j)/2*sinh((1+1j)*d/delta)- d)

F_R = F_B + F_J

print('F_B=',F_B)
print('F_J=',F_J)
print('F_R=',F_R)


F_R_bulk= 1/(fill_fact/F_R + (1-fill_fact)*mu0)
print('F_R_bulk=',F_R_bulk)

F_pow =sigma*omega**2*delta*d/8*(sinh(d/delta)-sin(d/delta))/(cosh(d/delta)-cos(d/delta))/fill_fact
print('F_pow', F_pow)

print('delta/d=',delta/d)


""" import numpy as np
X=np.linspace(0,0.003,100)
Y=np.ones_like(X)*(-0.0005)
B(mesh(X,Y)) 
plt.plot(X,B(mesh(X,Y)))
plt.show()"""

