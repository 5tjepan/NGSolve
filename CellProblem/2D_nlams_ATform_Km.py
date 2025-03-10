
from ngsolve import *
import netgen.gui
from netgen.occ import *

#air= Circle((0,0), 0.008).Face()

size=0.000092
n=9
d=0.001
ins=0.0001
#fill_fact=n*d/(n*d+(n-1)*ins)
fill_fact=3*d/(3*d + 2.5*ins)

w=(n*d+(n-1)*ins)/2/3
h=(n*d+(n-1)*ins)/2/3
n=3
air = MoveTo(-w, -h).Rectangle(2*w, 2*h).Face()
air.edges.name = 'rub'
air.faces.name="air"
air.faces.col = (1, 1, 0)

""" lamele = MoveTo(-w,-h).Line(d,0.0).Line(0,2*h).Line(-d,0).Close().Face()
for i in range(n-1):
    w-=(d+ins)
    lamele += MoveTo(-w,-h).Line(d,0.0).Line(0,2*h).Line(-d,0).Close().Face() """

lamele=[0]*n
#lamele = MoveTo(-w,-h).Rectangle(d,2*h).Face()
for i in range(n):
    lamele[i]= MoveTo(-w,-h).Rectangle(d,2*h).Face()
    lamele[i].faces.name=f'lamele{i}'
    lamele[i].edges.Min(X).name = f'l{i}'
    lamele[i].edges.Max(X).name = f'r{i}'
    lamele[i].edges.Min(Y).name = f'b{i}'
    lamele[i].edges.Max(Y).name = f't{i}'
    #air-=lamele[i]
    w-=(d+ins)
    
#print('lamele',lamele)
#lamele.faces.col = (0,0,1)
#lamele.faces.maxh=size
""" lamele.edges.Max(X).name = "r"
lamele.edges.Min(X).name = "l"
lamele.edges.Min(Y).name = "b"
lamele.edges.Max(Y).name = "t" """

izol=[air]
#core=lamele.append(air)
#core=Glue(lamele)
geo = Glue(lamele + izol)
mesh = Mesh(OCCGeometry(geo, dim=2).GenerateMesh(maxh=size, quad_dominated=True))

print(mesh.GetBoundaries())

##########################
##########################

#Tpot_dirich_top ="t0|t1|t2"
#Tpot_dirich="l0|r0|l1|r1|l2|r2|t0|t1|t2|rub"
Tpot_dirich="l0|r0|l1|r1|l2|r2|rub"
#Tpot_dirich+=Tpot_dirich_top
#fsU = HCurl(mesh, order=0, complex=True, nograds = False)
fsU = HCurl(mesh, order=0, dirichlet="l0|t0|t1|t2|r2|b2|b1|b0|rub", complex=True, nograds = False)
fsV = H1(mesh, order=1, dirichlet=Tpot_dirich, complex=True)
fes=fsU*fsV
mvp, csp = fes.TrialFunction()
alpha, theta = fes.TestFunction()

print('mesh.nv',mesh.nv)
print('mesh.nedge', mesh.nedge)
print('mesh.ne', mesh.ne)

sol = GridFunction(fes)
Apot, Tpot = sol.components

#bb=CF((0.1*y,-0.1*x))
#Apot.Set(bb, definedon=mesh.Boundaries('rub'))
#::::::::::::::::::::::::::::::::::::::

omega=2*pi*50 * 2
mu0 = 1.257e-6
mu = 1000 * mu0 *5
#rel = 1/mu
rel = mesh.MaterialCF({ 'air' : 1/mu0 }, default = 1/mu)

sigma = 2e6
#rho= 1/sigma 
rho = mesh.MaterialCF({ 'air' : 100 }, default = 1/sigma)


#Bfe=100*x+0.5
Bfe = -0.287862 /fill_fact
B0 = mesh.MaterialCF({ 'air' : Bfe/mu*mu0 }, default=Bfe)

term1 = rho * grad(csp)*grad(theta)*dx + 1j * omega*curl(mvp)*theta*dx + 0.001*mvp*alpha*dx
term2 = - 1j*omega* rel*curl(mvp)*curl(alpha)*dx + 1j*omega* csp*curl(alpha)*dx 
#term2 = - csp*curl(alpha)*dx + rel*curl(mvp)*curl(alpha)*dx 

a = BilinearForm(term1+term2)
a.Assemble()

force = -1j * omega*B0*theta*dx
f=LinearForm(force)
f.Assemble()

r = f.vec - a.mat * sol.vec
sol.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs()) * r

#PPPPPPPPPPPPPPPPPPPP

rot=CF( (0 , 1,  -1, 0), dims=(2,2) )
J = rot*grad(Tpot)
B=curl(Apot)+B0

Pow=0.5*rho*J*Conj(J)
print('PowLoss =', Integrate(Pow,mesh,order=5).real)


core = mesh.Materials('lamele')
Peddy=Integrate(Pow, mesh, order=5)
print('Peddy',Peddy)

Draw(Apot, mesh, "A")
Draw (curl(Apot), mesh, "Btilda")
Draw (J, mesh, "J")
Draw (B, mesh, "B")


int_Km = Integrate(CF(1), mesh)


b_sqr=B*B
b_bar= Integrate(B, mesh)/int_Km
int_bsqr = Integrate(b_sqr, mesh)
FF_avg=(1/int_Km * int_bsqr/(b_bar**2))
print('FF=1/Km * int_bsqr/bbar^2= rel* =',FF_avg)

#j_sqr=(J*Conj(J)).Norm()
j_sqr=J*J
int_jsqr = Integrate(j_sqr, mesh)
GG_avg =(1/int_Km * int_jsqr/(b_bar**2))
print('GG=1/Km * int_jsqr/jbar^2= rho* =',GG_avg)

nuzz=1/mu*FF_avg - 1j/omega/sigma *GG_avg
print('1/mu*FF_avg=',1/mu*FF_avg)
print('- 1j/omega/sigma *GG_avg=',- 1j/omega/sigma *GG_avg)
print('nuzz',nuzz)


jabs_sqr=(J.Norm())**2
int_jabs_sqr = Integrate(jabs_sqr, mesh)
GGabs_avg =(1/int_Km * int_jabs_sqr/abs(b_bar)**2)

Pcoeff=0.5/sigma * GGabs_avg #(1/int_Km * abs(int_jsqr)/((abs(b_bar))**2))
print('Pcoeff =', Pcoeff)
print('P_Km = Pcoeff*b_bar^2 * int_Km=', Pcoeff*abs(b_bar)**2 * int_Km)

area= Integrate(CF(1), mesh)
print('area',area)
print('Bavg=',Integrate(B, mesh)/area)


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
