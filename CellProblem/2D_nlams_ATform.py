
from ngsolve import *
from netgen.occ import *

size=0.00005
n=9
d=0.001
ins=0.0001
fill_fact=n*d/(n*d+(n-1)*ins)

w=(n*d+(n-1)*ins)/2
h=(n*d+(n-1)*ins)/2

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
w=(n*d+(n-1)*ins)/2

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
mesh = Mesh(OCCGeometry(geo, dim=2).GenerateMesh(maxh=size, quad_dominated=False))

print(mesh.GetBoundaries())

##########################
##########################
Apot_dirich= "l0|t0|t1|t2|t3|t4|t5|t6|t7|t8|r8|b0|b1|b2|b3|b4|b5|b6|b7|b8|rub"
Tpot_dirich= "l0|r0|l1|r1|l2|r2|l3|r3|l4|r4|l5|r5|l6|r6|l7|r7|l8|r8|rub"
Svi_dirich= "l0|r0|l1|r1|l2|r2|l3|r3|l4|r4|l5|r5|l6|r6|l7|r7|l8|r8|t0|t1|t2|t3|t4|t5|t6|t7|t8|b0|b1|b2|b3|b4|b5|b6|b7|b8|rub"

#fsU = HCurl(mesh, order=0, complex=True, nograds = False)
fsU = HCurl(mesh, order=0, dirichlet=Apot_dirich, complex=True, nograds = False)
fsV = H1(mesh, order=1, dirichlet=Svi_dirich, complex=True)
fes=fsU*fsV
mvp, csp = fes.TrialFunction()
alpha, theta = fes.TestFunction()

print('mesh.nv',mesh.nv)
print('mesh.nedge', mesh.nedge)
print('mesh.ne', mesh.ne)

sol = GridFunction(fes)
Apot, Tpot = sol.components

#bb=CF((0.1*y,-0.1*x))
#Apot.Set(bb, definedon=mesh.Boundaries(Apot_dirich))
#::::::::::::::::::::::::::::::::::::::

omega=314 
mu0 = 1.257e-6
mu = 1000 * mu0
#rel = 1/mu
rel = mesh.MaterialCF({ 'air' : 1/mu0 }, default = 1/mu)

sigma = 2e6
#rho= 1/sigma 
rho = mesh.MaterialCF({ 'air' : 100 }, default = 1/sigma)


term1 = rho * grad(csp)*grad(theta)*dx + 1j * omega*curl(mvp)*theta*dx + 0.1*mvp*alpha*dx
term2 = - 1j*omega* rel*curl(mvp)*curl(alpha)*dx + 1j*omega* csp*curl(alpha)*dx 
#term2 = - csp*curl(alpha)*dx + rel*curl(mvp)*curl(alpha)*dx 

a = BilinearForm(term1+term2)
a.Assemble()


Hfe = (50*x - 0.5/fill_fact)/mu
#Bfe = 0.5 /fill_fact
B0 = mesh.MaterialCF({ 'air' : Hfe*mu0 }, default=Hfe*mu)

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

#pppppppppppppppppppp

#f_Km = CF(1)
#f_Km = IfPos( (x+w)*(x-w), 0, IfPos((y+h)*(y+h/3), 0, 1))
f_Km = IfPos( (x-w/3)*(x+w/3), 0, IfPos((y-h/3)*(y+h/3), 0, 1))
Draw (f_Km, mesh, "f_Km")
int_Km = Integrate(f_Km, mesh)
print('int_Km',int_Km)


b_bar= Integrate(B*f_Km, mesh)/int_Km
print('b_bar=', b_bar)
print('abs(b_bar)=', abs(b_bar))

b_sqr=B*B
int_bsqr = Integrate(b_sqr*f_Km, mesh)
FF_avg=(1/int_Km * int_bsqr/(b_bar**2))
print('FF=1/Km * int_bsqr/bbar^2= rel* =',FF_avg)

#j_sqr=(J*Conj(J)).Norm()
j_sqr=J*J
int_jsqr = Integrate(j_sqr*f_Km, mesh)
GG_avg =(1/int_Km * int_jsqr/(b_bar**2))
print('GG=1/Km * int_jsqr/jbar^2= rho* =',GG_avg)

nuzz=1/mu*FF_avg - 1j/omega/sigma *GG_avg
print('nuzz',nuzz)

jabs_sqr=(J.Norm())**2
int_jabs_sqr = Integrate(jabs_sqr*f_Km, mesh)
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
