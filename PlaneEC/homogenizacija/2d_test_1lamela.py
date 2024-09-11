

from ngsolve import *
from netgen.occ import *


size=0.00002
n=1
d=0.001
ins=0.0001
fill_fact=n*d/(n*d+(n-1)*ins)
h=3*d


#!!!!!!!!!!!!!!::::

lam = MoveTo(-d/2,-h/2).Rectangle(d,h).Face()
lam.faces.name='lam'
lam.edges.Min(X).name = 'l'
lam.edges.Max(X).name = 'r'
lam.edges.Min(Y).name = 'b'
lam.edges.Max(Y).name = 't'
lam.faces.col = (0.2,0.3,1)

geo = lam


mesh = Mesh(OCCGeometry(geo, dim=2).GenerateMesh(maxh=size, quad_dominated=False))

print(mesh.GetBoundaries())

##########################
##########################


#fsU = HCurl(mesh, order=0, complex=True, nograds = False)
fsU = HCurl(mesh, order=0, dirichlet='l|t|r|b', complex=True, nograds = False)
fsV = H1(mesh, order=1, dirichlet='l|r', complex=True)
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

omega=314 *2
mu0 = 1.257e-6
mu = 1000 * mu0 *5
#rel = 1/mu
rel = mesh.MaterialCF({ 'air' : 1/mu0 }, default = 1/mu)

sigma = 2e6
#rho= 1/sigma 
rho = mesh.MaterialCF({ 'air' : 100 }, default = 1/sigma)

lam_ins = mesh.MaterialCF({ 'air' : 0 }, default = 1)

term1 = rho * grad(csp)*grad(theta)*dx + 1j * omega*curl(mvp)*theta*dx + 0.1*mvp*alpha*dx
term2 = - 1j*omega* rel*curl(mvp)*curl(alpha)*dx + 1j*omega* csp*curl(alpha)*dx 
#term2 = - csp*curl(alpha)*dx + rel*curl(mvp)*curl(alpha)*dx 

a = BilinearForm(term1+term2)
a.Assemble()


Hfe = (50*y*0 - 0.5/fill_fact)/mu
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

Pow=0.5*rho*J*Conj(J) #*lam_ins
print('PowLoss =', Integrate(Pow,mesh,order=5).real)


core = mesh.Materials('lamele')
Peddy=Integrate(Pow, mesh, order=5)
print('Peddy',Peddy)

Draw(Apot, mesh, "A")
Draw (curl(Apot), mesh, "Btilda")
Draw (J, mesh, "J")
Draw (B, mesh, "B")

#pppppppppppppppppppp

f_Km = CF(1)
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

#..............
jabs_sqr=(J.Norm())**2
int_jabs_sqr = Integrate(jabs_sqr*f_Km, mesh)
GGabs_avg =(1/int_Km * int_jabs_sqr/(abs(b_bar))**2)

Pcoeff=0.5/sigma * GGabs_avg #(1/int_Km * abs(int_jsqr)/((abs(b_bar))**2))
print('Pcoeff =', Pcoeff)
print('P_Km = Pcoeff*b_bar^2 * int_Km=', Pcoeff*abs(b_bar)**2 * int_Km)

area= Integrate(CF(1), mesh)
print('area',area)
print('Bavg=',Integrate(B, mesh)/area)

#######

#Draw(lam_ins,mesh,'lam_ins')
fill_fact= Integrate(lam_ins, mesh)/area

delta = sqrt(2/(omega*sigma*mu))
print('delta = ', delta)

F_B = sigma*d*1j*omega/(8*(sinh((1+1j)*d/2/delta))**2) * (delta*(1-1j)/2*sinh((1+1j)*d/delta)+ d)
F_J = sigma*d*1j*omega/(8*(sinh((1+1j)*d/2/delta))**2) * (delta*(1-1j)/2*sinh((1+1j)*d/delta)- d)

F_R = F_B + F_J

print('F_B=',F_B)
print('F_J=',F_J)
print('F_R=',F_R)


F_R_bulk= 1/(fill_fact/F_R + (1-fill_fact)*mu0)
print('F_R_bulk=',F_R_bulk)

""" gamma=(1+1j)/delta
loss1D = 0.5*rho* (b_bar * (-1j*omega*sigma*d/2) * sinh(gamma) ) """

Pow_analitik =sigma*omega**2*delta*d/8*(sinh(d/delta)-sin(d/delta))/(cosh(d/delta)-cos(d/delta)) * (abs(b_bar))**2
Watti=Integrate(Pow_analitik * f_Km*lam_ins, mesh)
print('Watti',Watti)

F_pow =sigma*omega**2*delta*d/8*(sinh(d/delta)-sin(d/delta))/(cosh(d/delta)-cos(d/delta))
print('F_pow', F_pow)


print('')

