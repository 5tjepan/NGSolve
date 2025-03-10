
from ngsolve import *
import netgen.gui
from netgen.occ import *

#air= Circle((0,0), 0.008).Face()

size=0.00005
n=9
d=0.001
ins=0.0001
fill_fact=n*d/(n*d+(n-1)*ins)
w=(n*d+(n-1)*ins)/2
h=(n*d+(n-1)*ins)/2

air = MoveTo(-1.2*w, -1.2*h).Rectangle(2.4*w,2.4*h).Face()
air.edges.name = 'rub'
air.faces.name="air"
air.faces.col = (1, 1, 0)

lamele = MoveTo(-w,-h).Line(d,0.0).Line(0,2*h).Line(-d,0).Close().Face()
for i in range(n-1):
    w-=(d+ins)
    lamele += MoveTo(-w,-h).Line(d,0.0).Line(0,2*h).Line(-d,0).Close().Face()

w=(n*d+(n-1)*ins)/2
h=(n*d+(n-1)*ins)/2

lamele.faces.name="lamele"
lamele.faces.col = (0,0,1)
lamele.faces.maxh=size

air -= lamele

geo = Glue([air, lamele])

""" lams = {}
for i in range(n):
    lams[f"lam{i}"] = MoveTo(-a,-h).Line(d,0.0).Line(0,2*h).Line(-d,0).Close().Face()
    a-=(d+ins)
    lams[f"lam{i}"].faces.maxh=size
    air-= lams[f"lam{i}"]

lista = [lams[key] for key in lams]
lista.append(air)
geo = Glue(lista) """


mesh = Mesh(OCCGeometry(geo, dim=2).GenerateMesh(maxh=0.0004, quad_dominated=False))

fes = HCurl(mesh, order=0, dirichlet="rub",  complex=True, nograds = False) #CMPLX
gfu = GridFunction(fes)

print('mesh.nv',mesh.nv)
print('mesh.nedge', mesh.nedge)
print('mesh.ne', mesh.ne)

u, v = fes.TnT()


omega=100*pi *2
mu0 = 1.257e-6
mu = mu0 * 1000 *5
rel = 1/mu
#rel = mesh.MaterialCF({ 'lamele' : 1/mu }, default = 1/mu0)
sigma = 2e6
lam = mesh.MaterialCF({ 'lamele' : 1 }, default = 0)
sig = IfPos( (x+w)*(x-w), 0, 1)


a = BilinearForm(fes)
a += 1*(1/mu0)*curl(u)*curl(v)*dx('air') + rel*curl(u)*curl(v)*dx('lamele') + \
    +1j*omega*sigma*u*v*dx('lamele') + 1e0*u*v*dx('air') 

a.Assemble()

#(((((((((((((())))))))))))))
Bfe=0.5
B0 = mesh.MaterialCF({ 'air' : Bfe/mu*mu0 }, default=Bfe)
bb=CF((0*y,-B0*x/1.2**2))
gfu.Set(bb, definedon=mesh.Boundaries('rub'))


#Jg=CF((0,8000*(x-2*w)/mu*sig)) 
#Jg=CF((0,50*(400*(x-w))/mu*sig)) 
Jg=CF((0,-50/mu*sig)) # *0
#Jg=CF((0,-50/mu*sig)) # *0
#Jg=CF((0,-50*rel*sig)) # *0
f = LinearForm(fes)
f += Jg*v*dx
f.Assemble()

#(((((((((((((())))))))))))))
#solvers.BVP(bf=a, lf=f, gf=gfu, pre=None, maxsteps=200, print=True)

r = f.vec - a.mat * gfu.vec
gfu.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs())*r


#pppppppppppppppppppppppppppppppppp
A=gfu

E = - 1j * omega * A
J = - 1j * omega * sigma*lam * A #+ Jg
B = curl(gfu)

core = mesh.Materials('lamele')
f_Km = IfPos( (x-w/3)*(x-w), 0, IfPos( (y+h/3)*(y-h/3), 0, 1))
f_core = IfPos( (x-w)*(x+w), 0, IfPos( (y-h)*(y+h), 0, 1))
reg_Km =f_Km
Draw(f_core, mesh, 'f_core')

area= Integrate(CF(1)*f_core, mesh)
print('area= ',area)
tok= Integrate(B*f_core, mesh)
print('tok= ',tok)
print('Bavg= ',tok/area)


Draw(A, mesh, "A")
Draw (B, mesh, "B")
Draw (J, mesh, "J")



Pow=0.5*E*Conj(J)
Peddy=Integrate(Pow*f_core, mesh, order=5)
print('Peddy',Peddy)

Draw (Pow, mesh, "Pow")

int_Km = Integrate(CF(1)*reg_Km, mesh)

#b_abs = B.Norm()
#b_sqr=(B*Conj(B)).Norm()
b_sqr=B*B
b_bar= Integrate(B*reg_Km, mesh)/int_Km
int_bsqr = Integrate(b_sqr*reg_Km, mesh)
FF_avg=(1/int_Km * int_bsqr/(b_bar**2))
print('FF=1/Km * int_bsqr/bbar^2= rel* =',FF_avg)

#j_sqr=(J*Conj(J)).Norm()
j_sqr=J*J
int_jsqr = Integrate(j_sqr*reg_Km, mesh)
GG_avg =(1/int_Km * int_jsqr/(b_bar**2))
print('GG=1/Km * int_jsqr/jbar^2= rho* =',GG_avg)

nuzz=rel*FF_avg - 1j/omega/sigma *GG_avg
print('nuzz',nuzz)

jabs_sqr=(J.Norm())**2
int_jabs_sqr = Integrate(jabs_sqr*reg_Km, mesh)
GGabs_avg =(1/int_Km * int_jabs_sqr/abs(b_bar)**2)

Pcoeff=0.5/sigma * GGabs_avg #(1/int_Km * abs(int_jsqr)/((abs(b_bar))**2))
print('Pcoeff =', Pcoeff)
print('P_Km = Pcoeff*b_bar^2 * int_Km=', Pcoeff*abs(b_bar)**2 * int_Km)


#Draw (sig, mesh, "sig")
#Draw (f_core, mesh, "f_core")
#Draw (f_Km, mesh, "f_Km")


delta = sqrt(2/(omega*sigma*mu))

F_B = sigma*d*1j*omega/(8*(sinh((1+1j)*d/2/delta))**2) * (delta*(1-1j)/2*sinh((1+1j)*d/delta)+ d)
F_J = sigma*d*1j*omega/(8*(sinh((1+1j)*d/2/delta))**2) * (delta*(1-1j)/2*sinh((1+1j)*d/delta)- d)

F_R = F_B + F_J

print('F_B=',F_B)
print('F_J=',F_J)
print('F_R=',F_R)


fill_fact=n*d/(n*d+(n-1)*ins)
F_R_bulk= 1/(fill_fact/F_R + (1-fill_fact)*mu0)
print('F_R_bulk=',F_R_bulk)

F_pow =sigma*omega**2*delta*d/8*(sinh(d/delta)-sin(d/delta))/(cosh(d/delta)-cos(d/delta))/fill_fact
print('F_pow', F_pow)


