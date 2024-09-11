
from ngsolve import *
from netgen.occ import *

#air= Circle((0,0), 0.008).Face()

size=0.0003
n=9
d=0.001
ins=0.0001
w=(n*d+(n-1)*ins)/2
h=(n*d+(n-1)*ins)/2

air = MoveTo(-1.2*w, -1.2*h).Rectangle(2.4*w,2.4*h).Face()
air.edges.name = 'rub'

lamele = MoveTo(-w,-h).Rectangle(2*w,2*h).Face()
lamele.faces.name="lamele"
lamele.faces.maxh=size
air -=lamele
air.faces.name="air"
air.faces.col = (1, 1, 0)

geo = Glue([air, lamele])


mesh = Mesh(OCCGeometry(geo, dim=2).GenerateMesh(maxh=0.003, quad_dominated=True))

fes = HCurl(mesh, order=0, dirichlet="rub",  complex=True, nograds = False) #CMPLX
gfu = GridFunction(fes)

print('mesh.nv',mesh.nv)
print('mesh.nedge', mesh.nedge)
print('mesh.ne', mesh.ne)

u, v = fes.TnT()

bb=CF((0.1*y,-0.1*x))
gfu.Set(bb, definedon=mesh.Boundaries('rub'))

omega=314
mu0 = 1.257e-6
mu = mu0 * 1000
rel = 1/mu
sigma = 2e6
sig = mesh.MaterialCF({ "lamele" : 1 }, default=None)
#sigma=CoefficientFunction( (10, 0,  0, 10), dims=(2,2) )
#rel=CoefficientFunction( ( 1200, 0, 0, 1200), dims=(2,2) )



a = BilinearForm(fes)
a += (1/mu0)*curl(u)*curl(v)*dx('air') + rel*curl(u)*curl(v)*dx('lamele') + \
    +1j*omega*sigma*u*v*dx('lamele') + 1j*1e-3*u*v*dx('air') 

a.Assemble()

sila=CF((0,0))
f = LinearForm(fes)
f += sila*v*dx
f.Assemble()

#solvers.BVP(bf=a, lf=f, gf=gfu, pre=None, maxsteps=200, print=True)

r = - a.mat * gfu.vec
gfu.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs())*r

A=gfu

E = - 1j * omega * A
J = - 1j * omega * sigma * A *sig
B = curl(gfu)

Pow=0.5*E*Conj(J) 
Peddy=Integrate(Pow, mesh, order=5)
print('Peddy',Peddy)

""" reg = mesh.Materials('lamele')
Pcoeff=7410
Pow=Pcoeff*(B*Conj(B)).Norm()
Peddy= Integrate(Pow, mesh, definedon=reg)
print('Peddy',Peddy) """

Draw(A, mesh, "A")
Draw (B, mesh, "B")
Draw (Pow, mesh, "Pow")



w=(n*d+(n-1)*ins)/2
f_Km = IfPos( (x-w/3)*(x+w/3), 0, IfPos( (y-h/3)*(y-h), 0, 1))
f_core = IfPos( (x-w)*(x+w), 0, IfPos( (y-h)*(y+h), 0, 1))
reg_Km =f_Km
#pozivom metode mesh.Materials() dobiva se objekt klase ngsolve.comp.Region
defon = mesh.Materials('lamele')
#int_Km = Integrate(CF(1)*reg_Km, mesh)
int_Km = Integrate(CF(1)*reg_Km, mesh, definedon=defon)

#b_abs = B.Norm()
#b_sqr=(B*Conj(B)).Norm()
b_abs=B
b_sqr=B*B
b_bar= Integrate(b_abs*reg_Km, mesh, definedon=defon)/int_Km
int_bsqr = Integrate(b_sqr*reg_Km, mesh, definedon=defon)
FF_avg=(1/int_Km * int_bsqr/(b_bar**2))
print('FF=1/Km * int_bsqr/bbar^2= rel* =',FF_avg)

#j_sqr=(J*Conj(J)).Norm()
j_sqr=J*J
int_jsqr = Integrate(j_sqr*reg_Km, mesh, definedon=defon)
GG_avg =(1/int_Km * int_jsqr/(b_bar**2))
print('GG=1/Km * int_jsqr/jbar^2= rho* =',GG_avg)

nuzz=rel*FF_avg - 1j/omega/sigma *GG_avg
print('nuzz',nuzz)


babs_bar= Integrate(B.Norm()*reg_Km, mesh, definedon=defon)/int_Km
jabs_sqr=(J.Norm())**2
int_jabs_sqr = Integrate(jabs_sqr*reg_Km, mesh, definedon=defon)
GGabs_avg =(1/int_Km * int_jabs_sqr/babs_bar**2)

Pcoeff=0.5/sigma * GGabs_avg #(1/int_Km * abs(int_jsqr)/((abs(b_bar))**2))
print('Pcoeff =', Pcoeff)
print('P_Km = Pcoeff*b_bar^2 * int_Km=', Pcoeff*babs_bar**2 * int_Km)

tok= Integrate(B.Norm()*f_core, mesh, definedon=defon)
print('tok',tok)

print('b_bar(Km)=',b_bar)


Draw (sig, mesh, "sig")
Draw (f_core, mesh, "f_core")
Draw (reg_Km, mesh, "reg_Km")

