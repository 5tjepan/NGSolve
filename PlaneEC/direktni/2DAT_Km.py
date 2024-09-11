from ngsolve import *
from netgen.occ import *

size=0.00008
d=0.001
h=0.01
core = MoveTo(-d/2, -h/2).Rectangle(d,h).Face()
#core = MoveTo(-0.08,-0.08).Line(0.16,0.0).Line(0,0.16).Line(-0.16,0).Close().Face()
#core.edges.name="interface"
#core.faces.maxh=0.0001
core.edges.Max(X).name = "r"
core.edges.Min(X).name = "l"
core.edges.Min(Y).name = "b"
core.edges.Max(Y).name = "t"
#air.edges.Min(X).maxh=0.1

Km = MoveTo(0.0, 0.0).Rectangle(d/2,h/10).Face()
Km.faces.name = "Km"
Km.edges.Max(X).name = "kmr"
Km.edges.Min(X).name = "kml"
Km.edges.Min(Y).name = "kmb"
Km.edges.Max(Y).name = "kmt"
Km.faces.maxh=size
core -= Km

core.faces.name="core"
core.faces.col = (1, 1, 0)  #colour

geo = Glue([Km, core])
#Draw(OCCGeometry(geo));

mesh = Mesh(OCCGeometry(geo, dim=2).GenerateMesh(maxh=size, quad_dominated=False))

fsU = HCurl(mesh, order=0, dirichlet="b|l|t|r|kmr", complex=True, nograds = False)
fsV = H1(mesh, order=1, dirichlet="l|r|kmr", complex=True)
#fsV = H1(mesh, order=1, dirichlet="interface", definedon="core", complex=True)
fes=fsU*fsV
mvp, csp = fes.TrialFunction()
alpha, theta = fes.TestFunction()

print('mesh.nv',mesh.nv)
print('mesh.nedge', mesh.nedge)
print('mesh.ne', mesh.ne)

sol = GridFunction(fes)
Apot, Tpot = sol.components

sig = mesh.MaterialCF({ "core|Km" : 1 }, default=0.000000001)
B0=sig*0.1

omega=314
mu0 = 1.257e-6
#rho = 0.00001

#rel = mesh.MaterialCF({ "core|Km" : 32000 }, default=(1/mu0))
#rho = mesh.MaterialCF({ "core|Km" : 6e-7 }, default=10)
rel = 32000
rho= 5e-7


#sigma=CoefficientFunction( (0, 0,  0, 10), dims=(2,2) )
#rel=CoefficientFunction( ( 1200, 0, 0, 1200), dims=(2,2) )

""" term1 = rho * grad(csp)*grad(theta)*dx('core') + 1j*omega*curl(mvp)*theta*dx('core')
term2 = rel*curl(mvp)*curl(alpha)*dx - csp*curl(alpha)*dx('core')
term3 = 1*mvp*alpha*dx('air') """


term1 = rho * grad(csp)*grad(theta)*dx + 1j * omega*curl(mvp)*theta*dx
term2 = - csp*curl(alpha)*dx + rel*curl(mvp)*curl(alpha)*dx 

a = BilinearForm(term1+term2)
a.Assemble()

force = -1j * omega*B0*theta*dx
f=LinearForm(force)
f.Assemble()

r = f.vec - a.mat * sol.vec
sol.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs()) * r

rot=CF( (0 , 1,  -1, 0), dims=(2,2) )
J = rot*grad(Tpot)
B=curl(Apot)

Power=0.5*rho*J*Conj(J)
print('Power losses =', Integrate(Power,mesh,order=5).real)




#pozivom metode mesh.Materials() dobiva se objekt klase ngsolve.comp.Region
Km_reg = mesh.Materials('Km')
int_Km = Integrate(CF(1), mesh, definedon =Km_reg)

#j_bar = Integrate(B0.Norm(), mesh=Km_reg)
#j_sqr=Integrate((B0*B0).Norm(), mesh=Km_reg)

j_abs = J.Norm()
j_sqr=(J*Conj(J)).Norm()
j_bar= Integrate(j_abs, mesh, definedon =Km_reg)/int_Km
int_jsqr = Integrate(j_sqr, mesh, definedon =Km_reg)
print('1/Km * int_jsqr/jbar^2 =',(1/int_Km * int_jsqr/(j_bar**2)))

t_abs = Tpot.Norm()
t_bar = Integrate(t_abs, mesh, definedon =Km_reg)/int_Km
int_tj=Integrate((Tpot*J).Norm(), mesh, definedon =Km_reg)
print('1/Km * int_tj/(jbar*tbar) =',(1/int_Km * int_tj* 1/(j_bar*t_bar)))

b_abs = B.Norm()
b_sqr=(B*Conj(B)).Norm()
b_bar= Integrate(b_abs, mesh, definedon =Km_reg)/int_Km
int_bsqr = Integrate(b_sqr, mesh, definedon =Km_reg)
print('1/Km * int_bsqr/bbar^2 =',(1/int_Km * int_bsqr/(b_bar**2)))




Draw (B0, mesh, "B0")
Draw (j_bar, mesh, "j_bar")
Draw (j_sqr, mesh, "j_sqr")

#Draw(rel,mesh, "rel")

Draw(Tpot, mesh, "T")
Draw(Apot, mesh, "A")
Draw (B, mesh, "B")
Draw (J, mesh, "J")

