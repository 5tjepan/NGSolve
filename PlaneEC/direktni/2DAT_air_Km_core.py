from ngsolve import *
from netgen.occ import *

air= Circle((0,0), 0.2).Face()
air.edges.name = "rub"

core = MoveTo(-0.1, -0.1).Rectangle(0.2,0.2).Face()
#core = MoveTo(-0.08,-0.08).Line(0.16,0.0).Line(0,0.16).Line(-0.16,0).Close().Face()
#core.edges.name="interface"
core.faces.maxh=0.01
core.edges.Max(X).name = "r"
core.edges.Min(X).name = "l"
core.edges.Min(Y).name = "b"
core.edges.Max(Y).name = "t"
#air.edges.Min(X).maxh=0.1
air = air - core
Km = MoveTo(0.06, -0.02).Rectangle(0.04,0.04).Face()
Km.faces.name = "Km"
Km.edges.Max(X).name = "d"
core -= Km

core.faces.name="core"
core.faces.col = (1, 1, 0)  #colour
air.faces.name="air"

geo = Glue([Km, core, air])
#Draw(OCCGeometry(geo));

mesh = Mesh(OCCGeometry(geo, dim=2).GenerateMesh(maxh=0.1, quad_dominated=False))

fsU = HCurl(mesh, order=0, dirichlet="rub", complex=True, nograds = False)
fsV = H1(mesh, order=1, dirichlet="b|l|t|r|d", complex=True)
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
B0=sig*CF(0.1)

omega=1
mu0 = 1.257e-6
#rho = 0.00001

rel = mesh.MaterialCF({ "core|Km" : 32000 }, default=(1/mu0))
rho = mesh.MaterialCF({ "core|Km" : 5e-7 }, default=10)
#sigma = 2e3
#rel=30000

#sigma=CoefficientFunction( (0, 0,  0, 10), dims=(2,2) )
#rel=CoefficientFunction( ( 1200, 0, 0, 1200), dims=(2,2) )

""" term1 = rho * grad(csp)*grad(theta)*dx('core') + 1j*omega*curl(mvp)*theta*dx('core')
term2 = rel*curl(mvp)*curl(alpha)*dx - csp*curl(alpha)*dx('core')
term3 = 1*mvp*alpha*dx('air') """


term1 = rho * grad(csp)*grad(theta)*dx + 1j * omega*curl(mvp)*theta*dx
term2 = - csp*curl(alpha)*dx + rel*curl(mvp)*curl(alpha)*dx 
print('omega',omega)

a = BilinearForm(term1+term2)
a.Assemble()

""" nova = BilinearForm(term2+term3)
nova.Assemble()

razlika=[]
for i in range(fes.ndof):
    for j in range(fes.ndof):
        #alista.append(a.mat[i,j])
        #novalista.append(nova.mat[i,j])
        if a.mat[i,j] != nova.mat[i,j]: razlika.append((i,j,a.mat[i,j],nova.mat[i,j]))
print('razlika',razlika) """

force = -1j * omega*B0*theta*dx
f=LinearForm(force)
f.Assemble()

r = f.vec - a.mat * sol.vec
sol.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs()) * r

rot=CF( (0 , 1,  -1, 0), dims=(2,2) )
J = rot*grad(Tpot)
B=curl(Apot)+B0

Power=0.5*rho*J*Conj(J)
print(Integrate(Power,mesh,order=5))

#Draw (B0, mesh, "B0")
#Draw(rel,mesh, "rel")

Draw(Tpot, mesh, "T")
Draw(Apot, mesh, "A")
Draw (B, mesh, "B")
Draw (J, mesh, "J")

""" maskaU=fsU.FreeDofs()
print('maskaU',maskaU)
maskaV=fsV.FreeDofs()
print('maskaV',maskaV) """
