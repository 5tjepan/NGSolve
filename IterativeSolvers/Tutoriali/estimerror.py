from ngsolve import *
from netgen.occ import *

import matplotlib.pyplot as plt
import sys
sys.argv = ["fun"]

def MakeGeometryOCC():
    base = Rectangle(1, 0.6).Face()
    chip = MoveTo(0.5,0.15).Line(0.15,0.15).Line(-0.15,0.15).Line(-0.15,-0.15).Close().Face()
    top = MoveTo(0.2,0.6).Rectangle(0.6,0.2).Face()
    base=base-chip

    base.faces.name="base"
    chip.faces.name="chip"
    top.faces.name="top"
    chip.faces.col=(1,0,0)
    top.faces.col=(0.5,0,1)
    #geo = base + chip + top
    geo = Glue([base,chip,top])
    Draw(OCCGeometry(geo, dim=2)) #Ako je koristen OCC, onda
    # se to mora naglasiti kod koristenja Draw i Mesh metode

    geo.edges.name='default'
    geo.edges.Min(Y).name='bot'

    chip.faces.maxh=0.05
    return OCCGeometry(geo,dim=2)

mesh=Mesh(MakeGeometryOCC().GenerateMesh(maxh=0.1))
Draw(mesh)
print(mesh.GetBoundaries())
print(mesh.GetMaterials())
fes = H1(mesh, order=2, dirichlet=[2,6,1])
u, v = fes.TnT()

#OVO NIJE HEAT EQUATION NEGO POISSON EQUATION FOR HEAT
#Za razliku od POISSON_EQ, HEAT_EQ ima d/dt tj vremsku evoluciju,
#POISSONOVA TOPLINSKA JEDNADZBA PREDSTAVLJA STACIONARNO STANJE (limes u vremenu)
#tj. konacnu razdiobu temeperature kada vise ne mjerimo promjenu temperature
#za zadane fiksirane (vremenski nepormjenjive) vrijednosti temerature
#za HEAT_source i granicu/granice domene (Dirichlet)
""" a1=1*grad(u)*grad(v)*dx('base')
a2=1000*grad(u)*grad(v)*dx('chip')
a3=10*grad(u)*grad(v)*dx('top')
a=BilinearForm(a1+a2+a3) """
lam=CoefficientFunction([1, 1000, 10])
a=BilinearForm(lam*grad(u)*grad(v)*dx)

f=LinearForm(1*v*dx('chip'))

c=Preconditioner(a, type='multigrid', inverse='sparsecholesky')

sol = GridFunction(fes)

def SolveBVP():
    fes.Update()  #nakon adaptacije mesha treba updejtat fes
    sol.Update()  #...i treba updejtat GridFunkciju za updejtani fes
    a.Assemble()  #i onda ponovno assemblirat STIFF matricu ...
    f.Assemble()
    inv = CGSolver(a.mat, c.mat)
    sol.vec.data=inv*f.vec

SolveBVP()
Draw(sol)

###ERROR ESTIMATOR:
DIVfes = HDiv(mesh, order=1)
gf_flux = GridFunction(DIVfes, "flux") #grid funkcija iz Hdiv
flux = lam * grad(sol) #izracunato (po dijelovima konst.) rjesenje za toplinu (sol pripada H1)

gf_flux.Set(flux) #ovakvom interpolacijom dobivamo Hdiv toplinski vektor (nema skokova normane komponente) 
error = 1/lam*(flux-gf_flux)*(flux-gf_flux)
Draw(error, mesh, 'error_representation')

Draw(flux,mesh,'flux')
Draw(gf_flux, mesh,'gf_flux')

eta2 = Integrate(error, mesh, VOL, element_wise=True)
print(eta2)

maxerr = max(eta2)
print ("maxerr = ", maxerr)

for el in mesh.Elements():
    mesh.SetRefinementFlag(el, eta2[el.nr] > 0.25*maxerr)
    # see below for vectorized alternative

#mesh.Refine()
#SolveBVP()
#Draw(sol)

#...automatiziranje adaptacije mi zasad nije vazno
#... to je lekcija za neka druga vremena