
from ngsolve import *
from netgen.occ import *

#air= Circle((0,0), 0.008).Face()

size= 0.03 #0.001003
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

bb=CF((0.121*y,-0.121*x))
#bb=CF((10,10))

#gfu.Set(bb, VOL_or_BND = BND)
gfu.Set(bb, definedon=mesh.Boundaries('rub'))

omega=314
mu0 = 1.257e-6


sigma = 2e6
sig = mesh.MaterialCF({ "lamele" : 1 }, default=None)
#sigma=CoefficientFunction( (10, 0,  0, 10), dims=(2,2) )
#rel=CoefficientFunction( ( 1200, 0, 0, 1200), dims=(2,2) )


rel_c = 830.040856545774+269.19114611986697*1j
rel_kut =1207.3430971214286+570.8626805860135*1j
rel_xbrid=1497.3752288165838+1476.8386120874934*1j
rel_ybrid=1504.7424295248843+1473.0015217785444*1j
#rel=987.0305289254093-790.0068241799272*1j
rel= IfPos((x-w/3)*(x+w/3), IfPos((y-h/3)*(y+h/3),rel_kut,rel_xbrid), \
           IfPos((y-h/3)*(y+h/3),rel_ybrid,rel_c))


Pcoeff_c = 40427.44309667297
Pcoeff_xbrid = 165002.40117130266
Pcoeff_kut = 92169.71814805517
Pcoeff_ybrid=165478.3309682206
#Pcoeff=124031.07139624856
Pcoeff= IfPos((x-w/3)*(x+w/3), IfPos((y-h/3)*(y+h/3),Pcoeff_kut,Pcoeff_xbrid), \
           IfPos((y-h/3)*(y+h/3),Pcoeff_ybrid,Pcoeff_c))



a = BilinearForm(fes)
a += (1/mu0)*curl(u)*curl(v)*dx('air') + rel*curl(u)*curl(v)*dx('lamele')+ 1j*1e-4*u*v*dx 

a.Assemble()

sila=CF((0,0))
f = LinearForm(fes)
f += sila*v*dx
f.Assemble()

#solvers.BVP(bf=a, lf=f, gf=gfu, pre=None, maxsteps=200, print=True)

r = - a.mat * gfu.vec
gfu.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs())*r

A=gfu
B = curl(gfu)

core = mesh.Materials('lamele')
Pow=Pcoeff*(B*Conj(B)).Norm()
Peddy= Integrate(Pow, mesh, definedon=core)
print('Peddy',Peddy)

tok= Integrate(B.Norm(), mesh, definedon=core)
print('tok',tok)

Draw(A, mesh, "A")
Draw (B, mesh, "B")
Draw (Pow, mesh, "Pow")
Draw(Pcoeff,mesh, "Pcoeff")


