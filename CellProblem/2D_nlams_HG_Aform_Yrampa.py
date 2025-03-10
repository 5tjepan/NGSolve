
from ngsolve import *
from netgen.occ import *

#air= Circle((0,0), 0.008).Face()

size= 0.003 #0.003 #0.0026
n=9
d=0.001
ins=0.0001
fill_fact=n*d/(n*d+(n-1)*ins)
w=(n*d+(n-1)*ins)/2
b=w

air = MoveTo(-1.2*w, -1.2*b).Rectangle(2.4*w,2.4*b).Face()
air.edges.name = 'rub'

lamele = MoveTo(-w,-b).Rectangle(2*w,2*b).Face()
lamele.faces.name="lamele"
lamele.faces.maxh=size
air -=lamele
air.faces.name="air"
air.faces.col = (1, 1, 0)

geo = Glue([air, lamele])


mesh = Mesh(OCCGeometry(geo, dim=2).GenerateMesh(maxh=size, quad_dominated=True))

fes = HCurl(mesh, order=0, dirichlet="rub",  complex=True, nograds = False) #CMPLX
gfu = GridFunction(fes)
""" gfu.Load('my_sol.sol')
Draw(curl(gfu),mesh, 'gfu')
print(nvghch) """

print('mesh.nv',mesh.nv)
print('mesh.nedge', mesh.nedge)
print('mesh.ne', mesh.ne)

u, v = fes.TnT()

#bb=CF((0.1*y,-0.1*x))
B0=0.5
bb=CF((0*y,-B0*x/1.2**2))
gfu.Set(bb, definedon=mesh.Boundaries('rub'))

omega=314
mu0 = 1.257e-6
mu = 1000*mu0

#rel_c = 866.8376706501075+56.8317640205614*1j 
#Pcoeff_c = 8916.5
#rel_fringe = 866.8890138503145+51.334309604750686*1j 
#Pcoeff_fringe = 8053.7 

rel_c = 873.4430072918316+56.14529407577265*1j 
Pcoeff_c = 8932.6
rel_up = 879.971086619287+50.14283221132355*1j 
Pcoeff_up = 8475.6
rel_down = 870.6241716528501+50.64335679870008*1j 
Pcoeff_down = 7832.5

rel= IfPos((y-b/3)*(y+b/3), IfPos(y, rel_up, rel_down), rel_c)
Pcoeff= IfPos((y-b/3)*(y+b/3), IfPos(y, Pcoeff_up, Pcoeff_down), Pcoeff_c)

#rel= IfPos((y-b/3)*(y+b/3), rel_fringe, rel_c)
#Pcoeff= IfPos((y-b/3)*(y+b/3), Pcoeff_fringe, Pcoeff_c)

Draw(rel,mesh,'rel')
#rel=866.6663925842802 + 56.35083173383855*1j
#Pcoeff=8839.83

sigma = 2e6
#sig = mesh.MaterialCF({ "lamele" : 1 }, default=0)
sig = IfPos( (y+b)*(y-b), 0, 1)
Draw(sig, mesh, 'sig')

a = BilinearForm(fes)
a += 100*(1/mu0)*curl(u)*curl(v)*dx('air') + rel*curl(u)*curl(v)*dx('lamele')+ 1j*1e-4*u*v*dx 

a.Assemble()

""" Jg=1e5 *10
generator=CF((0,Jg))
f = LinearForm(fes)
f += generator*v*dx
f.Assemble() """

sila=CF((50/mu*sig,0))
f = LinearForm(fes)
f += sila*v*dx
f.Assemble()
Draw(sila,mesh,'sila')
#solvers.BVP(bf=a, lf=f, gf=gfu, pre=None, maxsteps=200, print=True)

r = f.vec - a.mat * gfu.vec
gfu.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs())*r

A=gfu
B = curl(gfu)

core = mesh.Materials('lamele')
Pow=Pcoeff*(B*Conj(B)).Norm()
Peddy= Integrate(Pow, mesh, definedon=core)
print('Peddy',Peddy)

area= Integrate(CF(1), mesh, definedon=core)
print('area= ',area)

tok= Integrate(B, mesh, definedon=core)
print('tok= ',tok)

print('Bavg= ',tok/area)


Draw(A, mesh, "A")
Draw (B, mesh, "B")
Draw (Pow, mesh, "Pow")


#gfu.Save('my_sol.txt')

