
from ngsolve import *
from netgen.occ import *

#air= Circle((0,0), 0.008).Face()

size= 0.003 #0.001055 #0.00155 #0.00115 #0.003 #0.0026
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


omega=314
mu0 = 1.257e-6
mu = 1000*mu0*5


# OVO JE ZA f=500Hz sa ukljucenim Axn=0 među lamelama :::::::
""" rel_lc= 379.1588984758098+396.1929106123431*1j
Pcoeff_lc = 622444.7645075723
rel_lf= 357.2442019492663+341.60832766314024*1j
Pcoeff_lf = 536720.7197777092
rel_left = IfPos((y-b/3)*(y+b/3), rel_lf, rel_lc)
Pcoeff_left = IfPos((y-b/3)*(y+b/3), Pcoeff_lf, Pcoeff_lc)
#
rel_cc= 381.24098936307723+398.270994986128*1j
Pcoeff_cc = 625725.4562928587
rel_cf= 359.1896829767446+343.4120632595603*1j
Pcoeff_cf = 539580.6228634872
rel_centar = IfPos((y-b/3)*(y+b/3), rel_cf, rel_cc)
Pcoeff_centar = IfPos((y-b/3)*(y+b/3), Pcoeff_cf, Pcoeff_cc)
#
rel_rc= 383.37760232537755+400.3926236893971*1j
Pcoeff_rc = 628995.7658171107
rel_rf=361.23126028111864+345.2188423021463*1j
Pcoeff_rf = 542380.8425109958
rel_right = IfPos((y-b/3)*(y+b/3), rel_rf, rel_rc)
Pcoeff_right = IfPos((y-b/3)*(y+b/3), Pcoeff_rf, Pcoeff_rc) """
##..
#...

# OVO je za f=50Hz za Jg=CF((0,8000*(x-2*w)/mu*sig))
""" rel_lc= 179.21592422853905+55.89117073666316j
Pcoeff_lc = 8986.263556672902
rel_lf= 179.37819294000954+50.333755667288216*1j
Pcoeff_lf = 8096.075374714616
rel_left = IfPos((y-b/3)*(y+b/3), rel_lf, rel_lc)
Pcoeff_left = IfPos((y-b/3)*(y+b/3), Pcoeff_lf, Pcoeff_lc)
#
rel_cc= 180.14455520090058+55.61863996987342*1j
Pcoeff_cc = 9027.271129564964
rel_cf= 180.32094917281026+50.0841043559832*1j
Pcoeff_cf = 8133.743360181129
rel_centar = IfPos((y-b/3)*(y+b/3), rel_cf, rel_cc)
Pcoeff_centar = IfPos((y-b/3)*(y+b/3), Pcoeff_cf, Pcoeff_cc)
#
rel_rc= 181.09620028779895+55.340159421875654*1j
Pcoeff_rc = 9067.656229935186
rel_rf= 181.28746966059902+49.82867505438707*1j
Pcoeff_rf = 8170.749793098557
rel_right = IfPos((y-b/3)*(y+b/3), rel_rf, rel_rc)
Pcoeff_right = IfPos((y-b/3)*(y+b/3), Pcoeff_rf, Pcoeff_rc) """
##..
#...

# OVO JE ZA f=500Hz bez Axn=0 među lamelama :::::::
""" rel_lc= 374.26193596005743+390.86685033289854*1j
Pcoeff_lc = 614603.1
rel_lf= 352.6496870656935+336.9939395749558*1j
Pcoeff_lf = 530000.0
rel_left = IfPos((y-b/3)*(y+b/3), rel_lf, rel_lc)
Pcoeff_left = IfPos((y-b/3)*(y+b/3), Pcoeff_lf, Pcoeff_lc)
#
rel_cc= 374.39750693377187+390.7380798438004*1j
Pcoeff_cc = 614807.8839675991
rel_cf= 352.78926640626923+336.86699042988545*1j
Pcoeff_cf = 530230.9337417324
rel_centar = IfPos((y-b/3)*(y+b/3), rel_cf, rel_cc)
Pcoeff_centar = IfPos((y-b/3)*(y+b/3), Pcoeff_cf, Pcoeff_cc)
#
rel_rc= 374.7650775301808+390.4051017822036*1j
Pcoeff_rc = 615334.4407072379
rel_rf=353.2181309440552+336.5060348780718*1j
Pcoeff_rf = 530763.0720819999
rel_right = IfPos((y-b/3)*(y+b/3), rel_rf, rel_rc)
Pcoeff_right = IfPos((y-b/3)*(y+b/3), Pcoeff_rf, Pcoeff_rc)"""
#...
#... 

# OVO JE ZA f=50Hz za  ::::: Jg=CF((0,-50/mu*sig))
""" rel_lc= 177.62957994382518+56.35618858244098*1j
Pcoeff_lc = 8916.112363575401
rel_lf= 177.76735775201735+50.75991278172035*1j
Pcoeff_lf = 8031.724008024225
rel_left = IfPos((y-b/3)*(y+b/3), rel_lf, rel_lc)
Pcoeff_left = IfPos((y-b/3)*(y+b/3), Pcoeff_lf, Pcoeff_lc)
#
rel_cc= 178.10944429359006+56.215082816127854*1j
Pcoeff_cc = 8938.192452639347
rel_cf= 178.2544407155758+50.63070431444271*1j
Pcoeff_cf = 8051.996369428854
rel_centar = IfPos((y-b/3)*(y+b/3), rel_cf, rel_cc)
Pcoeff_centar = IfPos((y-b/3)*(y+b/3), Pcoeff_cf, Pcoeff_cc)
#
rel_rc= 179.3458308020795+55.852195376725945*1j
Pcoeff_rc = 8995.145373027179
rel_rf= 179.5101811806547+50.29792507804661*1j
Pcoeff_rf = 8104.130491338208
rel_right = IfPos((y-b/3)*(y+b/3), rel_rf, rel_rc)
Pcoeff_right = IfPos((y-b/3)*(y+b/3), Pcoeff_rf, Pcoeff_rc) """
#...
#... 


# OVO JE ZA f=500Hz za  ::::: Jg=CF((0,-50/mu*sig))
rel_lc= 374.26193596005743+390.86685033289854*1j
Pcoeff_lc = 614603.1435603083
rel_lf= 352.6496870656935+336.9939395749558*1j
Pcoeff_lf = 530000.0451694002
rel_left = IfPos((y-b/3)*(y+b/3), rel_lf, rel_lc)
Pcoeff_left = IfPos((y-b/3)*(y+b/3), Pcoeff_lf, Pcoeff_lc)
#
rel_cc= 374.39750693377187+390.7380798438004*1j
Pcoeff_cc = 614807.8839675991
rel_cf= 352.78926640626923+336.86699042988545*1j
Pcoeff_cf = 530230.9337417324
rel_centar = IfPos((y-b/3)*(y+b/3), rel_cf, rel_cc)
Pcoeff_centar = IfPos((y-b/3)*(y+b/3), Pcoeff_cf, Pcoeff_cc)
#
rel_rc= 374.7650775301808+390.4051017822036*1j
Pcoeff_rc = 615334.4407072379
rel_rf= 353.2181309440552+336.5060348780718*1j
Pcoeff_rf = 530763.0720819999
rel_right = IfPos((y-b/3)*(y+b/3), rel_rf, rel_rc)
Pcoeff_right = IfPos((y-b/3)*(y+b/3), Pcoeff_rf, Pcoeff_rc)
#...
#... 

# OVO JE SVIMA ZAJEDNICKO:
rel= IfPos((x-b/3)*(x+b/3), IfPos(x, rel_right, rel_left), rel_centar)
Pcoeff= IfPos((x-b/3)*(x+b/3), IfPos(x, Pcoeff_right, Pcoeff_left), Pcoeff_centar) 


Draw(rel,mesh,'rel')

sigma = 2e6
sig = IfPos( (x+b)*(x-b), 0, 1)
Draw(sig, mesh, 'sig')

a = BilinearForm(fes)
a += 100*(1/mu0)*curl(u)*curl(v)*dx('air') + rel*curl(u)*curl(v)*dx('lamele')+ 1j*1e-4*u*v*dx 
a.Assemble()


#bb=CF((0.1*y,-0.1*x))
B0=0.5
bb=CF((0*y,B0*x/1.2**2))
gfu.Set(bb, definedon=mesh.Boundaries('rub'))

#Jg=CF((0,8000*(x-2*w)/mu*sig)) 
#Jg=CF((0,50*(200*(x-w))/mu*sig)) 
Jg=CF((0,50/mu*sig))
f = LinearForm(fes)
f += Jg*v*dx
f.Assemble()
Draw(Jg,mesh,'Jg')

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

""" import numpy as np
import matplotlib.pyplot as plt
import sys
sys.argv = ["dummy"]

X = np.linspace(-0.5*d, 0.5*d, num=21)
Y = np.zeros_like(X)+0.00065

plt.plot(X, (B).imag(mesh(X, Y)))
plt.xlabel('x')
plt.show() """




postfes = H1(mesh, order=1, definedon='lamele', complex=True) #dirichlet=Svi_dirich
Blin = postfes.TrialFunction()
beta = postfes.TestFunction()


posta = BilinearForm( ( Blin*beta)*dx).Assemble()
postf = LinearForm(B*beta*dx).Assemble()

Bset = GridFunction(postfes)
Bpost = GridFunction(postfes)
Bset.Set(B)
Bpost.vec.data = posta.mat.Inverse(postfes.FreeDofs()) * postf.vec

#c = Preconditioner(posta, type="direct", inverse="umfpack")
#solvers.BVP(posta,postf, gfpost, c)

Draw(Bpost, mesh, 'Bpost')
Draw(Bset, mesh, 'Bset')

#print('integral = ', Integrate(Bpost, mesh))

