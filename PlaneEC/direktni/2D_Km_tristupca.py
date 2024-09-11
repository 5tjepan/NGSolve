from ngsolve import *
from netgen.occ import *

#HG = True
HG = False

if HG==True: size = 0.031
else: size = 0.0031

left = MoveTo(-0.045, -0.1).Rectangle(0.03,0.2).Face()
left.faces.name="left"
left.edges.Min(X).name = "ll"
left.edges.Min(Y).name = "lb"
left.edges.Max(Y).name = "lt"
left.faces.maxh= size

central = MoveTo(-0.015, -0.1).Rectangle(0.03,0.2).Face()
central.faces.name="central"
central.edges.Min(Y).name = "cb"
central.edges.Max(Y).name = "ct"
central.faces.maxh= size

right = MoveTo(0.015, -0.1).Rectangle(0.03,0.2).Face()
right.faces.name="right"
right.edges.Min(Y).name = "rb"
right.edges.Max(Y).name = "rt"
right.edges.Max(X).name = "rr"
right.faces.maxh= size

""" Km = MoveTo(0.02, -0.02).Rectangle(0.03,0.04).Face()
Km.faces.name = "Km"
Km.edges.Max(X).name = "d"
Km.faces.maxh= size
left -= Km
 """

left.faces.col = (1, 1, 0)  #colour
right.faces.col = (0, 1, 1) 

geo = Glue([left, central, right])
#Draw(OCCGeometry(geo));

mesh = Mesh(OCCGeometry(geo, dim=2).GenerateMesh(maxh=0.1, quad_dominated=True))

fsU = HCurl(mesh, order=0, dirichlet="lb|ll|lt|ct|rt|rr|rb|cb", complex=True, nograds = False)
fsV = H1(mesh, order=1, dirichlet="ll|rr", complex=True)
fes=fsU*fsV
mvp, csp = fes.TrialFunction()
alpha, theta = fes.TestFunction()

print('mesh.nv',mesh.nv)
print('mesh.nedge', mesh.nedge)
print('mesh.ne', mesh.ne)

sol = GridFunction(fes)
Apot, Tpot = sol.components

rot=CF( (0 , 1,  -1, 0), dims=(2,2) )
omega=314
mu0 = 1.257e-6
B0=CF(0.1)

if HG == True:
    
    omjer = mesh.MaterialCF({ "central" : 1}, default= 1)
    rho = mesh.MaterialCF({ "central" : 5e-7}, default= (5e-7))
    rel = mesh.MaterialCF({ "central" : 32000 }, default=32000)
    jedan = mesh.MaterialCF({ "central" : 1 }, default=1)

    term1 = rho * grad(csp)*grad(theta)*dx - 1j * jedan* omega*curl(mvp)*rot*grad(theta)*dx
    term2 = rel*curl(mvp)*curl(alpha)*dx + jedan* rot*grad(csp)*curl(alpha)*dx 

    a = BilinearForm(term1+term2)

    force = -1j * omega* B0*omjer * theta*dx
    f=LinearForm(force)

else:
   
    rho = mesh.MaterialCF({ "central" : 5e-7}, default= (5e-7))
    rel = mesh.MaterialCF({ "central" : 32000 }, default=32000)
    
    term1 = rho * grad(csp)*grad(theta)*dx + 1j * omega*curl(mvp)*theta*dx
    term2 = rel*curl(mvp)*curl(alpha)*dx - csp*curl(alpha)*dx

    a = BilinearForm(term1+term2)

    force = -1j * omega*B0*theta*dx
    f=LinearForm(force)

a.Assemble()
f.Assemble()

r = f.vec - a.mat * sol.vec
sol.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs()) * r


J = rot*grad(Tpot)
B=curl(Apot)

Power=0.5*rho*J*Conj(J)
print('Power losses =', Integrate(Power,mesh,order=5).real)




#pozivom metode mesh.Materials() dobiva se objekt klase ngsolve.comp.Region
Km_reg = mesh.Materials('right')
int_Km = Integrate(CF(1), mesh, definedon =Km_reg)

j_abs = J.Norm()
j_sqr=(J*Conj(J)).Norm()
j_bar= Integrate(j_abs, mesh, definedon =Km_reg)/int_Km
int_jsqr = Integrate(j_sqr, mesh, definedon =Km_reg)
print('GG=1/Km * int_jsqr/jbar^2= rho* =',(1/int_Km * int_jsqr/(j_bar**2)))

b_abs = B.Norm()
b_sqr=(B*Conj(B)).Norm()
b_bar= Integrate(b_abs, mesh, definedon =Km_reg)/int_Km
int_bsqr = Integrate(b_sqr, mesh, definedon =Km_reg)
print('FF=1/Km * int_bsqr/bbar^2= rel* =',(1/int_Km * int_bsqr/(b_bar**2)))

""" t_abs = Tpot.Norm()
t_bar = Integrate(t_abs, mesh, definedon =Km_reg)/int_Km
int_tb=Integrate((Tpot*B).Norm(), mesh, definedon =Km_reg)
print('1/Km * int_tj/(jbar*tbar) = jedan* =',(1/int_Km * int_tb* 1/(b_bar*t_bar))) """


int_jb=Integrate((J*B).Norm(), mesh, definedon =Km_reg)
print('QF=1/Km * int_jb/(jbar*bbar) = int_jb =',(1/int_Km * int_jb* 1/(b_bar*j_bar)))

int_t=Integrate((Tpot).Norm(), mesh, definedon =Km_reg)
print('Q=1/Km * int_t/jbar =',(1/int_Km * int_t/j_bar))


Draw (B0, mesh, "B0")
Draw (j_bar, mesh, "j_bar")
Draw (j_sqr, mesh, "j_sqr")

#Draw(rel,mesh, "rel")

Draw(Tpot, mesh, "T")
Draw(Apot, mesh, "A")
Draw (B, mesh, "B")
Draw (J, mesh, "J")

