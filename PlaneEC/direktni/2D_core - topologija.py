from ngsolve import *
from netgen.occ import *

#outer= Circle((0,0), 2).Face()
#outer.edges.name = 'rub'
#outer = MoveTo(-2, -2).Rectangle(4,4).Face()
#outer.edges.name="rub"

""" outer = MoveTo(0, -1).Rectangle(3,2).Face()
outer.edges.name="rub"
outer.edges.Max(X).name = "r"
outer.edges.Min(X).name = "l"
outer.edges.Min(Y).name = "b"
outer.edges.Max(Y).name = "t"
#outer.edges.Min(X).maxh=0.1 """

#inner = MoveTo(-0.8,-0.8).Line(1.6,0).Line(0,1.6).Line(-1.6,0).Close().Face()
#inner = MoveTo(-0.8,-0.8).Line(1.6,0.3).Line(0,1.6).Line(-1.6,-0.3).Close().Face()
inner = MoveTo(-0.8,-0.8).Line(1.6,-0.1).Line(0,1.8).Line(-1.6,-0.3).Close().Face()
#inner = MoveTo(0.0,-1.0).Line(1,1).Line(-1,1).Line(-1,-1).Close().Face()
#inner = MoveTo(-1,-0.5).Line(1.5,-0.5).Line(0.5,1.5).Line(-1.5,0.5).Close().Face()
#inner = MoveTo(0.0,-0.9).Line(1.1,0.9).Line(-1.1,0.9).Line(-1.0,-0.9).Close().Face()
#inner = MoveTo(0.0,-0.9).Line(1.1,0.9).Line(-1.1,0.9).Line(-1.1,-0.9).Close().Face()

inner.edges.name="interface"
inner.faces.maxh=0.51
#outer = outer - inner

inner.faces.name="inner"
inner.faces.col = (1, 1, 0)  #colour
#outer.faces.name="outer"

geo = inner
#geo = Glue([inner, outer])
#Draw(OCCGeometry(geo));


mesh = Mesh(OCCGeometry(geo, dim=2).GenerateMesh(maxh=0.51, quad_dominated=True))
fes = HCurl(mesh, order=0) #CMPLX
#fes = H1(mesh, order=1, dirichlet="rub",  complex=False) #CMPLX
gfu = GridFunction(fes)
print('mesh.nv',mesh.nv)
print('mesh.nedge', mesh.nedge)
print('mesh.ne', mesh.ne)


bb=CF((1*y,-1*x))
#bb=CF((10,10))

#gfu.Set(bb, VOL_or_BND = BND)
#gfu.Set(bb, definedon=mesh.Boundaries('rub'))



testu = GridFunction(fes)
testu.vec[:] = 0
k=27
testu.vec[k]=2  #[16]=2
#testu.vec[5]=2       #[29]=1.32
print('k=',k,'testu.vec[k]=',testu.vec[k])
Draw(testu)

B=curl(testu)
Draw(B,mesh,"B")



""" import numpy as np
dirich=GridFunction(fes)
dirich.Set(bb, definedon=mesh.Boundaries('rub'))

rows,cols,vals = a.mat.COO()
import scipy.sparse as sp
AA = sp.csr_matrix((vals,(rows,cols)))

S=AA.todense() #/795544.94828958

maska=fes.FreeDofs()
print(maska)
for dof in range(fes.ndof):
    if (not maska[dof]):
        S[dof,:]=0
        S[dof,dof]=1

print(S)
print(dirich.vec)
print(gfu.vec)

#determ=np.linalg.det(S)
Seigenvalues=np.linalg.eigvals(S)
print('eigenvaules = ',Seigenvalues)
 """
