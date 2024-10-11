from ngsolve import *

mesh = Mesh(unit_square.GenerateMesh(maxh=0.35))
fes = H1(mesh, order=2)
gfu1 = GridFunction(fes)
print('mesh.nv',mesh.nv)
print('mesh.nedge', mesh.nedge)

gfu1.vec[:]=0
gfu1.vec[16] = 1
Draw(gfu1, autoscale=False)

gfu = GridFunction(fes)
nodedofs=fes.GetDofNrs(NodeId(EDGE,32))
print('nodedofs tuple', nodedofs)
gfu.vec[:]=0
#gfu.vec[nodedofs[1]]=1
gfu.vec[nodedofs[0]]=-1
Draw(gfu)

for i in range(fes.ndof):
    print (i,":", fes.CouplingType(i))

fecurl=HCurl(mesh, order=0)
uc= GridFunction(fecurl, name='uc')

edgedofs = fecurl.GetDofNrs(NodeId(EDGE,17))
print('edgedofs =',edgedofs)
uc.vec[:]=0
#uc.vec[edgedofs[0]] = 1
#uc.vec[edgedofs[1]] = 1
uc.vec[17]=1
uc.vec[30]=0.2
Draw(uc)
rot=curl(uc)
Draw(rot, mesh,'rot')

vector=CF((2,0))
fes2 = HCurl(mesh, order=0)
gf2 = GridFunction(fes2)
gf2.Set(vector)
vec = gf2
Draw(vec, mesh, 'vec')




