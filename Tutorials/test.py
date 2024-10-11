from ngsolve import *

mesh = Mesh(unit_square.GenerateMesh(maxh=0.35))
#fes = H1(mesh, order=2)

fesL2 = L2(mesh, order=0)
gfL2 = GridFunction(fesL2)
gfL2.vec[:]=0
gfL2.vec[17]=1
vec = gfL2 * CF((1,0))
#Draw(vec, mesh, 'vec')

for i in range(fesL2.ndof):
    print (i,":", fesL2.CouplingType(i))

print('mesh.nv',mesh.nv)
print('mesh.nedge', mesh.nedge)


feh = HCurl(mesh, order=0)
gfh = GridFunction(feh)
#gfh.Set(vec)
gfh.vec[:]=0
gfh.vec[31] = 1
gfh.vec[32] = -1
hvec = gfh
Draw(hvec, mesh, 'hvec')