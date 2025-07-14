import netgen.gui
import ngsolve
from netgen.read_gmsh import ReadGmsh
from netgen.meshing import *
from ngsolve import *

import matplotlib.pyplot as plt
import sys
sys.argv = ["fun"]

ngsglobals.msg_level = 5

#Naziv 8BCC oznacava OsminuBoxCoreCoil
# Učitavanje Gmsh mreže
mesh_g = ReadGmsh("gmsh_BoxCoilCore_full/box_coil_boxH22.msh")
#mesh_g = ReadGmsh("BoxCoreCoil_eighth1.msh")

mesh = ngsolve.Mesh(mesh_g)
Draw(mesh)

print(mesh.GetBoundaries())
#----------------------


#test=H1(mesh, dirichlet='gamaB')
test=H1(mesh, dirichlet='tblr')
#test=H1(mesh, dirichlet='topleft')
gft=GridFunction(test)
gft.Set(2,BND)
Draw(gft,mesh,'gft')

