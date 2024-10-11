from ngsolve import *
import netgen.gui
from netgen.geom2d import unit_square
import tkinter as tk

def solve(value):
     global mesh
     fes = H1(mesh, order=2, dirichlet="bottom|right|left")
     u, v, gfu = *fes.TnT(), GridFunction(fes)
     a = BilinearForm(grad(u)*grad(v)*dx).Assemble()
     f = LinearForm((exp(-((x-0.5)**2+(y-0.5)**2)/float(value)))*v*dx).Assemble()
     gfu.vec.data = a.mat.Inverse(freedofs=fes.FreeDofs()) * f.vec
     Draw(gfu)

def run():
     mesh =  Mesh(unit_square.GenerateMesh(maxh=0.05))
     globals().update({"mesh":mesh})
     solve(1e-3)
     root = netgen.gui.win
     current_value = tk.DoubleVar()
     slider = tk.Scale(root,from_=1e-3,to=1e-1,resolution=1e-3,
orient="horizontal", variable=current_value, command=solve)
     slider.pack()
     input("stop")
     slider.destroy()

if __name__ == "__main__":
     run()