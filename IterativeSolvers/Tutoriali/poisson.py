from ngsolve import *
msh=Mesh(unit_square.GenerateMesh(maxh=0.2))
Draw(msh)
print("broj elemenata = ", msh.ne)
print("broj bridova = ", msh.nedge)
print("broj vrhova = ", msh.nv)

h1fes=H1(msh,order=3, dirichlet="left|top")
print("broj dofsa = ",h1fes.ndof)

u=h1fes.TrialFunction() # simbolicki objekt za izradu formulacije
v=h1fes.TestFunction()  # simbolicki objekt za izradu formulacije

stiff=BilinearForm(h1fes)
#print(stiff.mat) #javlja da matrica nije spremna; ne znam je li memorija alocirana
stiff+= grad(u)*grad(v)*dx 
stiff.Assemble()  #matrice sa stvarnim podacima u memoriji

force = LinearForm(h1fes)
#print("force",force.vec) printa nul-vektor; dakle, memorija je alocirana
force += x*v*dx #samo simbolicki izraz, zasad bez utjecaja na podatke
#print("force",force.vec) #...jer i dalje printa nulvektor
force.Assemble() #kreira vektor sa novim podacima u memoriji

#alternativno
#u, v = h1fes.TnT()
#stiff = BilinearForm(grad(u)*grad(v)*dx).Assemble()
#force = LinearForm(x*v*dx).Assemble()
#------------

#
gf=GridFunction(h1fes) #gridfunkcija sa stvarnim vrijednostima u memoriji (sve 0)
gf.vec.data = stiff.mat.Inverse(freedofs=h1fes.FreeDofs()) * force.vec
Draw(gf)

#print(gf.vec)


##### za kondenzaciju (eliminaciju lokalnih/unutarnjih/bubble dofsa)
dof_types={}
for i in range(h1fes.ndof):
    ctype = h1fes.CouplingType(i)
    if ctype in dof_types.keys():
        dof_types[ctype] += 1
    else:
        dof_types[ctype] = 1