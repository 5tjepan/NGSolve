from ngsolve import *

mesh = Mesh(unit_cube.GenerateMesh(maxh=1))
Draw(mesh)

""" for v in mesh.vertices:  #mesh.vertices je ?iterator? [i instanca klase MeshNodeRange?] unutar modula ngsolve.comp
    print('v =',v, 'v.point=',v.point)
    print(type(v)) #class 'ngsolve.comp.MeshNode'

#ako je mesh.vertices bio iterator, onda mu se nakon citanja vise ne moze pristupiti
#indeks je stigao do kraja iteratora
#sljedeca print naredba printa rezultat
print('mesh.vertices[3].point =',mesh.vertices[3].point)
#ali razlog za to je sto je mesh.vertices dohvatljiv i kao argument objekta mesh
#mada mi nije jasno kako to da je ?indeksabilan?

for el in mesh.Elements(VOL):  #Elements() je metoda unutar klase ngsolve.comp.Mesh (tj. klasa Mesh u modulu ngsolve.comp) koja vraca iterator koji je instanca klase ngsolve.comp.ElementRange 
    print('type(el)',type(el)) #class 'ngsolve.comp.Ngs_Element'
    print('vertices', el.vertices)
    print('edges', el.edges) """


vvv= NodeId(VERTEX,3) # v je dakle instanca klase NodeId u modulu ngsolve.comp
print(type(vvv))
print ("tip Nodea = ", vvv.type, "vvv.nr =", vvv.nr)

""" meshv=mesh[v]  #meshv je instanca klase MeshNode() u modulu ngsolve.comp koja je dijete klase NodeId
print('type',type(meshv))
print('point', meshv.point)
print(meshv.edges)
print(meshv.elements)
print(meshv.faces) """

c = NodeId(CELL, 11)  #ide do 11 jer ima 12 elemenata [0,...,11]
meshc = mesh[c]  #meshc je klase MeshNode()
print('tip za cell',type(meshc))  #...printa se class 'ngsolve.comp.MeshNode'
print(type(meshc.faces))  #meshc.faces je obicna instanca klase tuple
#print('point', meshc.point) #samo za VERTEX node
#print(meshc.edges)  #samo za VERTEX node
#print(meshc.elements)  #samo za VERTEX node
print(meshc.faces)

ed = NodeId(FACE, 1)
meshf = mesh[ed]
print('type',type(meshf))
print(meshf.edges)  #dostupno,...meshf.edges je opet NodeId objekt
#print(meshf.elements)  #samo za VERTEX node
#print(meshf.faces) #nisu dostupne za FACE node
#print('point', meshc.point) #samo za VERTEX node


ei = ElementId(BND,0)
print(type(ei)) # printa se <class 'ngsolve.comp.ElementId'>

meshel = mesh[ei]
type(meshel)
print ("meshel = \n ", meshel) #printa se <ngsolve.comp.Ngs_Element>
print ("vertices =", meshel.vertices)

brid = ElementId(BBND,3) #BBND je druga ko-dimenzija - a u 3D modelu je to onda 1D objekt
print(type(brid))  #printa se ngsolve.comp.ElementId
#ako zelim dobit topologiju ovog brida onda
print(mesh[brid].vertices)

print('...DOFS...:')
fes = H1(mesh, order=2) #instanca klase H1 u modulu ngsolve.comp
print(mesh.edges)
for edge in mesh.edges:
    print ("type = ", type(edge))
    print ("dofs = ", fes.GetDofNrs(edge))
    # metoda fes.GetDofNrs(ElementId) je nacin za dohvacanje
    #rednih brojeva za DOFs koji pripadaju ElementId entitetu
    #npr. ako je ElementId nekakav VOL i koristimo fes.order=1
    #dobit cemo prazan tuple jer su dofsi samo u vrhovima teraedra
#ako zelim dobit sve FE elemente, onda koristim sljedeci pristup
#al on mi se ne svidja jer ne znam kako dohvatiti dofs konkretnog elementa
#bez potrebe da iteriram kroz cijeli iterator

for el in fes.Elements(VOL):  # metoda Elements() unutar klase H1 -> vraca iterator
    print(type(el)) #printa se <class 'ngsolve.comp.FESpaceElement'>. klasa FESpaceElement je dijete klase Ngs_Element
    print (el.dofs) #za npr. fes.order=2, za svaki element (VOL - kodimenzija 0) printaju se liste sa po 10 clanova (tetraedar=4vrha+6bridova) 
    print(fes.GetFE(el))
    print(mesh[el])
    print(mesh[el].nr)

