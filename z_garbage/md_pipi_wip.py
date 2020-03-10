import numpy, MDAnalysis
from MDAnalysis.analysis import contacts
universe = MDAnalysis.Universe('md.gro','md_fitted.xtc')

#get the residue numbers
piSystems=universe.select_atoms("resname PHE TRP TYR HIS HSD HSE HSP and name CA").resids
###  In [18]: print(piSystems)
###  [  8  14  22  25  45  68  71  84  85  94  96 107 119 128 136 153 154 167 185 186 192 204 207 211 223 238 242 264 273 274 282 285 289 295 300 302 314 317 319]
###  In [21]: print(piSystems[1])
###  14


i = 0
while i < len(piSystems):
  #print(piSystems[i].resname)
  if piSystems[i].resname == "TRP":
    print("this was a real sexy TRP")
  else:
    print("trash")
  i = i + 1




i = 0
while i < len(piSystems):
  #print(piSystems[i].resname)
  if piSystems[i].resname == "TRP":
    print("this was a real sexy TRP")
  else:
    print("trash")
  i = i + 1
     ...:         





i = 0
while i < len(piSystems):
  if piSystems[i].resname == "TRP":
    print(piSystems[i].resname+" "+piSystems[i].resid)
  elif piSystems[i].resname !== "TRP":
    print("das ist nein TRP!!!")
  i = i + 1

i,j = 0,0

while i < len(piSystems):
  j = i+1

  if a[]




i = 0
while i < len(a):
    item = a[i]
    #res = universe.select_atoms("resid "+str(a[i])+" and name CG").positions
    #print(res)
    picoord1=[universe.select_atoms("resid "+str(a[i])+" and name CG").positions,universe.select_atoms("resid "+str(a[i])+" and name CD1").positions,universe.select_atoms("resid "+str(a[i])+" and name CD2").positions,universe.select_atoms("resid "+str(a[i])+" and name CE1").positions,universe.select_atoms("resid "+str(a[i])+" and name CE2").positions,universe.select_atoms("resid "+str(a[i])+" and name CZ").positions]
    vector1=picoord1[0]-picoord1[5]
    vector2=picoord1[1]-picoord1[4]
    vector3=picoord1[2]-picoord1[3]
    prod1 = numpy.cross(vector1,vector2)
    prod2 = numpy.cross(vector1,vector3)
    prod3 = numpy.cross(vector2,vector3)
    normalvector1=(prod1+prod2+prod3)/3
    print(list(normalvector1))
    i += 1
                                                                                                                                                         

while i < len(a):
  j = i + 1

  if a[i] == "aa":
    test1 = "1aa"
    #print(test1)
  elif a[i] == "bb":
    test1 = "1bb"
    #print(test1)
  elif a[i] == "cc":
    test1 = "1cc"
    #print(test1)
  else:
    test1 = "1else"
    #print(test1) 

  while j < len(a):
    if a[j] == "aa":
      test2 = "2aa"
      #print(test1)
    elif a[j] == "bb":
      test2 = "2bb"
      #print(test2)
    elif a[j] == "cc":
      test2 = "2cc"
      #print(test2)
    else:
      test2 = "2else"
      #print(test2) 
    
    print(test1," ",test2)
      
    
    #elif a[i]!="trp":
    #  print("different")
    j = j +1
  i = i + 1

In [4]: a = 33
   ...: b = 33
   ...: if b > a:
   ...:     print("b > a")
   ...: elif b < a:
   ...:     print("b < a")
   ...: elif b == a:
   ...:     print("b = a")
   ...:
   ...:
b = a




In [9]: while i < len(a):
   ...:     j = i + 1
   ...:     while j < len(a):
   ...:         print(a[i]," ",a[j])
   ...:         j = j +1
   ...:     i = i + 1
   ...:




PHEsystems = universe.select_atoms("resname PHE").residues

pheresidues = universe.select_atoms("resname PHE and name CA")
#In [86]: print(pheresidues.resids)
#[ 22  84  94 107 153 154 167 186 192 204 211 274 282 285 295 302]
i = 0
while i < len(pheresidues):
    print(pheresidues[i].resids)
    i = i +1

In [11]: for acidicAtom in acidic:
         acidic_pos = acidicAtom.position

print(res1)
print(list(res))
print(list(res1))
res1 = universe.select_atoms('resnum 1')
print(list(res1))
for acidicResidue in acidic:
    acidic_post = acidicResidue.cog()
for acidicAtom in acidic:
    acidic_pos = acidicAtom.position
for acidicResidue in acidic:
    acidic_pos = acidicResidue.position
print(list(acidic))
history

In [28]: temp = universe.select_atoms("resid 5")

residue_coms =  numpy.array([r.atoms.center_of_mass()
  for r in universe.select_atoms("protein and resname PHE").residues])