import numpy, MDAnalysis
from MDAnalysis.analysis import contacts
universe = MDAnalysis.Universe('md.gro','md_fitted.xtc')

#get the residue numbers
piSystems=universe.select_atoms("resname PHE TRP TYR HIS HSD HSE HSP and name CA").resids
###  In [18]: print(piSystems)
###  [  8  14  22  25  45  68  71  84  85  94  96 107 119 128 136 153 154 167 185 186 192 204 207 211 223 238 242 264 273 274 282 285 289 295 300 302 314 317 319]
###  In [21]: print(piSystems[1])
###  14

print(piSystems)

i = 0
while i < len(piSystems):
    print(piSystems[i])
    t = universe.select_atoms("resid "+str(piSystems[i])).resids
    print(t)
    i = i + 1