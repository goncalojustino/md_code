import numpy, MDAnalysis
from MDAnalysis.analysis import contacts
universe = MDAnalysis.Universe('md.gro','md_fitted.xtc')

piSystems=universe.select_atoms("resname PHE TRP TYR HIS HSD HSE HSP").residues
print(piSystems)

i = 0
while i < len(piSystems):
    print(piSystems[i])         # this lists the residues
    print(piSystems[i].resname) # this gives only the names
    print(piSystems[i].resid)   # this gives only the numbers
    i = i + 1

piCoords1=numpy.zeros(3)
piCoords2=numpy.zeros(3)
centroid1=numpy.zeros(3)
centroid2=numpy.zeros(3)
normalVector1=numpy.zeros(3)
normalVector2=numpy.zeros(3)

for i in range(len(piSystems)):
    j = i + 1
    while j < len(piSystems):
        print(piSystems[i].resname," ",piSystems[j].resname)
		# this gives all pairs
        #get centroid for ring1; get centroid for ring2; compute centroid distance; compute angle
		#below this lines comes from rip-md
		#BELOW THIS LINE STAYS
		if piSystems[i].resname = "TYR":
			print("TYR")
        j = j + 1
	