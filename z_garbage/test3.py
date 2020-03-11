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
        print(piSystems[i].resname," ",piSystems[j].resname) # this gives all pairs
        #get centroid for ring1; get centroid for ring2; compute centroid distance; compute angle
		#below this lines comes from rip-md
        if piSystems[i].resname=="PHE" or piSystems[i].resname=="TYR" or piSystems[i].resname=="TYO":
            piCoords1=[piSystems[i].atoms.CG.position,piSystems[i].atoms.CD1.position,piSystems[i].atoms.CD2.position,piSystems[i].atoms.CE1.position,piSystems[i].atoms.CE2.position,piSystems[i].atoms.CZ.position]
			#computing normal vector
            vector1=piCoords1[0]-piCoords1[5]; vector2=piCoords1[1]-piCoords1[4]; vector3=piCoords1[2]-piCoords1[3]
			prod1=numpy.cross(vector1,vector2); prod2=numpy.cross(vector1,vector3); prod3=numpy.cross(vector2,vector3)	
			normalVector1=(prod1+prod2+prod3)/3
        elif piSystems[i].resname=="TRP":
			piCoords1=[piSystems[i].atoms.CG.position,piSystems[i].atoms.CD1.position,piSystems[i].atoms.CD2.position,piSystems[i].atoms.NE1.position,piSystems[i].atoms.CE2.position,piSystems[i].atoms.CE3.position,piSystems[i].atoms.CZ2.position,piSystems[i].atoms.CZ3.position,piSystems[i].atoms.CH2.position]
			#computing normal vector
			vector1=piCoords1[3]-piCoords1[7]
			vector2=piCoords1[0]-piCoords1[8]
			vector3=piCoords1[1]-((piCoords1[7]+piCoords1[8])/2)
			prod1=numpy.cross(vector1,vector2)
			prod2=numpy.cross(vector1,vector3)
			prod3=numpy.cross(vector2,vector3)
			normalVector1=(prod1+prod2+prod3)/3
		elif piSystems[i].resname=="HIS" or piSystems[i].resname=="HSD" or piSystems[i].resname=="HSE" or piSystems[i].resname=="HSP":
			piCoords1=[piSystems[i].atoms.CG.position,piSystems[i].atoms.ND1.position,piSystems[i].atoms.CE1.position,piSystems[i].atoms.NE2.position,piSystems[i].atoms.CD2.position]
			#computing normal vector			
			vector1=piCoords1[1]-piCoords1[3]
			vector2=piCoords1[0]-piCoords1[2]
			vector3=piCoords1[0]-((piCoords1[2]+piCoords1[3])/2)
			prod1=numpy.cross(vector1,vector2)
			prod2=numpy.cross(vector1,vector3)
			prod3=numpy.cross(vector2,vector3)
			normalVector1=(prod1+prod2+prod3)/3		
	    #only if the first ring is complete we will compute this with the second one
	    if str(piCoords1)!=str(numpy.zeros(3)):
            for coord in piCoords1:
                if str(centroid1)==str(numpy.zeros(3)):
                    centroid1=coord
                else:
				    centroid1+=coord
            centroid1=centroid1/len(piCoords1)

		#geting coords for the second ring
		if piSystems[j].resname=="PHE" or piSystems[j].resname=="TYR" or piSystems[i]=="TYO":
			piCoords2=[piSystems[j].atoms.CG.position,piSystems[j].atoms.CD1.position,piSystems[j].atoms.CD2.position,piSystems[j].atoms.CE1.position,piSystems[j].atoms.CE2.position,piSystems[j].atoms.CZ.position]
			#computing normal vector
			vector1=piCoords2[0]-piCoords2[5]
			vector2=piCoords2[1]-piCoords2[4]
			vector3=piCoords2[2]-piCoords2[3]
			prod1=numpy.cross(vector1,vector2)
			prod2=numpy.cross(vector1,vector3)
			prod3=numpy.cross(vector2,vector3)
			normalVector2=(prod1+prod2+prod3)/3
		elif piSystems[j].resname=="TRP":
			piCoords2=[piSystems[j].atoms.CG.position,piSystems[j].atoms.CD1.position,piSystems[j].atoms.CD2.position,piSystems[j].atoms.NE1.position,piSystems[j].atoms.CE2.position,piSystems[j].atoms.CE3.position,piSystems[j].atoms.CZ2.position,piSystems[j].atoms.CZ3.position,piSystems[j].atoms.CH2.position]
			#computing normal vector
			vector1=piCoords2[3]-piCoords2[7]
			vector2=piCoords2[0]-piCoords2[8]
			vector3=piCoords2[1]-((piCoords2[7]+piCoords2[8])/2)
			prod1=numpy.cross(vector1,vector2)
			prod2=numpy.cross(vector1,vector3)
			prod3=numpy.cross(vector2,vector3)
			normalVector2=(prod1+prod2+prod3)/3
		elif piSystems[j].resname=="HIS" or piSystems[j].resname=="HSD" or piSystems[j].resname=="HSE" or piSystems[j].resname=="HSP":
			piCoords2=[piSystems[j].atoms.CG.position,piSystems[j].atoms.ND1.position,piSystems[j].atoms.CE1.position,piSystems[j].atoms.NE2.position,piSystems[j].atoms.CD2.position]
			#computing normal vector			
			vector1=piCoords2[1]-piCoords2[3]
			vector2=piCoords2[0]-piCoords2[2]
			vector3=piCoords2[0]-((piCoords2[2]+piCoords2[3])/2)
			prod1=numpy.cross(vector1,vector2)
			prod2=numpy.cross(vector1,vector3)
			prod3=numpy.cross(vector2,vector3)
			normalVector2=(prod1+prod2+prod3)/3
		#computing centroid for the second ring
		if str(piCoords2)!=str(numpy.zeros(3)):
			for coord in piCoords2:
				if str(centroid2)==str(numpy.zeros(3)):
					centroid2=coord
				else:
					centroid2+=coord
            centroid2/=len(piCoords2)

		#computing distance to know if it can be an interaction
		dist=distance.euclidean(centroid1,centroid2)
		if dist<=Dist:
			#now we define the dihedral angle and  two distances defined as n and p
			#n represents the distance in A between the origin of the orthogonal system (the ring centroid) and the projection of the mass center of the coupled ring j on the Z-axis of the so defined orthogonal system
			#p is the distance in A between the origin of the orthogonal system and the projection of the mass center of the ring j on the XY plane.
			#calculating angle between 2 normal vector we can get the dihedral angle
			dihedral=Angle.angle_twoVectors(normalVector1, normalVector2)
			#to get n and p we will change the plane of j-ring and we will use pitagoras theorem
			traslated_centroid=centroid2*1#*1 is to get a copy of coords and not a copy of the pointer to the coords
			traslated_centroid[2]=centroid1[2]
				n=distance.euclidean(centroid1, traslated_centroid)
				p=math.sqrt((centroid2[2]-traslated_centroid[2])**2)

				orientation="undefined"
				#now we will define the spatial position
				if (dihedral>=0 and dihedral<30) or (dihedral>=150 and dihedral<180):
					orientation="Parallel orientation"
				if (dihedral>=30 and dihedral<150):
					if p<3.5:
						orientation="T-orientation with the edge to face"
					elif p>=3.5 and n<3:
						orientation="T-orientation with the face to the edge"
					else: #p>=3.5 and n >=3
						orientation="L-orientation"

        #BELOW THIS LINE STAYS
        j = j + 1
