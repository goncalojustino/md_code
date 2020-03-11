import numpy
import MDAnalysis
from numpy import linalg
from numpy import (array, dot, arccos, clip)
from numpy.linalg import norm

from MDAnalysis.analysis import contacts

universe = MDAnalysis.Universe('md2.gro','md2_fitted.xtc')

#options
catpi_max_distance = 7.0
#Angles in degrees ranges between the cation mass center and the vector normal to the aromatic ring plane
catpi_angles = [0, 60, 120, 180]
#Angles in degrees ranges that define the planar spatial configuration
catpi_planar_angles = [0,30,150,180]
#Angles in degrees ranges that define the oblique spatial configuration
catpi_oblique_angles = [30,60,120,150]
#Angles in degrees ranges that define the orthogonal spatial configuration
catpi_ortho_angles = [60,120]



###########################################
#get aromatic and cation residues
###########################################
piSystems = universe.select_atoms("resname PHE TRP TYR TYO HIS HSD HSE HSP HIT").residues
cations = universe.select_atoms("resname ARG LYS LYSX HSP" or "resname HIS and name HD1" or "resname HIS and name HE2").residues

###########################################
#do the actual thingies
###########################################
for i in range(len(piSystems)):
	for j in range(len(cations)):
		piCoords1=numpy.zeros(3)
		centroid1=numpy.zeros(3)
		normalVector1=numpy.zeros(3)

		###########################################
		# get aromatic coordinates
		###########################################

		if piSystems[i].resname == "PHE" or piSystems[i].resname == "TYR" or piSystems[i].resname == "TYO":
			piCoords1=[piSystems[i].atoms.CG.position,piSystems[i].atoms.CD1.position,piSystems[i].atoms.CD2.position,piSystems[i].atoms.CE1.position,piSystems[i].atoms.CE2.position,piSystems[i].atoms.CZ.position]
			#computer normal vector to ring plane
			vector1=piCoords1[0]-piCoords1[5]
			vector2=piCoords1[1]-piCoords1[4]
			vector3=piCoords1[2]-piCoords1[3]
			prod1=numpy.cross(vector1,vector2)
			prod2=numpy.cross(vector1,vector3)
			prod3=numpy.cross(vector2,vector3)
			normalVector1=(prod1+prod2+prod3)/3
		elif piSystems[i].resname=="TRP":
			piCoords1=[piSystems[i].atoms.CG.position,piSystems[i].atoms.CD1.position,piSystems[i].atoms.CD2.position,piSystems[i].atoms.NE1.position,piSystems[i].atoms.CE2.position,piSystems[i].atoms.CE3.position,piSystems[i].atoms.CZ2.position,piSystems[i].atoms.CZ3.position,piSystems[i].atoms.CH2.position]
			vector1=piCoords1[3]-piCoords1[7]
			vector2=piCoords1[0]-piCoords1[8]
			vector3=piCoords1[1]-((piCoords1[7]+piCoords1[8])*2)
			prod1=numpy.cross(vector1,vector2)
			prod2=numpy.cross(vector1,vector3)
			prod3=numpy.cross(vector2,vector3)
			normalVector1=(prod1+prod2+prod3)/3
		elif piSystems[i].resname=="HIS" or piSystems[i].resname=="HSD" or piSystems[i].resname=="HSE" or piSystems[i].resname=="HSP" or piSystems[i].resname=="HIT":
			piCoords1=[piSystems[i].atoms.CG.position,piSystems[i].atoms.ND1.position,piSystems[i].atoms.CE1.position,piSystems[i].atoms.NE2.position,piSystems[i].atoms.CD2.position]
			vector1=piCoords1[1]-piCoords1[3]
			vector2=piCoords1[0]-piCoords1[2]
			vector3=piCoords1[0]-((piCoords1[2]+piCoords1[3])/2)
			prod1=numpy.cross(vector1,vector2)
			prod2=numpy.cross(vector1,vector3)
			prod3=numpy.cross(vector2,vector3)
			normalVector1=(prod1+prod2+prod3)/3
		#if there were coordinates, proceed
		if str(piCoords1) != str(numpy.zeros(3)):
			for coord in piCoords1:
				if str(centroid1) == str(numpy.zeros(3)):
					centroid1 = coord
				else:
					centroid1+=coord
			centroid1 = centroid1/len(piCoords1)
			#print(centroid1)

		###########################################
		# get cation coordinates
		###########################################
		pointOfCharge=numpy.zeros(3)
		if cations[j].resname=="ARG":
			try:
				pointOfCharge=cations[j].NH1.position
			except:
				pass
		elif cations[j].resname=="LYS or LYSX":
			try:
				pointOfCharge=cations[j].NZ.position
			except:
				pass
		elif cations[j].resname=="HSP" or cations[j].resname=="HIS":
			try:
				pointOfCharge=cations[j].NE2.postition
			except:
				pass
		if str(pointOfCharge)!=str(numpy.zeros(3)):
			if str(normalVector1)!=str(numpy.zeros(3)):
				
				#computing distance between the charge and the aromatic centroid
				dist = numpy.linalg.norm(pointOfCharge-centroid1)
				
				if dist <= catpi_max_distance:
					with open('CationPi_Interactions.dat','a') as filet:
						#calculate the angle between 2 vectors, the vector between centroid of the ring and cation and the normal vecto
						aux_vector=numpy.subtract(pointOfCharge, centroid1)
						c = dot(aux_vector,normalVector1)/(norm(normalVector1)*norm(aux_vector))
						angle = arccos(clip(c,-1,1)) * 57.2957795

						if (angle >= catpi_angles[0] and angle <= catpi_angles[1]) or (angle >= catpi_angles[2] and angle <= catpi_angles[3]):
							print(piSystems[i].resname,piSystems[i].resid," - ",cations[j].resname,cations[j].resid," cation - pi interaction with distance", dist, "and normal to centroid-charge angle of ",angle, "degrees",file=filet)
						
						# compute dihedrals between the CA of the cation,
						# the mass center of the cation, the centroid of the ring
						# and the atom closer to the point of charge

						# classify interactions based on dihedrals
						last_less_dist = 1000
						closer_atom = None
						for pos in piCoords1:
								aux = numpy.linalg.norm(pointOfCharge-pos)
								if aux < last_less_dist:
										last_less_dist = aux
										closer_atom = pos

						cationCA = cations[j].CA.position
						
						aux_vector2 = numpy.subtract(cationCA, pointOfCharge)
						aux_vector3 = numpy.subtract(closer_atom, centroid1)
						c = dot(aux_vector2,aux_vector3)/(norm(aux_vector2)*norm(aux_vector3))
						dihedral = arccos(clip(c,-1,1)) * 57.2957795

						if (dihedral >= catpi_planar_angles[0] and dihedral <= catpi_planar_angles[1]) or (dihedral >= catpi_planar_angles[2] and dihedral <= catpi_planar_angles[3]):
								print(piSystems[i].resname,piSystems[i].resid," - ",cations[j].resname,cations[j].resid," cation - pi interaction of the PLANAR type with distance", dist, ", normal to centroid-charge angle of ",angle, "degrees and dihedral",dihedral, " degrees",file=filet)
						elif (dihedral >= catpi_oblique_angles[0] and dihedral <= catpi_oblique_angles[1]) or (dihedral >= catpi_oblique_angles[2] and dihedral <= catpi_oblique_angles[3]):
								print(piSystems[i].resname,piSystems[i].resid," - ",cations[j].resname,cations[j].resid," cation - pi interaction of the OBLIQUE type with distance", dist, ", normal to centroid-charge angle of ",angle, "degrees and dihedral",dihedral, " degrees",file=filet)
						elif (dihedral >= catpi_ortho_angles[0] and dihedral <= catpi_ortho_angles[1]):
								print(piSystems[i].resname,piSystems[i].resid," - ",cations[j].resname,cations[j].resid," cation - pi interaction of the ORTHO type with distance", dist, ", normal to centroid-charge angle of ",angle, "degrees and dihedral",dihedral, " degrees",file=filet)