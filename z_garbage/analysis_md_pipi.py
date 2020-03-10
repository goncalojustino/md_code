import numpy
import MDAnalysis
from numpy import linalg
from numpy import (array, dot, arccos, clip)
from numpy.linalg import norm

from MDAnalysis.analysis import contacts

universe = MDAnalysis.Universe('md2.gro','md2_fitted.xtc')

piSystems = universe.select_atoms("resname PHE TRP TYR TYO HIS HSD HSE HSP HIT").residues
#print(list(piSystems))

#i = 0
#while i < len(piSystems):
#  print(piSystems[i])             # this lists the residues
#  print(piSystems[i].resname)     # this gives the names
#  print(piSystems[i].resid)       # this gives the numbers
#  i = i + 1


for i in range(len(piSystems)):
  piCoords1=numpy.zeros(3)
  piCoords2=numpy.zeros(3)
  centroid1=numpy.zeros(3)
  centroid2=numpy.zeros(3)
  normalVector1=numpy.zeros(3)
  normalVector2=numpy.zeros(3)

  j = i + 1
  while j < len(piSystems):
    with open('Pi-Pi non-interactions.dat','a') as filebad:
      with open('Pi-Pi interactions.dat','a') as filegood:
        #print(piSystems[i].resname," ",piSystems[j].resname) #this gives all pairs
        #RIP-MD PI-PI BLOCK START
        #getting coordinates for the first ring
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
        #getting coordinates for the second ring
        if piSystems[j].resname == "PHE" or piSystems[j].resname == "TYR" or piSystems[j].resname == "TYO":
          piCoords2=[piSystems[j].atoms.CG.position,piSystems[j].atoms.CD1.position,piSystems[j].atoms.CD2.position,piSystems[j].atoms.CE1.position,piSystems[j].atoms.CE2.position,piSystems[j].atoms.CZ.position]
          #computer normal vector to ring plane
          vector1=piCoords2[0]-piCoords2[5]
          vector2=piCoords2[1]-piCoords2[4]
          vector3=piCoords2[2]-piCoords2[3]
          prod1=numpy.cross(vector1,vector2)
          prod2=numpy.cross(vector1,vector3)
          prod3=numpy.cross(vector2,vector3)
          normalVector2=(prod1+prod2+prod3)/3
        elif piSystems[j].resname=="TRP":
          piCoords2=[piSystems[j].atoms.CG.position,piSystems[j].atoms.CD1.position,piSystems[j].atoms.CD2.position,piSystems[j].atoms.NE1.position,piSystems[j].atoms.CE2.position,piSystems[j].atoms.CE3.position,piSystems[j].atoms.CZ2.position,piSystems[j].atoms.CZ3.position,piSystems[j].atoms.CH2.position]
          vector1=piCoords2[3]-piCoords2[7]
          vector2=piCoords2[0]-piCoords2[8]
          vector3=piCoords2[1]-((piCoords2[7]+piCoords2[8])*2)
          prod1=numpy.cross(vector1,vector2)
          prod2=numpy.cross(vector1,vector3)
          prod3=numpy.cross(vector2,vector3)
          normalVector2=(prod1+prod2+prod3)/3
        elif piSystems[j].resname=="HIS" or piSystems[j].resname=="HSD" or piSystems[j].resname=="HSE" or piSystems[j].resname=="HSP" or piSystems[j].resname=="HIT":
          piCoords2=[piSystems[j].atoms.CG.position,piSystems[j].atoms.ND1.position,piSystems[j].atoms.CE1.position,piSystems[j].atoms.NE2.position,piSystems[j].atoms.CD2.position]
          vector1=piCoords2[1]-piCoords2[3]
          vector2=piCoords2[0]-piCoords2[2]
          vector3=piCoords2[0]-((piCoords1[2]+piCoords2[3])/2)
          prod1=numpy.cross(vector1,vector2)
          prod2=numpy.cross(vector1,vector3)
          prod3=numpy.cross(vector2,vector3)
          normalVector2=(prod1+prod2+prod3)/3
        #if there were coordinates, proceed
        if str(piCoords2) != str(numpy.zeros(3)):
          for coord in piCoords2:
            if str(centroid2) == str(numpy.zeros(3)):
              centroid2 = coord
            else:
              centroid2+=coord
          centroid2 = centroid2/len(piCoords1)  
          #print(centroid2)
          #compute distance between centroids
          dist = numpy.linalg.norm(centroid1-centroid2)
          #compute dihedral (angle between vector1 and vector2)
          c = dot(normalVector1,normalVector2)/(norm(normalVector1)*norm(normalVector2))
          dihedral = arccos(clip(c,-1,1)) * 57.2957795

          #n represents the distance in A between the origin of the orthogonal system (the ring centroid) and the projection of the mass center of the coupled ring j on the Z-axis of the so defined orthogonal system
          #p is the distance in A between the origin of the orthogonal system and the projection of the mass center of the ring j on the XY plane.
          #to get n and p we will change the plane of j-ring and we will use pitagoras theorem
          translated_centroid=centroid2*1     #*1 is to get a copy of coords and not a copy of the pointer to the coords
          translated_centroid[2]=centroid1[2]
          n = numpy.linalg.norm(centroid1-translated_centroid)
          p=numpy.sqrt((centroid2[2]-translated_centroid[2])**2)

          orientation="undefined"

          if dist < 4.4 and ((dihedral>=0 and dihedral<30) or (dihedral>=150 and dihedral<180)):
            orientation="parallel"
          elif dist < 5.5:
            if (dihedral>=30 and dihedral<150):
              if p<3.5:
                orientation="T-orientation with the edge to face"
              elif p>=3.5 and n<3:
                orientation="T-orientation with the face to the edge"
              else: #p>=3.5 and n >=3
                orientation="L-orientation"
          #minority report
          if orientation == "undefined":
            print(piSystems[i].resname,piSystems[i].resid,piSystems[j].resname,piSystems[j].resid,": distance =",dist,", dihedral =", dihedral,file=filebad)
            #dummy = 3
          else:
            print(piSystems[i].resname,piSystems[i].resid,piSystems[j].resname,piSystems[j].resid,": ",orientation,", centroid distance =",dist,", dihedral = ",dihedral,file=filegood)
        #RIP-MD PI-PI BLOCK END
        j = j + 1


