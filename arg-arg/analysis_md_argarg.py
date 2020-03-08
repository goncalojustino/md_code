#Zhang 2013 ArgArg geometrical parameters
#works fine for ONE structure
#20200304 GCJ

#CZ-CZ distance
#alpha - plane angle defined from NE-NH1-NH2 vectors
#beta - side chain orientation; angle from two CZ−NE vectors indicating the relative orientation of the two guanidinium group
#theta - angular displacement; angle between normal[i] and CZ[i]-to-CZ[j] vectors

#cutoffs for analysis
#CZ-CZ < 5 A / 0.5 nm
#angles - no cutoff


import numpy
import MDAnalysis
from numpy import linalg
from numpy import (array, dot, arccos, clip)
from numpy.linalg import norm

from MDAnalysis.analysis import contacts

u = MDAnalysis.Universe('md.gro','md.xtc')

#options
argarg_max_distances = 5.0

CZs = u.select_atoms("resname ARG and name CZ")
arginines = u.select_atoms("resname ARG and name CZ").residues

for i in range(len(CZs)):
    j = i + 1
    while j < len(CZs):
        distance = numpy.linalg.norm(CZs[i].position-CZs[j].position)
        #the following gives the same results
        distance2 = MDAnalysis.analysis.distances.calc_bonds(CZs[i].position, CZs[j].position,box=u.dimensions)
        
        print(CZs[i].resname,CZs[i].resid,"-",CZs[j].resname,CZs[j].resid," - ", distance,distance2)
        j = j + 1



for i in range(len(arginines)):
    arg1coords=numpy.zeros(3)
    arg2coords=numpy.zeros(3)
    centroid1=numpy.zeros(3)
    centroid2=numpy.zeros(3)
    normalVector1=numpy.zeros(3)
    normalVector2=numpy.zeros(3)

    j = i + 1
    while j < len(arginines):
        arg1coords = [arginines[i].atoms.NE.position,arginines[i].atoms.NH1.position,arginines[i].atoms.NH2.position]
        vector1 = arg1coords[0]-arg1coords[1]
        vector2 = arg1coords[0]-arg1coords[2]
        vector3 = arg1coords[2]-arg1coords[1]
        prod1=numpy.cross(vector1,vector2)
        prod2=numpy.cross(vector1,vector3)
        prod3=numpy.cross(vector2,vector3)    
        normalVector1=(prod1+prod2+prod3)/3
        if str(arg1coords) != str(numpy.zeros(3)):
          for coord in arg1coords:
            if str(centroid1) == str(numpy.zeros(3)):
              centroid1 = coord
            else:
              centroid1+=coord
          centroid1 = centroid1/len(arg1coords)

        betav1 = arginines[i].atoms.CZ.position - arginines[i].atoms.NE.position 
        #folowing zhang 2013 definition
        thetav1 = normalVector1

        arg2coords = [arginines[j].atoms.NE.position,arginines[j].atoms.NH1.position,arginines[j].atoms.NH2.position]
        vector1 = arg1coords[0]-arg1coords[1]
        vector2 = arg1coords[0]-arg1coords[2]
        vector3 = arg1coords[2]-arg1coords[1]
        prod1=numpy.cross(vector1,vector2)
        prod2=numpy.cross(vector1,vector3)
        prod3=numpy.cross(vector2,vector3)    
        normalVector2=(prod1+prod2+prod3)/3
        if str(arg2coords) != str(numpy.zeros(3)):
          for coord in arg2coords:
            if str(centroid1) == str(numpy.zeros(3)):
              centroid1 = coord
            else:
              centroid1+=coord
          centroid1 = centroid1/len(arg1coords)     

        betav2 = arginines[j].atoms.CZ.position - arginines[j].atoms.NE.position    
        thetav2 = arginines[i].atoms.CZ.position - arginines[j].atoms.CZ.position    

        #compute dihedral (angle between vector1 and vector2)
        c = dot(normalVector1,normalVector2)/(norm(normalVector1)*norm(normalVector2))
        alpha = arccos(clip(c,-1,1)) * 57.2957795

        d = dot(betav1,betav2)/(norm(betav1)*norm(betav2))
        beta = arccos(clip(d,-1,1)) * 57.2957795

        e = dot(thetav1,thetav2)/(norm(thetav1)*norm(thetav2))
        theta = arccos(clip(e,-1,1)) * 57.2957795

        print(arginines[i].resname,arginines[i].resid,"-",arginines[j].resname,arginines[j].resid," - ", alpha, beta, theta)



        j = j + 1
