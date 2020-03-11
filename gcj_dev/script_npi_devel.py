import numpy
import MDAnalysis
from numpy import linalg
from numpy import (array, dot, arccos, clip)
from numpy.linalg import norm

from MDAnalysis.analysis import contacts

universe = MDAnalysis.Universe('md.gro','md_fitted.xtc')

#options
max_distance = 3.3
#max n pi distance is the sum of the van der Waals radii of O and C

residues = universe.select_atoms("name O").residues



for i in range(len(residues)):
  for j in range(len(residues)):
    distance = numpy.linalg.norm(residues[i].atoms.C.position - residues[j].atoms.O.position)
    if distance < max_distance and i != j+1 and j != i+1 and i != j:

        # compute liltheta:
        a = residues[i].atoms.O.position #this is the O that defines the npi angle and the pyramidalization degree
        a = residues[i].atoms.O.position #this is the O that defines the npi angle and the pyramidalization degree
        b = residues[i].atoms.C.position
        c = residues[j].atoms.O.position
        ba = a - b
        bc = c - b
        cosine = numpy.dot(ba,bc) / (numpy.linalg.norm(ba)*numpy.linalg.norm(bc))
        angle = numpy.arccos(cosine)
        liltheta = numpy.degrees(angle)

        if liltheta > 98 and liltheta < 120: #makes sense in 109 +/- 10:
          #compute lildel
          d = residues[i].atoms.CA.position #CA of the self residue
          e = residues[i+1].atoms.N.position #N of the next residue
          #a,d,e define the plane to which the distance is computed
          #define the vectors
          v1 = e - a
          v2 = d - a
          #crossproduct vector is normal to plane
          cp = numpy.cross(v1,v2)
          pla,plb,plc = cp
          pld = numpy.dot(cp,e)
          #print('The equation is {0}x + {1}y + {2}z = {3}'.format(pla, plb, plc, pld))
          #compute distance from point b to plane
          lildel = numpy.abs(pla*b[0]+plb*b[1]+plc*b[2]-pld)/numpy.sqrt(pla**2+plb**2+plc**2)
          # print(lildel)

          #compute bigtheta
          normal1=([pla,plb,plc])#this is normal to base plane, already computed
          #compute new plane with atoms b, e and e
          #define the vectors
          v1 = e - b
          v2 = d - b
          #crossproduct vector is normal to plane
          cp = numpy.cross(v1,v2)
          pla2,plb2,plc2 = cp
          pld2 = numpy.dot(cp,e)
          normal2=([pla2,plb2,plc2]) # this is normal to second plane

          cosine = numpy.dot(normal1,normal1) / (numpy.linalg.norm(normal1)*numpy.linalg.norm(normal2))
          angle = numpy.arccos(cosine)
          bigtheta = numpy.degrees(angle)
          # print(bigtheta)

          print("from O in",residues[j].resname,residues[j].resid,"to C in",residues[i].resname,residues[i].resid,":")
          print("C-O distance",distance,"O-C-O angle",liltheta,"displacement from plane",lildel,"pyramidalization angle",bigtheta)



