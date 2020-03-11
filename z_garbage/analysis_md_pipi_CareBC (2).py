import MDAnalysis
from MDAnalysis.analysis import contacts
import numpy
import MDAnalysis
from numpy import linalg
from numpy import (array, dot, arccos, clip)
from numpy.linalg import norm
import matplotlib.pyplot as plt
import sys


def pysistems_function(elemento, piCoords, centroid):
  if elemento.resname == "PHE" or elemento.resname == "TYR" or elemento.resname == "TYO":
    piCoords = [elemento.atoms.CG.position, elemento.atoms.CD1.position, elemento.atoms.CD2.position,
                 elemento.atoms.CE1.position, elemento.atoms.CE2.position, elemento.atoms.CZ.position]
    # computer normal vector to ring plane
    vector1 = piCoords[0] - piCoords[5]
    vector2 = piCoords[1] - piCoords[4]
    vector3 = piCoords[2] - piCoords[3]
    normalVector = (numpy.cross(vector1, vector2) +
                     numpy.cross(vector1, vector3) +
                     numpy.cross(vector2, vector3)) / 3

  elif elemento.resname == "TRP":
    piCoords = [elemento.atoms.CG.position, elemento.atoms.CD1.position, elemento.atoms.CD2.position,
                 elemento.atoms.NE1.position, elemento.atoms.CE2.position, elemento.atoms.CE3.position,
                 elemento.atoms.CZ2.position, elemento.atoms.CZ3.position, elemento.atoms.CH2.position]
    vector1 = piCoords[3] - piCoords[7]
    vector2 = piCoords[0] - piCoords[8]
    vector3 = piCoords[1] - ((piCoords[7] + piCoords[8]) * 2)
    normalVector = (numpy.cross(vector1, vector2) +
                     numpy.cross(vector1, vector3) +
                     numpy.cross(vector2, vector3)) / 3

  elif elemento.resname == "HIS" or elemento.resname == "HSD" or elemento.resname == "HSE" or elemento[
    i].resname == "HSP" or elemento.resname == "HIT":
    piCoords = [elemento.atoms.CG.position, elemento.atoms.ND1.position, elemento.atoms.CE1.position,
                 elemento.atoms.NE2.position, elemento.atoms.CD2.position]
    vector1 = piCoords[1] - piCoords[3]
    vector2 = piCoords[0] - piCoords[2]
    vector3 = piCoords[0] - ((piCoords[2] + piCoords[3]) / 2)
    normalVector = (numpy.cross(vector1, vector2) +
                     numpy.cross(vector1, vector3) +
                     numpy.cross(vector2, vector3)) / 3

  else:
    normalVector = ""

  if str(piCoords) != str(numpy.zeros(3)):
    for coord in piCoords:
      if str(centroid) == str(numpy.zeros(3)):
        centroid = coord
      else:
        centroid += coord
    centroid = centroid / len(piCoords)

  return normalVector, centroid



if __name__ == '__main__':
  universe = MDAnalysis.Universe('md.gro', "md.xtc")

  piSystems = universe.select_atoms("resname PHE TRP TYR TYO HIS HSD HSE HSP HIT").residues
  ids = dict()
  for ts in universe.trajectory: # for the trajectory (to count each time)

    for index, element in enumerate(piSystems):
      piCoords1 = numpy.zeros(3)
      piCoords2 = numpy.zeros(3)
      centroid1 = numpy.zeros(3)
      centroid2 = numpy.zeros(3)
      normalVector1 = numpy.zeros(3)
      normalVector2 = numpy.zeros(3)

      if (index+1) < len(piSystems):
        with open('Pi-Pi non-interactions.dat','a+') as noninter:
          with open('Pi-Pi interactions.dat', 'a+') as interactions:
            first_ring = pysistems_function(element, piCoords1, centroid1)
            normalVector1 = first_ring[0]
            centroid1 = first_ring[1]

            second_ring = pysistems_function(piSystems[index+1], piCoords2, centroid2)
            normalVector2 = second_ring[0]
            centroid2 = second_ring[1]

            distance = numpy.linalg.norm(centroid1 - centroid2)
            twodimension_array = dot(normalVector1, normalVector2) / (norm(normalVector1) * norm(normalVector2))
            dihedral = arccos(clip(twodimension_array, -1, 1)) * 57.2957795

            translated_centroid = centroid2 * 1  # *1 is to get a copy of coords and not a copy of the pointer to the coords
            translated_centroid[2] = centroid1[2]
            normalization = numpy.linalg.norm(centroid1 - translated_centroid)
            sqr_root = numpy.sqrt((centroid2[2] - translated_centroid[2]) ** 2)

            orientation = "undefined"

            if distance < 4.4 and ((dihedral >= 0 and dihedral < 30) or (dihedral >= 150 and dihedral < 180)):
              orientation = "parallel"
            elif distance < 5.5:
              if (dihedral >= 30 and dihedral < 150):
                if sqr_root < 3.5:
                  orientation = "T-orientation with the edge to face"
                elif sqr_root >= 3.5 and normalization < 3:
                  orientation = "T-orientation with the face to the edge"
                else:  # p>=3.5 and n >=3
                  orientation = "L-orientation"

            # minority report
            if orientation == "undefined":
              noninter.write("Time: {time}\n"
                             "{elm1_name} {elm1_id} {elm2_name} {elm2_id} : "
                             "distance= {distance} , dihedral= {dihedral}\n".format(time=universe.trajectory.time,
                                                                                    elm1_name=element.resname,
                                                                                  elm1_id=element.resid,
                                                                                  elm2_name=piSystems[index+1].resname,
                                                                                  elm2_id=piSystems[index+1].resid,
                                                                                  distance=distance,
                                                                                  dihedral=dihedral))
              # dummy = 3
            else:          
              elename = str("{elm1_name} {elm1_id} {elm2_name} {elm2_id}".format(elm1_name=element.resname, elm1_id=element.resid, elm2_name=piSystems[index+1].resname, elm2_id=piSystems[index+1].resid))
              if elename not in ids.keys():
                 ids[elename] = [universe.trajectory.time, distance]
              else:
                 ids[elename].append(universe.trajectory.time)
                 ids[elename].append(distance)
              interactions.write("Time: {time}\n"
                                 "{elm1_name} {elm1_id} {elm2_name} {elm2_id} : "
                                 "distance= {distance} , dihedral= {dihedral}\n".format(time=universe.trajectory.time,
                                                                                        elm1_name=element.resname,
                                                                                      elm1_id=element.resid,
                                                                                      elm2_name=piSystems[
                                                                                        index + 1].resname,
                                                                                      elm2_id=piSystems[index + 1].resid,
                                                                                      distance=distance,
                                                                                      dihedral=dihedral))

  for key, value in ids.items():
      times = list()
      dists = list()
      for i in range(0, len(value)-1, 2):
          times.append(value[i])
          dists.append(value[i+1])
      
      plt.figure()
      plt.plot(times, dists, 'r--', lw=2)
      plt.title(key)
      plt.xlabel("time (ps)")
      plt.ylabel("distance")
      plt.savefig("Distance {title} .pdf".format(title=key))
      



