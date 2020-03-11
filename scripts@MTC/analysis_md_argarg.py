import numpy
import MDAnalysis
from numpy import linalg
from numpy import (array, dot, arccos, clip)
from numpy.linalg import norm

from MDAnalysis.analysis import contacts

universe = MDAnalysis.Universe('md2.gro',)

#options
argarg_max_distances = 5.0

argCZs = universe.select_atoms("resname ARG and name CZ")
with open('ArgArg.dat','a') as filet:
    print("ArgArg Max Distance is defined as ", argarg_max_distances,file=filet)
    for i in range(len(argCZs)):
        j = i + 1
        while j < len(argCZs):
            distance = numpy.linalg.norm(argCZs[i].position-argCZs[j].position)
            print(argCZs[i].resname,argCZs[i].resid,"-",argCZs[j].resname,argCZs[j].resid," interaction with a distance of ", distance,file=filet)
            j = j + 1

