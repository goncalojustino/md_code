#reorganized
#plotter code from Katie
#working on full trajectory for distances below cut-off

import numpy
import MDAnalysis
from numpy import linalg
from numpy import (array, dot, arccos, clip)
from numpy.linalg import norm
import math

from MDAnalysis.analysis import contacts

def plotter(distances, pos):
    import os
    os.makedirs("distances", exist_ok=True)
    import matplotlib.pyplot as plt
    ax = plt.subplot()
    ax.plot(distances[:, 0], distances[:, 1], 'r--', lw=2, label=r"$R_G$")
    ax.set_xlabel("time (ps)")
    ax.set_ylabel(r"Pair distance ($\AA$)")
    #ax.set_ylim(0,16)
    #ax.hlines(5,0,1000)
    title = "distances/" + arginines[i].resname + str(arginines[i].resid) + " (" + arginines[i].resname + ")" \
            + '-' + arginines[j].resname + str(arginines[j].resid) + " (" + arginines[j].resname + ")"
    ax.figure.savefig(title)
    plt.draw()  # draws graph in the file
    plt.cla()  # clears the axis for next use

u = MDAnalysis.Universe('md.gro','md.xtc')

arginines = u.select_atoms("resname ARG and name CZ").residues

for i in range(len(arginines)):
    j = i + 1
    while j < len(arginines):
        distance = numpy.linalg.norm(arginines[i].atoms.CZ.position - arginines[j].atoms.CZ.position)
        if distance < 10:
            #print(arginines[i].resid, arginines[j].resid, distance)
            distances = []
            alphas = []
            betas = []
            thetas = []
            for ts in u.trajectory:
                distances.append((u.trajectory.time, numpy.linalg.norm(arginines[i].atoms.CZ.position - arginines[j].atoms.CZ.position)))
                arg1coords = [arginines[i].atoms.NE.position,arginines[i].atoms.NH1.position,arginines[i].atoms.NH2.position]
                vector1 = arg1coords[0]-arg1coords[1]
                vector2 = arg1coords[0]-arg1coords[2]
                vector3 = arg1coords[2]-arg1coords[1]
                prod1=numpy.cross(vector1,vector2)
                prod2=numpy.cross(vector1,vector3)
                prod3=numpy.cross(vector2,vector3)    
                normalVector1=(prod1+prod2+prod3)/3
                betav1 = arginines[i].atoms.CZ.position - arginines[i].atoms.NE.position
                thetav1 = normalVector1
                arg2coords = [arginines[j].atoms.NE.position,arginines[j].atoms.NH1.position,arginines[j].atoms.NH2.position]
                vector1 = arg2coords[0]-arg2coords[1]
                vector2 = arg2coords[0]-arg2coords[2]
                vector3 = arg2coords[2]-arg2coords[1]
                prod1=numpy.cross(vector1,vector2)
                prod2=numpy.cross(vector1,vector3)
                prod3=numpy.cross(vector2,vector3)    
                normalVector2=(prod1+prod2+prod3)/3
                betav2 = arginines[j].atoms.CZ.position - arginines[j].atoms.NE.position    
                thetav2 = arginines[i].atoms.CZ.position - arginines[j].atoms.CZ.position   
                c = dot(normalVector1,normalVector2)/(norm(normalVector1)*norm(normalVector2))
                alpha = arccos(clip(c,-1,1)) * (180/math.pi)
                alphas.append((u.trajectory.time, alpha))
                d = dot(betav1,betav2)/(norm(betav1)*norm(betav2))
                beta = arccos(clip(d,-1,1)) * (180/math.pi)
                betas.append((u.trajectory.time, beta))
                e = dot(thetav1,thetav2)/(norm(thetav1)*norm(thetav2))
                theta = arccos(clip(e,-1,1)) * (180/math.pi)
                thetas.append((u.trajectory.time, theta))

            distances = numpy.array(distances)
            alphas = numpy.array(alphas)
            betas = numpy.array(betas)
            thetas = numpy.array(thetas)

            plotter(alphas, i)

        j = j + 1

