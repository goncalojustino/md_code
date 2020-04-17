'''
20200314
 GCJ
  plot beauty
DONE
'''

#reorganized
#plotter code from Katie
#working on full trajectory for distances below cut-off

import numpy
import MDAnalysis
from numpy import (dot, arccos, clip)
from numpy.linalg import norm
import math

def alpha_ploter(distances, my_position):
    '''
    This function will create some plots for different alphas
    :param distances: a list containing all the distances of an atom pair between all frames
    :param my_position: the number of the my_ position of the atoms of the AtomGroup
    :return: a folder containing the graphs of all the bonds that in the end the distance of
            bond are below 15.0 A
    '''
    import os
    os.makedirs("α_Plane_Angle", exist_ok=True)
    import matplotlib.pyplot as plt
    ax = plt.subplot()
    ax.plot(distances[:, 0], distances[:, 1], '-b', lw=2, label=r"$R_G$")
    ax.set_xlabel("time (ps)",fontsize=16)
    ax.set_ylabel(r"α (°)",fontsize=16)
    ax.set_ylim(0,180)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.set_title('Plane Angle (α):' + arginines[my_position].resname + str(arginines[my_position].resid) \
            + ' - ' + arginines[temporary].resname + str(arginines[temporary].resid),fontsize=16)
    title = "α_Plane_Angle/" + arginines[my_position].resname + str(arginines[my_position].resid) + "_(" + arginines[my_position].resname + ")" \
            + '-' + arginines[temporary].resname + str(arginines[temporary].resid) + "_(" + arginines[temporary].resname + ")"
    ax.figure.savefig(title)
    plt.draw()  # draws graph in the file
    plt.cla()  # clears the axis for next use

def theta_ploter(distances, my_position):
    '''
    This function will create some plots for different theta
    :param distances: a list containing all the distances of an atom pair between all frames
    :param my_position: the number of the my_ position of the atoms of the AtomGroup
    :return: a folder containing the graphs of all the bonds that in the end the distance of
            bond are below 15.0 A
    '''
    import os
    os.makedirs("θ_Angular_Displacement", exist_ok=True)
    import matplotlib.pyplot as plt
    ax = plt.subplot()
    ax.plot(distances[:, 0], distances[:, 1], '-b', lw=2, label=r"$R_G$")
    ax.set_xlabel("time (ps)",fontsize=16)
    ax.set_ylabel(r"θ (°)",fontsize=16)
    ax.set_ylim(0,180)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.set_title('Angular Displacement (θ) ' + arginines[my_position].resname + str(arginines[my_position].resid) \
            + ' - ' + arginines[temporary].resname + str(arginines[temporary].resid),fontsize=18)
    title = "θ_Angular_Displacement/" + arginines[my_position].resname + str(arginines[my_position].resid) + "_(" + arginines[
        my_position].resname + ")" \
            + '-' + arginines[temporary].resname + str(arginines[temporary].resid) + "_(" + arginines[
                temporary].resname + ")"
    ax.figure.savefig(title)
    plt.draw()  # draws graph in the file
    plt.cla()  # clears the axis for next use

def beta_ploter(distances, my_position):
    '''
    This function will create some plots for different beta
    :param distances: a list containing all the distances of an atom pair between all frames
    :param my_position: the number of the my_ position of the atoms of the AtomGroup
    :return: a folder containing the graphs of all the bonds that in the end the distance of
            bond are below 15.0 A
    '''
    import os
    os.makedirs("β_Side_Chain_Angle", exist_ok=True)
    import matplotlib.pyplot as plt
    ax = plt.subplot()
    ax.plot(distances[:, 0], distances[:, 1], '-b', lw=2, label=r"$R_G$")
    ax.set_xlabel("time (ps)",fontsize=16)
    ax.set_ylabel(r"β (°)",fontsize=16)
    ax.set_ylim(0,180)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.set_title('Side-chain angle (β) '  + arginines[my_position].resname + str(arginines[my_position].resid) \
            + ' - ' + arginines[temporary].resname + str(arginines[temporary].resid),fontsize=18)
    title = "β_Side_Chain_Angle/" + arginines[my_position].resname + str(arginines[my_position].resid) + "_(" + arginines[my_position].resname + ")" \
            + '-' + arginines[temporary].resname + str(arginines[temporary].resid) + "_(" + arginines[temporary].resname + ")"
    ax.figure.savefig(title)
    plt.draw()  # draws graph in the file
    plt.cla()  # clears the axis for next use

def distance_plotter(distances, my_position):
    '''
    This function will create some plots for different distances
    :param distances: a list containing all the distances of an atom pair between all frames
    :param my_position: the number of the my_ position of the atoms of the AtomGroup
    :return: a folder containing the graphs of all the bonds that in the end the distance of
            bond are below 15.0 A
    '''
    import os
    os.makedirs("distances_Cζ_Cζ", exist_ok=True)
    import matplotlib.pyplot as plt
    ax = plt.subplot()
    ax.plot(distances[:, 0], distances[:, 1], '--b', lw=2, label=r"$R_G$")
    ax.set_xlabel("time (ps)",fontsize=16)
    ax.set_ylabel(r"CZ-CZ distance ($\AA$)",fontsize=16)
    ax.set_ylim(0,20)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.hlines(5,0,len(distances)+5)
    ax.set_title('Cζ distance ' + arginines[my_position].resname + str(arginines[my_position].resid) + \
        ' - ' + arginines[temporary].resname + str(arginines[temporary].resid),fontsize=18)
    title = "distances_Cζ_Cζ/" + arginines[my_position].resname + str(arginines[my_position].resid) + "_(" + arginines[my_position].resname + ")" \
            + '-' + arginines[temporary].resname + str(arginines[temporary].resid) + "_(" + arginines[temporary].resname + ")"
    ax.figure.savefig(title)
    plt.draw()  # draws graph in the file
    plt.cla()  # clears the axis for next use

if __name__ == '__main__':

    u = MDAnalysis.Universe("./md003/md3.gro",
                            "./md003/md3.xtc")

    arginines = u.select_atoms("resname ARG and name CZ").residues

    for my_position in range(len(arginines)):
        temporary = my_position + 1
        while temporary < len(arginines):
            distance = numpy.linalg.norm(arginines[my_position].atoms.CZ.position - arginines[temporary].atoms.CZ.position)
            if distance < 15:
                distances = []
                alphas = []
                betas = []
                thetas = []
                for trajectories in u.trajectory:
                    distances.append((u.trajectory.time, numpy.linalg.norm(arginines[my_position].atoms.CZ.position - arginines[temporary].atoms.CZ.position)))
                    arg1coords = [arginines[my_position].atoms.NE.position,arginines[my_position].atoms.NH1.position,arginines[my_position].atoms.NH2.position]
                    vector1 = arg1coords[0]-arg1coords[1]
                    vector2 = arg1coords[0]-arg1coords[2]
                    vector3 = arg1coords[2]-arg1coords[1]
                    prod1=numpy.cross(vector1,vector2)
                    prod2=numpy.cross(vector1,vector3)
                    prod3=numpy.cross(vector2,vector3)
                    normalVector1=(prod1+prod2+prod3)/3
                    betav1 = arginines[my_position].atoms.CZ.position - arginines[my_position].atoms.NE.position
                    thetav1 = normalVector1
                    arg2coords = [arginines[temporary].atoms.NE.position,arginines[temporary].atoms.NH1.position,arginines[temporary].atoms.NH2.position]
                    vector1 = arg2coords[0]-arg2coords[1]
                    vector2 = arg2coords[0]-arg2coords[2]
                    vector3 = arg2coords[2]-arg2coords[1]
                    prod1=numpy.cross(vector1,vector2)
                    prod2=numpy.cross(vector1,vector3)
                    prod3=numpy.cross(vector2,vector3)
                    normalVector2=(prod1+prod2+prod3)/3
                    betav2 = arginines[temporary].atoms.CZ.position - arginines[temporary].atoms.NE.position
                    thetav2 = arginines[my_position].atoms.CZ.position - arginines[temporary].atoms.CZ.position
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


                alpha_ploter(alphas, my_position)
                theta_ploter(thetas, my_position)
                beta_ploter(betas, my_position)
                distance_plotter(distances, my_position)

            temporary = temporary + 1
