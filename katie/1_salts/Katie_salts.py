'''
20200314
 GCJ
  plot beauty
  basic extension to include HIS N atoms (as per VMD plugin that has HIS Ns as basic)
DONE
'''


# ----IMPORTS----#
import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis import contacts


# ---STROPMI----#

# ----FUNCTIONS----#
def plotter(distances, my_position):
    import os
    os.makedirs("distances_salt", exist_ok=True)  # Creates folder if it doesn't exist already
    import matplotlib.pyplot as plt
    '''
    Plots the distances along the trajectory
    :param distances: 
        a list containing all the distances of an atom pair between all frames
    :param my_position: 
        the number of the my_positionition of the atoms of the AtomGroup
    :return: 
        a folder containing the graphs of all the bonds that in the end the distance of bond
        is below 4.0 A
    '''
    ax = plt.subplot()
    ax.plot(distances[:, 0], distances[:, 1], '-b', lw=2, label=r"$R_G$")
    ax.set_xlabel("time (ps)", fontsize=16)
    ax.set_ylabel(r"Distance ($\AA$)", fontsize=16)
    ax.set_ylim(0, 16)
    ax.hlines(5, 0, len(distances))
    ax.tick_params(axis='both', which='major', labelsize=12)
    title = "distances_salt/" + acids[my_position].resname + str(acids[my_position].resid) + acids[
        my_position].name \
            + '-' + basics[my_position].resname + str(basics[my_position].resid) + basics[
                my_position].name
    title2 = acids[my_position].resname + str(acids[my_position].resid) + "  " + acids[
        my_position].name \
             + ' - ' + basics[my_position].resname + str(basics[my_position].resid) + "  " + basics[
                 my_position].name
    ax.set_title(title2, fontsize=18)
    ax.figure.savefig(title)
    plt.draw()  # draws graph in the file
    plt.cla()  # clears the axis for next use


# ---- SNIOTCNUF----#

if __name__ == '__main__':
    MAX_DISTANCE = 5.0

    universe = mda.Universe("./md003/md3.gro",
                            "./md003/md3.xtc")  # creates the universe with the gro and xtc made in gromacs

    # selects all the acids
    acid = universe.select_atoms(
        "(resname ASP and name OD*) or (resname GLU and name OE*)")

    # selects all the basics
    basic = universe.select_atoms(
        "(resname ARG and name NH*) or (resname LYS* and name NZ*) or (resname HIS and name ND1) or (resname HIS and name NE2)")

    # create a matrix with all distances between all positions
    distance = mda.analysis.distances.distance_array(acid.positions, basic.positions)

    my_position = [[], []]

    # filtrates the matrix to keep only those with max_distace or less
    for i, line in enumerate(distance):
        for j, column in enumerate(line):
            if column < MAX_DISTANCE and column != 0.0:
                my_position[0].append(i)  # acids
                my_position[1].append(j)  # basics
                print(acid[i].resname, acid[i].resid, acid[i].name, '-', basic[j].resname,
                      basic[j].resid, basic[j].name, 'distance=', round(column, 2))

    # picks up all the residues that had those distances
    # pos => (<pos of acid res>, <pos of basic res>) in the original AtomGroups
    indexes_a = my_position[0]
    indexes_b = my_position[1]

    acids = acid[indexes_a]
    basics = basic[indexes_b]

    global temp
    temp = []

    # calculates the distance over time
    for i in range(len(acids)):
        distances = []
        for trajectory in universe.trajectory:
            distances.append((
                universe.trajectory.time, np.linalg.norm(acids[i].position - basics[i].position)))
        distances = np.array(distances)
        plotter(distances, i)  # plots for every bond
