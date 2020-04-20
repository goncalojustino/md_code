# IMPORTS
import math

import numpy as np

from distance_fe3 import connections, create_universe_pdb, universe_pdb, \
    create_universe, universe_type


def calculate_tau(list1):
    '''
    This function will calculate the Tau through the list that contains all the angles calculated
    in the previous function, it will check which are the two largest existing values in the list
    to make this calculation, thus obtaining a unique value
    :param list1: A list with the values of all formed angles
    :return: The tau value calculated between the two largest angles
    '''
    largest = list1[0]
    largest2 = None
    for item in list1[1:]:
        if item > largest:
            largest2 = largest
            largest = item
        elif largest2 == None or largest2 < item:
            largest2 = item
    tau_six_abs = (360 - largest - largest2) / 180
    tau_six_normf = (360 - math.degrees(np.arccos(5 / 7))) / 180
    tau_six = tau_six_abs / tau_six_normf
    print('tau_six: ', tau_six)
    return tau_six


def aminoacids_selection(resname, iron, basic, resid):
    '''
    This function will check in all the names previously selected if they fit in the written
    restrictions, thus checking their positions to be able to calculate vectors and the angles
     between all relevant amino acids
    :param resname: it's all the names from each residues
    :param iron: all residues with FE3
    :param basic: all residues without Carbon or hidrogen
    :param resids: it's all the ids from each residues
    :return: angles formed by each set of three amino acids, one of them being FE3
    '''
    values = [[], []]
    for j, i in enumerate(resname):
        if i == 'ASP':
            asp = basic[resid[j] - 3].atoms.OD1.position
            values[1].append(asp)
        elif i == 'TYO' or i == 'TYR':
            tyo1 = basic[resid[j] - 3].atoms.OH.position
            values[1].append(tyo1)
        elif i == 'TYO' or i == 'TYR':
            tyo2 = basic[resid[j] - 3].atoms.OH.position
            values[1].append(tyo2)
        elif i == 'HIT' or i == 'HIS':
            hit = basic[resid[j] - 3].atoms.NE2.position
            values[1].append(hit)
        elif i == 'COA':
            coa1 = basic[resid[j] - 9].atoms.O1.position
            values[1].append(coa1)
            coa2 = basic[resid[j] - 9].atoms.O2.position
            values[1].append(coa2)
    fe = iron.atoms.FE.position
    values[0].append(fe)
    angles = []
    
    # Here we will calculate the angles betwen all selected aminoacids
    for position in range(len(values[1])):
        position1 = position + 1
        if position != position1:
            while position1 < len(values[1]):
                vector1 = values[1][position] - values[0][0]
                vector2 = values[1][position1] - values[0][0]
                angle_vectors = np.dot(vector1, vector2) / (
                        np.linalg.norm(vector1) * np.linalg.norm(vector2))
                angle = np.arccos(angle_vectors)
                angles.append(np.degrees(angle))
                position1 = position1 + 1

    return angles


if __name__ == '__main__':
    # UNIVERSE = create_universe('./md003/md3.gro','./md003/md3.xtc')   #Use if using a gro/xtc file
    # CUT_OFF, IRON, BASICS = universe_type(UNIVERSO)     #Use if using a gro/xtc file
    UNIVERSE = create_universe_pdb('./md003/1A8E.pdb')  # Use if using a pdb file
    CUT_OFF, IRON, BASIC = universe_pdb(UNIVERSE)  # Use if using a pdb file
    RESIDS, RESNAMES, NAMES = connections(CUT_OFF, IRON, BASIC, UNIVERSE)
    BASICS = BASIC.residues
    IRONS = IRON.residues
    ANGLES = aminoacids_selection(RESNAMES, IRON, BASICS, RESIDS)
    calculate_tau(ANGLES)
