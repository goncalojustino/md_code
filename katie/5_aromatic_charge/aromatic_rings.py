#IMPORTS
import os
import MDAnalysis as mda
import numpy as np
from numpy import linalg as la
from Look_for_point_charges import dicti, temp

centroid_list = []
all_centroids = []
all_normavectors = []
centroids = {}
normal_vectors = {}


def normal_centroid_angle(distance, name_arom, name_charge, normalVector, vetor_centroid):
    """
    This function will calculate the angles made between the vector normal to the aromatic ring and
    the vector done from the distance between the centroid and the correspondingly charged point,
    It will write everything in a file in the following format:
        Aromatic Aminoacid: <name_arom> - Charged Aminoacid: <name_charge> Distance: <distance>
        Angle: <angles>°
    :param distance: distance between centroid and charged amino acid
    :param name_arom: name of the amino acid with aromatic ring
    :param name_charge: name of the charged amino acid
    :param normalVector: normal vector value
    :param centroid: vector value between centroid and the load
    """
    with open('angle_normal_centroid.txt', 'a+') as file:
        angles = [np.rad2deg(np.arccos((np.dot(normalVector, vetor_centroid)) /
                                       (la.norm(normalVector) * la.norm(vetor_centroid))))]
        my_sentence = 'Aromatic Aminoacid: ' + str(name_arom) + ' - ' + 'Charged Aminoacid: ' + \
                      str(name_charge) + ' Distance:' + str(distance) + ' Angle:' + str(angles)[
                                                                                    1:-1] + \
                      '\xc2\xb0' + '\n'
        file.write(my_sentence + '\n')
    file.close()


def cut_off(distance, name_arom, name_charge, normalVector, subtraction):
    '''
    This function checks if the distance is less than the desired parameters if so, the angle
    between the centroid and the normal vector will be calculated
    :param distance: represents the distance from an aromatic amino acid to a charged one
    :param name_arom: is the name of the aromatic amino acid used
    :param name_charge: is the name of the loaded amino acid used
    :param normalVector: normal vector value
    :param subtraction: vector value between centroid and the load
    '''
    if distance < 6.0:
        normal_centroid_angle(distance, name_arom, name_charge, normalVector, subtraction)


def search(my_number):
    '''
    This function will remove all values from the dictionary that comes from the file
    'Look_for_point_charge.py'
    :param nume: this is the number used to create the universe, so you can search your residue
    :return: the name of the load
    '''
    my_valeus = []
    for value in dicti.values():
        my_valeus.append(value)
    for values1 in my_valeus:
        for values2 in values1:
            if str(my_number) in values2:
                name_charge = values2
                return name_charge


def distances(list_centroids, name_arom, normalVector):
    '''
    Creates two lists to store the keys and the values of those keys so you can use them, you will
    subsequently create a 'phrase' that contains the elements necessary to create a universe each
    Whenever possible, there is a 'search' function for the name of the loaded amino acids
    (values of key) when we verify that the universe exists the code will calculate the distance
    to pass to the cut-off function
    :param list_centroids: these are the values of the centroids of the aromatic amino acids
    :param name_arom: this is the name of the aromatic amino acids already selected
    '''
    my_keys = []
    my_valeus = []
    for key in dicti.keys():
        my_keys.append(key)
    for valeu in dicti.values():
        my_valeus.append(valeu)
    for num in temp:
        for position in range(len(my_keys)):
            my_sentence = "resid " + str(num) + " and name " + str(my_keys[position])
            ind = universe.select_atoms(my_sentence).positions
            name_charge = search(num)
            if str(ind) != "[]":
                if len(np.array(ind)) == 1:
                    subtraction = np.subtract(np.array(ind[0]), np.array(list_centroids))
                    distance = np.linalg.norm(subtraction)
                    cut_off(distance, name_arom, name_charge, normalVector, subtraction)
                if len(np.array(ind)) == 2:
                    subtraction = np.subtract(np.array(ind[1]), np.array(list_centroids))
                    distance = np.linalg.norm(subtraction)
                    cut_off(distance, name_arom, name_charge, normalVector, subtraction)


def delete_existing_files(interactions):
    '''
    This function will delete the file if it already exists, if it does not exist it will create a new one
    :param interactions: these are the selected residues of the universe
    '''
    if os.path.exists("aromaticos.txt"):
        os.remove("aromaticos.txt")
        if os.path.exists("angle_normal_centroid.txt"):
            os.remove("angle_normal_centroid.txt")
            coordinates_phe_tyr_tyo(interactions)
            coordinates_his_hsd_hse_hsp_hit(interactions)
            coordinates_trp(interactions)
        else:
            coordinates_phe_tyr_tyo(interactions)
            coordinates_his_hsd_hse_hsp_hit(interactions)
            coordinates_trp(interactions)
    else:
        coordinates_phe_tyr_tyo(interactions)
        coordinates_his_hsd_hse_hsp_hit(interactions)
        coordinates_trp(interactions)


def six_points(aromatic_coords, aminoacids):
    '''
    This function will treat all amino acids that contain a hexose by adding the information
    desired to an 'aromaticos.txt' outpur file
    :param aromatic_coords: coordinates of the required amino acids belonging to hexose
    :param aminoacids: necessary information about the amino acid in question
    '''
    with open('aromatics.txt', 'a+') as my_file:
        dicti1 = {}
        dicti2 = {}
        vector1 = aromatic_coords[0] - aromatic_coords[5]
        vector2 = aromatic_coords[1] - aromatic_coords[4]
        vector3 = aromatic_coords[2] - aromatic_coords[3]
        meanX = (aromatic_coords[0] + aromatic_coords[5]) / 2
        meanY = (aromatic_coords[1] + aromatic_coords[4]) / 2
        meanZ = (aromatic_coords[2] + aromatic_coords[3]) / 2
        centroid = (meanX + meanY + meanZ) / 3
        prod1 = np.cross(vector1, vector2)
        prod2 = np.cross(vector1, vector3)
        prod3 = np.cross(vector2, vector3)
        normalVector = (prod1 + prod2 + prod3) / 3
        name = str(aminoacids.resid) + str(aminoacids.resname)
        dicti2.update({'aromatic ring: ':'6' , 'centroid ': centroid,
                       'normal vector': normalVector})
        dicti1.update({name: dicti2})
        centroids.update({name: centroid})
        all_centroids.append(centroid)
        normal_vectors.update({name: normalVector})
        all_normavectors.append(normalVector)
        my_file.write(str(dicti1) + '\n')
        distances(centroid, name, normalVector)

    my_file.close()


def five_points(aromatic_coords, aminoacids):
    '''
    This function will treat all amino acids containing a pentose by adding the desired information to an out_ my_file of 'aromaticos.txt'
    : stop aromatic_coords: coordinates of the required amino acids belonging to pentose
    : stop aminoacids: necessary information about the amino acid in question
    '''
    with open('aromaticos.txt', 'a+') as my_file:
        dicti1 = {}
        dicti2 = {}
        vector1 = aromatic_coords[1] - aromatic_coords[3]
        vector2 = aromatic_coords[0] - aromatic_coords[2]
        vector3 = aromatic_coords[0] - ((aromatic_coords[2] + aromatic_coords[3]) / 2)
        meanX = (aromatic_coords[1] + aromatic_coords[3]) / 2
        meanY = (aromatic_coords[0] + aromatic_coords[2]) / 2
        meanZ = (aromatic_coords[0] + ((aromatic_coords[2] + aromatic_coords[3]) / 2)) / 2
        centroid = (meanX + meanY + meanZ) / 3
        prod1 = np.cross(vector1, vector2)
        prod2 = np.cross(vector1, vector3)
        prod3 = np.cross(vector2, vector3)
        normalVector = (prod1 + prod2 + prod3) / 3
        name = str(aminoacids.resid) + str(aminoacids.resname)
        dicti2.update({'aromatic ring: ': '5' , 'centroid ': centroid,
                       'normal vector': normalVector})
        dicti1.update({name: dicti2})
        centroids.update({name: centroid})
        all_centroids.append(centroid)
        normal_vectors.update({name: normalVector})
        all_normavectors.append(normalVector)
        my_file.write(str(dicti1) + '\n')
        distances(centroid, name, normalVector)
    my_file.close()


def nine_points(aromatic_coords, aminoacids):
    '''
    This function will treat all amino acids that contain a hexose by adding the information
    desired to an 'aromaticos.txt' outpur file
    :param aromatic_coords: coordinates of the required amino acids belonging to hexose
    :param aminoacids: necessary information about the amino acid in question
    '''
    with open('aromaticos.txt', 'a+') as my_file:
        dicti1 = {}
        dicti2 = {}
        vector1 = aromatic_coords[3] - aromatic_coords[7]
        vector2 = aromatic_coords[0] - aromatic_coords[8]
        vector3 = aromatic_coords[1] - ((aromatic_coords[7] + aromatic_coords[8]) * 2)
        meanX = (aromatic_coords[3] + aromatic_coords[7]) / 2
        meanY = (aromatic_coords[0] + aromatic_coords[8]) / 2
        meanZ = (aromatic_coords[1] + ((aromatic_coords[7] + aromatic_coords[8]) / 2)) / 2
        centroid = (meanX + meanY + meanZ) / 3
        prod1 = np.cross(vector1, vector2)
        prod2 = np.cross(vector1, vector3)
        prod3 = np.cross(vector2, vector3)
        normalVector = (prod1 + prod2 + prod3) / 3
        name = str(aminoacids.resid) + str(aminoacids.resname)
        dicti2.update({'aromatic ring: ': '9', 'centroid ': centroid,
                       'normal vector': normalVector})
        dicti1.update({name: dicti2})
        centroids.update({name: centroid})
        all_centroids.append(centroid)
        normal_vectors.update({name: normalVector})
        all_normavectors.append(normalVector)
        my_file.write(str(dicti1) + '\n')
        distances(centroid, name, normalVector)
    my_file.close()


def create_universe(gro, xtc):
    '''
    Will create the universe for my_files with the extension gro and xtc
    :param gro: input .gro file
    :param xtc: input .xtc file
    :return: the created universe
    '''
    global universe
    universe = mda.Universe(gro, xtc)
    return universe


def tipo_universo(universe):
    '''
    This function will create the type of universe defining the iron and the other connections
    that interest us
    :param universe: all the information that is on gro and xtc
    :return: the information processed like,all the interactions we choose on the universe
    '''

    interactions = universe.select_atoms("(resname PHE) or (resname TRP)  or (resname TYR)"
                                         " or (resname TYO) or (resname HIS) or (resname HSD)"
                                         " or (resname HSE) or (resname HSP) or "
                                         "(resname HIT)").residues
    return interactions


def coordinates_phe_tyr_tyo(interactions):
    '''
    This function will check all the positions (in this case containing 6 atoms) of the atoms that
    belong by adding them to a list and sending it to another function where it will be treated
    :param interactions: all the information with was chosen on gro and xtc
    '''
    resname = ["PHE", "TYR", "TYO"]
    for pos in resname:
        for aa in interactions:
            if pos in str(aa):
                coordinates = []
                cg = aa.atoms.CG.position
                cd1 = aa.atoms.CD1.position
                cd2 = aa.atoms.CD2.position
                ce1 = aa.atoms.CE1.position
                ce2 = aa.atoms.CE2.position
                cz = aa.atoms.CZ.position
                coordinates.append(cg)
                coordinates.append(cd1)
                coordinates.append(cd2)
                coordinates.append(ce1)
                coordinates.append(ce2)
                coordinates.append(cz)
                six_points(coordinates, aa)


def coordinates_his_hsd_hse_hsp_hit(interactions):
    '''
    This function will check all the positions (in this case containing 5 atoms) of the atoms that
    belong by adding them to a list and sending it to another function where it will be treated
    :param interactions: all the information with was chosen on gro and xtc
    '''
    resname = ["HIS", "HSD", "HSE", "HSP", "HIT"]
    for pos in resname:
        for aa in interactions:
            if pos in str(aa):
                coordinates = []
                cg = aa.atoms.CG.position
                nd1 = aa.atoms.ND1.position
                ce1 = aa.atoms.CE1.position
                ne2 = aa.atoms.NE2.position
                cd2 = aa.atoms.CD2.position
                coordinates.append(cg)
                coordinates.append(nd1)
                coordinates.append(ce1)
                coordinates.append(ne2)
                coordinates.append(cd2)
                five_points(coordinates, aa)


def coordinates_trp(interactions):
    '''
    this function will check all positions (in this case it may contain 5, 6 or 9 atoms) of
    atoms that they belong to by adding them to a list and sending it to another function
    where this will be treated
    :param interactions: all the information with was chosen on gro and xtc
    '''
    for aa in interactions:
        if "TRP" in str(aa):
            coordinates = []
            cg = aa.atoms.CG.position
            cd1 = aa.atoms.CD1.position
            ne1 = aa.atoms.NE1.position
            ce2 = aa.atoms.CE2.position
            cd2 = aa.atoms.CD2.position
            coordinates.append(cg)
            coordinates.append(cd1)
            coordinates.append(ne1)
            coordinates.append(ce2)
            coordinates.append(cd2)
            five_points(coordinates, aa)
    for aa in interactions:
        if "TRP" in str(aa):
            coordinates = []
            cd2 = aa.atoms.CD2.position
            ce3 = aa.atoms.CE3.position
            cz3 = aa.atoms.CZ3.position
            ch2 = aa.atoms.CH2.position
            cz2 = aa.atoms.CZ2.position
            ce2 = aa.atoms.CE2.position
            coordinates.append(cd2)
            coordinates.append(ce3)
            coordinates.append(cz3)
            coordinates.append(ch2)
            coordinates.append(cz2)
            coordinates.append(ce2)
            six_points(coordinates, aa)
    for aa in interactions:
        if "TRP" in str(aa):
            coordinates = []
            cg = aa.atoms.CG.position
            cd1 = aa.atoms.CD1.position
            ne1 = aa.atoms.NE1.position
            ce2 = aa.atoms.CE2.position
            cd2 = aa.atoms.CD2.position
            ce3 = aa.atoms.CE3.position
            cz3 = aa.atoms.CZ3.position
            ch2 = aa.atoms.CH2.position
            cz2 = aa.atoms.CZ2.position
            coordinates.append(cg)
            coordinates.append(cd1)
            coordinates.append(ne1)
            coordinates.append(ce2)
            coordinates.append(cd2)
            coordinates.append(ce3)
            coordinates.append(cz3)
            coordinates.append(ch2)
            coordinates.append(cz2)
            nine_points(coordinates, aa)


def main():
    '''
    It is a function that when called cher will run all the code
    '''
    UNIVERSO = create_universe("./md003/md3.gro",
                               "./md003/md3.xtc")

    INTERACTIONS = tipo_universo(UNIVERSO)

    delete_existing_files(INTERACTIONS)

    list_centroids = []
    list_centroids.append(centroid_list)


if __name__ == '__main__':
    UNIVERSO = create_universe("./md003/md3.gro",
                               "./md003/md3.xtc")

    INTERACTIONS = tipo_universo(UNIVERSO)

    delete_existing_files(INTERACTIONS)

    list_centroids = []
    list_centroids.append(centroid_list)
