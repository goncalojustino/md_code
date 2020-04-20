##IMPORTS##
import os
import statistics as sta

import MDAnalysis as mda
import numpy as np


def distance_ploter(distances, my_position, mean):
    '''
    It will create a folder that will contain all the graphics needed by the user
    with the parameters made earlier.
    These graphs will have distances over time of the connection movement indicated in the image title
    There will also be a line that represents the average of the distances
    :param distances: alpha distance over time
    :param my_position: position of the aminoacids
    :param mean: the mean of all alpha angles
    '''
    os.makedirs("Distances", exist_ok=True)
    import matplotlib.pyplot as plt
    # zorder:priority mark between graphs, the higher the number the more the front of the screen is
    plt.plot(distances[:, 0], distances[:, 1], 'r--', lw=2, label=r"$Distance variation$",
             zorder=1)
    plt.hlines(mean, 0, len(distances), 'k', lw=2, label=r"$Mean$",
               zorder=2)  # plt.hlines(<valores>, xmininmo, xmaximo, cor e tipo de linha, label

    plt.xlabel("time (ps)",fontsize=16)
    plt.ylabel(r"Pair distance ($\AA$)", fontsize=16)
    plt.legend(title='Parameter where:')
    plt.title("Distances " + str(iron[0].resname) + "( " + str(iron[0].resid) + " ) - " + str(
        basic[my_position].resname) + " ( " + str(basic[my_position].resid) + " )")
    title = "Distances/" + str(iron[0].resname) + "(" + str(iron[0].resid) + ")-" + str(
        basic[my_position].resname) + "(" + str(basic[my_position].resid) + ")-" + str(mean).replace(".",
                                                                                      ",")
    plt.savefig(title)
    plt.draw()  # draws graph in the file
    plt.cla()  # clears the axis for next use


def create_universe(gro, xtc):
    '''
    Will create the universe for gro and xtc files
    :param gro: gro file
    :param xtc: xtc file
    :return: a universe with all the information that is on gro and xtc
    '''
    universe = mda.Universe(gro, xtc)
    return universe


def create_universe_pdb(pdb):
    '''
    Will create the universe for pdb files
    :param pdb: pdb file
    :return: a universe with all the information that is on pdb
    '''
    universe = mda.Universe(pdb)
    return universe


def universe_type(universe):
    '''
    This function will create the type of universe defining the iron and the other connections
    that interest us
    :param universe: all the information that is on gro and xtc
    :return: the information processed like,all the type of FE3 on the universe and the cut-off
    '''
    global iron, basic, cut_off
    cut_off = 2.8

    iron = universe.select_atoms(
        "(resname FE3 and name FE) ")

    basic = universe.select_atoms(
        "name * and (not name C* and not name H*)")


    return cut_off, iron, basic


def universe_pdb(universe):
    '''
    Will create the type of universe defining the iron and the other links that interest us using pdb information
    :param universe: all the information that is on pdb file
    :return: the information processed like,all the tip of FE3 on the universe and the cut-off
    '''
    global iron, basic, tentativa_max_distances
    tentativa_max_distances = 2.8

    iron = universe.select_atoms(
        "(resname FE  or resname FE3 and name FE) ")

    basic = universe.select_atoms(
        "name * and (not name C* and not name H*)")


    return tentativa_max_distances, iron, basic


def connections(cut_off, iron, basic, universe):
    '''
    It will create a file that will have the iron connected to another atom in the universe with the
    desired distance
    :param cut_off: it's the max distance betwen
    :param iron: all residues with FE3
    :param basic: all residues without Carbon or hidrogen
    :param universe: all the information that is on gro and xtc
    :return: all the names and ids for each lap
    '''
    resids = []
    resname = []
    name = []
    with open('ligações_existentes_com_FE3.txt', 'w+') as ficheiro:
        for new_position in range(len(basic)):
            distance = np.linalg.norm(
                iron[0].position - basic[new_position].position)  # distancia de posições
            if distance < cut_off and distance != 0:  # Irá definoir dizerq eu a
                # distancia entre os atmos terá de quer menor que a distancia que tem no espaço, se
                # isso falhar, então será 'descartada' essa opção
                distances = []
                valores = []
                for ts in universe.trajectory:
                    distances.append((universe.trajectory.time, np.linalg.norm(
                        iron[0].position - basic[
                            new_position].position)))  # procurar em todos os que passaram a trajetoria ao longo do tempo
                    valores.append(np.linalg.norm(iron[0].position - basic[new_position].position))
                mean = sta.mean(valores)  # média que será colocada no gráfico
                frase = str(iron[0].resname) + ' ' + str(iron[0].resid) + ' ' + str(
                    iron[0].name) + ' - ' + str(basic[new_position].resname) + ' ' + str(
                    basic[new_position].resid) + ' ' + str(
                    basic[new_position].name) + ' -- Distancia: ' + str(
                    distance) + ' [Mean = ' + str(mean) + ']'
                ficheiro.write(frase + '\n')
                distances = np.array(distances)
                distance_ploter(distances, new_position, mean)
                resids.append(basic[new_position].resid)
                resname.append(basic[new_position].resname)
                name.append(basic[new_position].name)

    ficheiro.close()
    return resids, resname, name


def delete_existing_files_1(resids, resnames, names):
    '''
    This function will check if the file exists or not, if there is one, it will be deleted to be able to be
     written again, otherwise it will just run all the rest of the code
    :param resids: it's all the ids from each residues
    :param resnames: it's all the names from each residues
    :param names: it's all the ids from each residues
    '''
    if os.path.exists("lines_useds_file_gro.txt"):
        os.remove("lines_useds_file_gro.txt")
        file_user_lines(resids, resnames, names)
    else:
        file_user_lines(resids, resnames, names)


def file_user_lines(resids, resnames, names):
    '''
    Write in a separate file with the name 'lines_useds_file_gro' containing all the lines used
    previously in the code
    :param resids: it's all the ids from each residues
    :param resnames: it's all the names from each residues
    :param names: it's all the ids from each residues
    '''
    new_keys = []
    new_keys.append('339FE3     FE')
    for localization in range(len(resids)):
        if len(names[localization]) == 2:
            name_key = ' ' + str(resids[localization]) + str(resnames[localization]) + '     ' + str(names[localization])
            new_keys.append(name_key)
        if len(names[localization]) == 3:
            name_key = ' ' + str(resids[localization]) + str(resnames[localization]) + '    ' + str(names[localization])
            new_keys.append(name_key)
        if len(names[localization]) == 1:
            name_key = ' ' + str(resids[localization]) + str(resnames[localization]) + '      ' + str(names[localization])
            new_keys.append(name_key)
    with open('./md003/md3.gro') as f:
        linhas = f.readlines()
        for new_lines in linhas:
            for my_position in range(len(new_keys)):
                if new_keys[my_position] in new_lines:
                    with open('lines_useds_file_gro.txt', 'a+') as file:
                        file.write(new_lines)
                    file.close()
    f.close()


if __name__ == '__main__':
    # universo = create_universe('./protein_md001_A_2/md2.gro',
    #                            './protein_md001_A_2/md2.xtc')

    UNIVERSO = create_universe('./md003/md3.gro',
                               './md003/md3.xtc')

    CUT_OFF, IRON, BASICS = universe_type(UNIVERSO)

    RESIDS, RESNAMES, NAMES = connections(CUT_OFF, IRON, BASICS, UNIVERSO)
    delete_existing_files_1(RESIDS, RESNAMES, NAMES)
