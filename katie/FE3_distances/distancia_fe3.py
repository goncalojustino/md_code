##IMPORTS##
import MDAnalysis as mda
import numpy as np


def plotter_alpha(distances, pos):
    print(distances)
    print(pos)
    import os
    os.makedirs("Distancias", exist_ok=True)
    import matplotlib.pyplot as plt
    ax = plt.subplot()
    ax.plot(distances[:, 0], distances[:, 1], 'r--', lw=2, label=r"$R_G$")
    ax.set_xlabel("time (ps)")
    ax.set_ylabel(r"Pair distance ($\AA$)")
    # ax.set_ylim(0, 5)
    # ax.hlines(2.8, 0, 1000)
    ax.set_title("Distancias " + str(acid[0].resname) + "( " + str(acid[0].resid) + " ) - " + str(
        basic[pos].resname) + " (" + str(basic[pos].resid) + " )")
    title = "Distancias/" + str(acid[0].resname) + "(" + str(acid[0].resid) + ")-" + str(
        basic[pos].resname) + "(" + str(basic[pos].resid) + ")"
    ax.figure.savefig(title)
    plt.draw()  # draws graph in the file
    plt.cla()  # clears the axis for next use


def create_universe(gro, xtc):
    universe = mda.Universe(gro, xtc)
    return universe


def tipo_universo(universe):
    global acid, basic
    tentativa_max_ddistances = 2.8

    acid = universe.select_atoms(
        "(resname FE3 ) ")  # criação do universo com FE3

    basic = universe.select_atoms(
        "not resname SOL or not name C* or not name H*")  # criação do universo os restantes residuos

    return tentativa_max_ddistances, acid, basic


def criacao_ficheiro(tentativa_max_ddistances, acid, basic, universe):
    with open('Saltbridges_acid_list.txt', 'w+') as ficheiro:
        for posicao in range(len(basic)):
            distance = np.linalg.norm(
                acid[0].position - basic[posicao].position)  # distancia de posições
            if distance < tentativa_max_ddistances and distance != 0:
                lista = []
                resids = []
                resname = []
                distances = []
                frase = str(acid[0].resname) + ' ' + str(acid[0].resid) + ' - ' + str(
                    basic[posicao].resname) + ' ' + str(basic[posicao].resid) + ' Distance: ' + str(
                    distance)
                print(frase)
                ficheiro.write(frase + '\n')
                resids.append(basic[posicao].resid)
                resname.append(basic[posicao].resname)
                for ts in universe.trajectory:
                    distances.append((universe.trajectory.time, np.linalg.norm(
                        acid[0].position - basic[posicao].position)))
                    # distancia = universe.trajectory.time, acid[0].position - basic[posicao].position
                    lista.append(distances)
                    # lista = np.array(lista)
                    # print(distances)
                distances = np.array(distances)
                plotter_alpha(distances, posicao)

    ficheiro.close()
    return resids, resname


def cenas(universe, acid, resids):
    with open('Distancias.txt', 'w+') as file:
        listsa = []
        for i in resids:
            i = str(i)
            frase = "(resid " + i + " )"
            atom = universe.select_atoms(frase)
            for ts in universe.trajectory:
                listsa.append((universe.trajectory.time, np.linalg.norm(
                    acid[0].position - atom[0].position)))
            dis = np.linalg.norm(acid[0].position - atom[0].position)  # distancia de posições
            # listsa.append(dis)
            frase = str(acid[0].resname) + ' ' + str(acid[0].resid) + ' - ' + str(
                atom[0].type) + ' ' + str(
                atom[0].resid) + ' Distance: ' + str(dis)
        plotter_alpha(listsa, i)
        file.write(frase + '\n')
    file.close()


if __name__ == '__main__':
    universo = create_universe('./protein_md001_A_2/md2.gro',
                               './protein_md001_A_2/md2.xtc')
    valor, acidos, basicos = tipo_universo(universo)
    resids, resnames = criacao_ficheiro(valor, acidos, basicos, universo)
    # cenas(universo,acidos, resids)
