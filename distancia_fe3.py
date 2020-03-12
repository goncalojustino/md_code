##IMPORTS##
import MDAnalysis as mda
import numpy as np


def create_universe(gro,xtc):
    universe = mda.Universe(gro,xtc)
    return universe

def tipo_universo (universe):

    tentativa_max_ddistances = 2.8

    acid = universe.select_atoms(
        "(resname FE3 ) ")  #criação do universo com FE3


    basic = universe.select_atoms(
        "not resname SOL or not name C* or not name H*") #criação do universo os restantes residuos

    return tentativa_max_ddistances, acid, basic


def criacao_ficheiro(tentativa_max_ddistances, acid, basic):
    lista = []
    resids = []
    resname = []

    with open('Saltbridges_acid_list.txt', 'w+') as ficheiro:

        for posicao in range(len(basic)):
            distance = np.linalg.norm(acid[0].position - basic[posicao].position) #distancia de posições
            if distance < tentativa_max_ddistances and distance != 0:
                lista.append(distance)
                resids.append(basic[posicao].resid)
                resname.append(basic[posicao].resname)
                frase = str(acid[0].resname) + ' ' + str(acid[0].resid) + ' - ' + str(basic[posicao].resname) + ' ' + str(basic[posicao].resid) + ' Distance: ' + str(distance)
                ficheiro.write(frase + '\n')
    ficheiro.close()
    return  resids,resname

def cenas(universe,acid, resids):
    with open('Distancias.txt', 'w+') as file:
        listsa = []
        for i in resids:
            i = str(i)
            frase = "(resid "+ i +" )"
            atom = universe.select_atoms(frase)
            dis = np.linalg.norm(acid[0].position - atom[0].position)  # distancia de posições
            listsa.append(dis)
            frase = str(acid[0].resname) + ' ' + str(acid[0].resid) + ' - ' + str(
                atom[0].type) + ' ' + str(
                atom[0].resid) + ' Distance: ' + str(dis)
            file.write(frase + '\n')
    file.close()





if __name__ == '__main__':
    universo = create_universe('./protein_md001_A_2/md2.gro',
                               './protein_md001_A_2/md2.xtc')
    valor, acidos, basicos = tipo_universo(universo)
    resids, resnames = criacao_ficheiro(valor, acidos, basicos)
    cenas(universo,acidos, resids)


