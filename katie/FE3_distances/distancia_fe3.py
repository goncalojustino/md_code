##IMPORTS##
import os
import statistics as sta

import MDAnalysis as mda
import numpy as np


def plotter_alpha(distances, pos, media):
    '''Irá criar uma pasta que irá conter todos os graficos necessarios para o utilizador de acordo
    com os parametros feitos antriormente.
    Estes gráficos irão ter as distancias ao longo do tempo do movimento da ligação indicada no titulo da imagem
    Irá ter também uma linha que representa a média das distancias'''
    os.makedirs("Distancias33", exist_ok=True)
    import matplotlib.pyplot as plt
    # zorder: marca a prioridade entre gráficos, quanto maior o numero mais a frente fica do ecra
    plt.plot(distances[:, 0], distances[:, 1], 'r--', lw=2, label=r"$Variação da Distancia$",
             zorder=1)
    plt.hlines(media, 0, len(distances), 'k', lw=2, label=r"$Media$",
               zorder=2)  # plt.hlines(<valores>, xmininmo, xmaximo, cor e tipo de linha, label

    plt.xlabel("time (ps)")
    plt.ylabel(r"Pair distance ($\AA$)")
    plt.legend(title='Parameter where:')
    plt.title("Distancias " + str(ferro[0].resname) + "( " + str(ferro[0].resid) + " ) - " + str(
        basic[pos].resname) + " ( " + str(basic[pos].resid) + " )")
    title = "Distancias33/" + str(ferro[0].resname) + "(" + str(ferro[0].resid) + ")-" + str(
        basic[pos].resname) + "(" + str(basic[pos].resid) + ")-" + str(media).replace(".",
                                                                                      ",")
    # havia graficos a
    plt.savefig(title)
    plt.draw()  # draws graph in the file
    plt.cla()  # clears the axis for next use


def create_universe(gro, xtc):
    universe = mda.Universe(gro, xtc)
    return universe


def tipo_universo(universe):
    '''irá criar o tipo de universo definindo o ferro e as outras ligações que nos interessão'''
    global ferro, basic
    tentativa_max_ddistances = 2.8

    ferro = universe.select_atoms(
        "(resname FE3 and name FE) ")  # criação do universo com FE3

    basic = universe.select_atoms(
        "name * and (not name C* and not name H*)")  # criação do universo os restantes residuos

    return tentativa_max_ddistances, ferro, basic


def ligações(tentativa_max_ddistances, ferro, basic, universe):
    resids = []
    resname = []
    name = []
    '''Irá criar um ficheiro que terá o ferro ligado a um outro atomo do universo com a 
    distancia desejada'''
    with open('Ligações_com_FE3.txt', 'w+') as ficheiro:
        for posicao in range(len(basic)):
            distance = np.linalg.norm(
                ferro[0].position - basic[posicao].position)  # distancia de posições
            if distance < tentativa_max_ddistances and distance != 0:  # Irá definoir dizerq eu a
                # distancia entre os atmos terá de quer menor que a distancia que tem no espaço, se
                # isso falhar, então será 'descartada' essa opção
                lista = []
                distances = []
                valores = []
                for ts in universe.trajectory:
                    distances.append((universe.trajectory.time, np.linalg.norm(
                        ferro[0].position - basic[
                            posicao].position)))  # procurar em todos os que passaram a trajetoria ao longo do tempo
                    lista.append(distances)
                    valores.append(np.linalg.norm(ferro[0].position - basic[posicao].position))
                media = sta.mean(valores)  # média que será colocada no gráfico
                frase = str(ferro[0].resname) + ' ' + str(ferro[0].resid) + ' ' + str(
                    ferro[0].name) + ' - ' + str(basic[posicao].resname) + ' ' + str(
                    basic[posicao].resid) + ' ' + str(
                    basic[posicao].name) + ' -- Distancia: ' + str(
                    distance) + ' [Media = ' + str(media) + ']'
                ficheiro.write(frase + '\n')
                distances = np.array(distances)
                plotter_alpha(distances, posicao, media)
                resids.append(basic[posicao].resid)
                resname.append(basic[posicao].resname)
                name.append(basic[posicao].name)

    ficheiro.close()
    return resids, resname, name


def linhas_usadas_de_ficheiro(resids, resnames, names):
    os.remove('Distanciasas1.txt')
    chaves = []
    for loc in range(len(resids)):
        if len(names[loc]) == 2:
            chave = ' ' + str(resids[loc]) + str(resnames[loc]) + '     ' + str(names[loc])
            chaves.append(chave)
        if len(names[loc]) == 3:
            chave = ' ' + str(resids[loc]) + str(resnames[loc]) + '    ' + str(names[loc])
            chaves.append(chave)
        if len(names[loc]) == 1:
            chave = ' ' + str(resids[loc]) + str(resnames[loc]) + '      ' + str(names[loc])
            chaves.append(chave)
    chaves.append('339FE3     FE')

    with open('./md003/md3.gro') as f:
        linhas = f.readlines()
        for lines in linhas:
            for pos in range(len(chaves)):
                if chaves[pos] in lines:
                    with open('Distanciasas1.txt', 'a+') as file:
                        file.write(lines + '\n')
                    file.close()
    f.close()


if __name__ == '__main__':
    # universo = create_universe('./protein_md001_A_2/md2.gro',
    #                            './protein_md001_A_2/md2.xtc')

    UNIVERSO = create_universe('./md003/md3.gro',
                               './md003/md3.xtc')

    VALOR, FERROS, BASICOS = tipo_universo(UNIVERSO)
    RESIDS, RESNAMES, NAMES = ligações(VALOR, FERROS, BASICOS, UNIVERSO)
    linhas_usadas_de_ficheiro(RESIDS, RESNAMES, NAMES)
