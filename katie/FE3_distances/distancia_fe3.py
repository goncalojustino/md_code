##IMPORTS##
import statistics as sta

import MDAnalysis as mda
import numpy as np


def plotter_alpha(distances, pos, media):
    '''Irá criar uma pasta que irá conter todos os graficos necessarios para o utilizador de acordo
    com os parametros feitos antriormente.
    Estes gráficos irão ter as distancias ao longo do tempo do movimento da ligação indicada no titulo da imagem
    Irá ter também uma linha que representa a média das distancias'''
    # print(distances)
    # print(len(distances))
    # print(pos)
    import os
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
                                                                                      ",")  # havia graficos a
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


def criacao_ficheiro(tentativa_max_ddistances, ferro, basic, universe):
    '''Irá criar um ficheiro que terá o ferro ligado a um outro atomo do universo com a distancia desejada'''
    with open('Saltbridges_ferro_list.txt', 'w+') as ficheiro:
        for posicao in range(len(basic)):
            distance = np.linalg.norm(
                ferro[0].position - basic[posicao].position)  # distancia de posições
            if distance < tentativa_max_ddistances and distance != 0:  # Irá definoir dizerq eu a distancia entre os atmos terá de quer menor que a distancia que tem no espaço, se isso falhar, então será 'descartada' essa opção
                lista = []
                resids = []
                resname = []
                name = []
                distances = []
                valores = []
                resids.append(basic[posicao].resid)
                resname.append(basic[posicao].resname)
                name.append(basic[posicao].name)
                print(name)
                for ts in universe.trajectory:
                    distances.append((universe.trajectory.time, np.linalg.norm(
                        ferro[0].position - basic[
                            posicao].position)))  # procurar em todos os que passaram a trajetoria ao longo do tempo
                    lista.append(distances)
                    valores.append(np.linalg.norm(ferro[0].position - basic[posicao].position))
                media = sta.mean(valores)  # média que será colocada no gráfico
                frase = str(ferro[0].resname) + ' ' + str(ferro[0].resid) + ' ' + str(
                    ferro[0].name) + ' - ' + str(basic[posicao].resname) + ' ' + str(
                    basic[posicao].resid) + ' ' + str(basic[posicao].name) + ' -- Distancia: ' + str(
                    distance) + ' [Media = ' + str(media) + ']'
                ficheiro.write(frase + '\n')
                distances = np.array(distances)
                plotter_alpha(distances, posicao, media)
                print('plot feito')

    ficheiro.close()
    return resids, resname


# def cenas(universe, ferro, resids):
#     with open('Distanciasas.txt', 'w+') as file:
#         listsa = []
#         for i in resids:
#             i = str(i)
#             frase = "(resid " + i + " )"
#             atom = universe.select_atoms(frase)
#             for ts in universe.trajectory:
#                 listsa.append((universe.trajectory.time, np.linalg.norm(
#                     ferro[0].position - atom[0].position)))
#             dis = np.linalg.norm(ferro[0].position - atom[0].position)  # distancia de posições
#             # listsa.append(dis)
#             frase = str(ferro[0].resname) + ' ' + str(ferro[0].resid) + ' - ' + str(
#                 atom[0].type) + ' ' + str(
#                 atom[0].resid) + ' Distance: ' + str(dis)
#         plotter_alpha(listsa, i)
#         file.write(frase + '\n')
#     file.close()


if __name__ == '__main__':
    # universo = create_universe('./protein_md001_A_2/md2.gro',
    #                            './protein_md001_A_2/md2.xtc')

    universo = create_universe('./md003/md3.gro',
                               './md003/md3.xtc')

    valor, ferroos, basicos = tipo_universo(universo)
    resids, resnames = criacao_ficheiro(valor, ferroos, basicos, universo)
