from MDAnalysis.analysis import contacts
import numpy as np
import MDAnalysis as mda
from numpy.linalg import norm
from numpy import linalg as la


def criacao_universo(gro, xtc):
    '''Esta função cria o universo'''
    universe = mda.Universe(gro, xtc)
    return universe


def tipo_universo(universe):
    '''Esta função ira criar as defenições do universo
        Definindo desta foorma acidos e basicos e a distancia maxima que estes podem ter'''
    tentativa_max_distances = 15.0

    acid = universe.select_atoms(
        "(resname ARG and name CZ) ")  # Seleção dos atomos Cz em argininas

    basic = universe.select_atoms(
        "(resname ARG and name CZ )")  # Seleção dos atomos Cz em argininas

    return tentativa_max_distances, acid, basic


def criacao_ficheiro(tentativa_max_distances, acid, basic):
    '''Craiação de um fiucheiro com tudo
        Desde acidos e basicos
        Aqui vamos encontrar o que achamos necessario
        Neste caso tratam-se das ligaçoes arg-arg com uma distancia menor que 15 amestrongs'''

    distancias = []  # distania entre CZ
    resids = [[], []]  # posição [0] sao basicos posicao [1] acidos
    resname = [[], []]  # posição [0] sao nomes de um composto posicao [1] sao nomes de outro composto

    with open('distances_ARG_ARG.txt', 'w+') as ficheiro:  # ira escrever o ficheiro 'distances_ARG_ARG'
        for posicao in range(len(basic)):  # irá verificar quantos basicos existem
            for position in range(len(acid)):  # irá verificar quantos acidos existem
                distance = np.linalg.norm(
                    acid[position].position - basic[posicao].position)  # distancia de posições Cz Cz

                if distance < tentativa_max_distances and distance != 0 and posicao >= position:  # verificar as distancias entre Cz's
                    distancias.append(distance)
                    resids[0].append(basic[posicao].resid)  # Adicionar na posição [0] os Cz convenientes
                    resids[1].append(acid[position].resid)  # Adicionar na posição [1] os Cz convenientes
                    resname[0].append(basic[posicao].resname)
                    resname[1].append(acid[position].resname)
        for numero_basic in range(len(resids[1])):  # ira correr a lista na posição[1]
            for numero_acid in range(len(resids[0])):  # ira correr a lista na posição[0]
                if numero_acid == numero_basic:  # verifica se as posições são iguais, se forem pode selecionar postriormente o nome dos residuos
                    for name_basic in resname[1]:
                        for name_acid in resname[0]:
                            ids = str(name_basic) + ' ' + str(resids[1][numero_basic]) + ' - ' + str(
                                name_acid) + ' ' + str(resids[0][numero_acid]) + ' Distancia: ' + str(
                                distancias[numero_basic]) + '\n'
                    ficheiro.write(ids)  # escreve no ficheiro as ligações pretinetes
    print('resids: ',resids)
    ficheiro.close()
    return resids, distancias, resname


def saber_normal(acid,resids):
    '''Esta função tem como proposito  encontrar os vetores normais dos residuos a volta do CZ das argininas selecionadas postriormente
    ira tambem calcular a media de cada uma destas'''

    medias = [[], []]
    saber_pos_normal = [[], []]
    CZ_valores = [[], []]
    for i, numero in enumerate(resids):
        media = 0
        for selecao in numero:
            stringA = "resid " + str(selecao) + " and name NH1"
            stringB = "resid " + str(selecao) + " and name NE"
            stringC = "resid " + str(selecao) + " and name NH2"
            stringCZ = "resid " + str(selecao) + " and name CZ"


            A = universe.select_atoms(stringA).center_of_geometry()
            B = universe.select_atoms(stringB).center_of_geometry()
            C = universe.select_atoms(stringC).center_of_geometry()
            CZ = universe.select_atoms(stringCZ).center_of_geometry()

            print(selecao)


            Az = acid[0].residue.atoms.NE.position

            print('A',Az)




            arginines = universe.select_atoms("resname ARG and name CZ").residues

            CZ_valores[0].append(selecao)
            CZ_valores[1].append(CZ)



            # fazer os vetores das coordenadas

            AB = B - A
            AC = C - A
            BC = C - B

            # multipçica-se os vetores existentes

            product1 = np.cross(AB, AC)
            product2 = np.cross(AB, BC)
            product3 = np.cross(AC, BC)

            # Calcular a média

            media = (product1 + product2 + product3) / 3  # vetor normal
            print(selecao, 'media: ',media)



            medias[i].append(media)
            saber_pos_normal[i].append(selecao)

            print(stringA, stringB, stringC)
            print('A/B/C', A,B,C)
    print('medias: ',medias)

    # print(saber_pos_normal)


    return medias, CZ_valores,saber_pos_normal


def calcular_alpha(medias):
    '''Calcular o alpha entre 2 vetores de arginina

        Formula:
            cos(alpha)=(vetor u * vetor v)/(norma vetor u * nroma vetor v)<=>
            alpha=arccos((vetor u * vetor v)/(norma vetor u * nroma vetor v))'''


    alpha = []
    alpha1 = []
    alpha2 = []

    for i in medias[1]:
        alpha1.append(i)
    for j in medias[0]:
        alpha2.append(j)



    for contagem in range(len(alpha1)):
        alpha.append(np.rad2deg(np.arccos(((np.dot(alpha1[contagem], alpha2[contagem])) /
                                           (la.norm(alpha1[contagem]) * la.norm(alpha2[contagem]))))))

    print('Valores de alpha', alpha)
    return alpha


def criação_theta(CZ_lista, medias, saber_pos_normal):
    '''Cruia uma função que ira obter o valor de theta.
        Este valor obtem-se a partir do vetor normal ao plano e do vetor da distancia entre CZ '''

    theta1 = []
    theta2 = []
    r = [[],[]]
    theta = []
    for num1, CZ1 in enumerate(CZ_lista[1]):
        for num2, CZ2 in enumerate(CZ_lista[1]):

            if CZ1[0] != CZ2[0] and CZ1[1] != CZ2[1] and CZ1[2] != CZ2[2]:
                discancia_cz = CZ_lista[1][num1] - CZ_lista[1][num2]
                r[0].append(CZ_lista[0][num1])
                frase = 'hi: ', str(CZ_lista[0][num1]), 'hi2: ', str(CZ_lista[0][num2]), 'Distancia: ', str(discancia_cz)
                r[1].append(discancia_cz)
                # print(saber_pos_normal[0][0])

    for i in medias[1]:
        theta1.append(i)
    for j in medias[0]:
        theta2.append(j)
        # print(theta2)

    for contagem in range(len(theta1)):
        theta.append(np.rad2deg(np.arccos(((np.dot(theta1[contagem], r[1][contagem])) / (
                    la.norm(theta1[contagem]) * la.norm(r[1][contagem]))))))






    # to do:
    # arranjar forma de colocar algo do estilo:
    # Angulo theta entre resid x e resid y

    print('Valores de theta', theta)
    return theta

def ficheir_final(resids,resname,alpha,theta):
    with open('ficheiro_final.txt', 'w+') as file:  # ira escrever o ficheiro 'distances_ARG_ARG'
        for numero_basic in range(len(resids[1])):  # ira correr a lista na posição[1]
            for numero_acid in range(len(resids[0])):  # ira correr a lista na posição[0]
                if numero_acid == numero_basic:  # verifica se as posições são iguais, se forem pode selecionar postriormente o nome dos residuos
                    for name_basic in resname[1]:
                        for name_acid in resname[0]:
                            for posa in range(len(alpha)):
                                for post in range(len(theta)):
                                    ids = str(name_basic) + ' ' + str(resids[1][numero_basic]) + ' - ' + str(
                                        name_acid) + ' ' + str(resids[0][numero_acid]) + ' Distancia: ' + str(
                                        distancias[numero_basic]) + ' Angulo aplha: ' + str(alpha[posa]) + ' Angulo Theta: '\
                                        + str(theta[post]) + '\n'

            file.write(ids)  # escreve no ficheiro as ligações pretinetes

    file.close()
    return resids, distancias


if __name__ == '__main__':
    universe = criacao_universo('../../MD simulation/md_0_1.gro',
                                '../../MD simulation/md_0_1.xtc')

    valor, acidos, basicos = tipo_universo(universe)
    resids, distancias, resname = criacao_ficheiro(valor, acidos, basicos)
    medias, CZ, saber_pos_normal = saber_normal(acidos, resids)
    alpha = calcular_alpha(medias)
    theta = criação_theta(CZ, medias, saber_pos_normal)
    ficheir_final(resids, resname, alpha, theta)


