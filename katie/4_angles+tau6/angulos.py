import math

import numpy as np

from distancia_fe3 import ligacoes, create_universe_pdb, tipo_universo_pdb


def cenas(list1):
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
    print("Largest element is:", largest)
    print("Second Largest element is:", largest2)


def wiii(resname, ferro, basico, resid):
    valores = [[], []]
    for j, i in enumerate(resname):
        if i == 'ASP':
            asp = basico[resid[j] - 3].atoms.OD1.position
            valores[1].append(asp)
        elif i == 'TYO' or i == 'TYR':
            tyo1 = basico[resid[j] - 3].atoms.OH.position
            valores[1].append(tyo1)
        elif i == 'TYO' or i == 'TYR':
            tyo2 = basico[resid[j] - 3].atoms.OH.position
            valores[1].append(tyo2)
        elif i == 'HIT' or i == 'HIS':
            hit = basico[resid[j] - 3].atoms.NE2.position
            valores[1].append(hit)
        elif i == 'COA':
            coa1 = basico[resid[j] - 9].atoms.O1.position
            valores[1].append(coa1)
            coa2 = basico[resid[j] - 9].atoms.O2.position
            valores[1].append(coa2)
    fe = ferro.atoms.FE.position
    valores[0].append(fe)

    angulos = []
    for pos in range(len(valores[1])):
        pos1 = pos + 1
        while pos1 < len(valores[1]):
            vector1 = valores[1][pos] - valores[0][0]
            vector2 = valores[1][pos1] - valores[0][0]
            angle_vectors = np.dot(vector1, vector2) / (
                    np.linalg.norm(vector1) * np.linalg.norm(vector2))
            angle = np.arccos(angle_vectors)
            angulos.append(np.degrees(angle))
            pos1 = pos1 + 1
    cenas(angulos)


if __name__ == '__main__':
    # UNIVERSO = create_universe('./md003/md3.gro','./md003/md3.xtc')   #Usar caso se use um Xtc e Gro
    # VALOR, FERROS, BASICOS = tipo_universo_pdb(UNIVERSO)      #Usar caso se use um pdb
    UNIVERSO = create_universe_pdb('./md003/1A8E.pdb')      #Usar caso se use um pdb
    VALOR, FERROS, BASICOS = tipo_universo_pdb(UNIVERSO)        #Usar caso se use um pdb
    RESIDS, RESNAMES, NAMES = ligacoes(VALOR, FERROS, BASICOS, UNIVERSO)
    BASICOS = BASICOS.residues
    FERROS = FERROS.residues
    wiii(RESNAMES, FERROS, BASICOS, RESIDS)
