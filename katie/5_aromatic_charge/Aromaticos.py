import os
import MDAnalysis as mda
import numpy as np
from numpy import linalg as la
from Look_for_point_charges import dicti, temp

centroid_list = []  # esta lista sera criada aqui para poder ser utilizada no codigo inteiro
                    # sem a reeniciar

def angulo_normal_centroid(distance, nome_arom, name_carga, normal, vetor_centroid):
    '''
    Esta funcao ira calcular os angulos feitos entre o vetor normal ao anel aromático e o vetor
    feito a partir da distancia entre o centroid e ao ponto com carga correspondente,
    Ira escrever tudo num ficheiro como seguinte formato:
        Aromatic Aminoacid: <nome_arom> - Charged Aminoacid: <name_carga> Distance: <distance>
        Angle: <angles>°
    :param distance: distancia entre centroide e o aminoacido com carga
    :param nome_arom: nome do aminoacido com anel aromatico
    :param name_carga: nome do aminoacido com carga
    :param normal: valor do vetor normal
    :param centroid: valor do vetor entre centroid e a carga
    '''
    with open('angle_normal_centroid.txt', 'a+') as file:
        angles = []
        angles.append(np.rad2deg(np.arccos((np.dot(normal, vetor_centroid)) /
                                           (la.norm(normal) * la.norm(vetor_centroid)))))
        frase ='Aromatic Aminoacid: ' + str(nome_arom) + ' - '+ 'Charged Aminoacid: ' + \
               str(name_carga) + ' Distance:' + str(distance) + ' Angle:' + str(angles)[1:-1]+\
               '\xc2\xb0'+'\n'
        file.write(frase + '\n')
    file.close()


def cut_off(distance, nome_arom, name_carga, normalVector, subtraction):
    '''
    Esta função irá verificar se a distancia é menor do que os parametros desejados e dar um print
    dos aminoacidos em que esta foi feita com o seu respetivo valor
    :param distance: representa a distancia de um aminoacido aromatico para um com carga
    :param nome_arom: é o nome do aminoacido aromatico utilizado
    :param name_carga: é o nome do aminoacido carregado utilizado
    '''
    if distance < 6.0:
        angulo_normal_centroid(distance, nome_arom, name_carga, normalVector, subtraction)


def procura(nume):
    '''
    Esta função irá retirar todos os valores do dicionário que vem do ficheiro
    'Look_for_point_charge.py'
    :param nume: trata-se do numero utilizado para criar o universo, podendo desta forma procurar
    o seu residuo
    :return: o nome da carga
    '''
    place = []
    for valeu in dicti.values():
        place.append(valeu)
    for i in place:
        for j in i:
            if str(nume) in j:
                name_carga = j
                return name_carga


def distances(lista_centroides, nome_arom, normalVector):
    '''
    cria duas listas para guardar as chaves e os valores dessas chaves para as poder utilizar, irá
    postriormente criar uma 'frase' que contem os elementos necessarios para cria um universo cada
    vez que possivel indo buscar há funcao 'procura' o nome dos aminoacidos carregados (valores da
    key) ao verificarmos que o universo existe o codigo irá calcular a distancia a passar para a
    função cut-off
    :param lista_centroides: trata-se dos valores dos centroides dos aminoacidos aromáticos
    :param nome_arom: trata-se do nome dos aminoacido aromaticos já selecoionados
    '''
    chave = []
    place = []
    for key in dicti.keys():
        chave.append(key)
    for valeu in dicti.values():
        place.append(valeu)
    for num in temp:
        for position in range(len(chave)):
            frase = "resid " + str(num) + " and name " + str(chave[position])
            ind = universe.select_atoms(frase).positions
            name_carga = procura(num)
            if str(ind) != "[]":
                if len(np.array(ind)) == 1:
                    subtraction = np.subtract(np.array(ind[0]), np.array(lista_centroides))
                    distance = np.linalg.norm(subtraction)
                    cut_off(distance, nome_arom, name_carga, normalVector, subtraction)
                if len(np.array(ind)) == 2:
                    subtraction = np.subtract(np.array(ind[1]), np.array(lista_centroides))
                    distance = np.linalg.norm(subtraction)
                    cut_off(distance, nome_arom, name_carga, normalVector, subtraction)


def eliminar_ficheiros_existentes(interactions):
    '''
    Esta função irá eliminar o ficheiro caso ja exista, se não existir irá criar um novo
    :param interactions: trata-se dos residuos selecionados do universo
    '''
    if os.path.exists("aromaticos.txt"):
        os.remove("aromaticos.txt")
        os.remove("angle_normal_centroid.txt")
        coordinates_phe_tyr_tyo(interactions)
        coordinates_his_hsd_hse_hsp_hit(interactions)
        coordinates_trp(interactions)
    else:
        coordinates_phe_tyr_tyo(interactions)
        coordinates_his_hsd_hse_hsp_hit(interactions)
        coordinates_trp(interactions)


def six_points(aromatic_coords, aminoacids):
    '''
    Esta função irá tratar todos os aminoacidos que contem uma hexose adicionando a informação
    desejada a um ficheiro de outpur 'aromaticos.txt'
    :param aromatic_coords: coordenadas dos aminoacidos necessarios e pertencentes a hexose
    :param aminoacids: informação necessária sobre o aminoacido em questão
    '''
    with open('aromaticos.txt', 'a+') as ficheiro:
        dicti1 = {}
        dicti2 = {}
        vector1 = aromatic_coords[0] - aromatic_coords[5]
        vector2 = aromatic_coords[1] - aromatic_coords[4]
        vector3 = aromatic_coords[2] - aromatic_coords[3]
        medioX = (aromatic_coords[0] + aromatic_coords[5]) / 2
        medioY = (aromatic_coords[1] + aromatic_coords[4]) / 2
        medioZ = (aromatic_coords[2] + aromatic_coords[3]) / 2
        centroid = (medioX + medioY + medioZ) / 3
        prod1 = np.cross(vector1, vector2)
        prod2 = np.cross(vector1, vector3)
        prod3 = np.cross(vector2, vector3)
        normalVector = (prod1 + prod2 + prod3) / 3
        nome = str(aminoacids.resid) + str(aminoacids.resname)
        dicti2.update({'anel aromatico: ': aminoacids.resname, 'centroide ': centroid,
                       'normal vector': normalVector})
        dicti1.update({nome: dicti2})
        ficheiro.write(str(dicti1) + '\n')
        distances(centroid, nome, normalVector)

    ficheiro.close()


def five_points(aromatic_coords, aminoacids):
    '''
    Esta função irá tratar todos os aminoacidos que contem uma pentose adicionando a informação desejada a um ficheiro de outpur 'aromaticos.txt'
    :param aromatic_coords: coordenadas dos aminoacidos necessarios e pertencentes a pentose
    :param aminoacids: informação necessária sobre o aminoacido em questão
    '''
    with open('aromaticos.txt', 'a+') as ficheiro:
        dicti1 = {}
        dicti2 = {}
        vector1 = aromatic_coords[1] - aromatic_coords[3]
        vector2 = aromatic_coords[0] - aromatic_coords[2]
        vector3 = aromatic_coords[0] - ((aromatic_coords[2] + aromatic_coords[3]) / 2)
        medioX = (aromatic_coords[1] + aromatic_coords[3]) / 2
        medioY = (aromatic_coords[0] + aromatic_coords[2]) / 2
        medioZ = (aromatic_coords[0] + ((aromatic_coords[2] + aromatic_coords[3]) / 2)) / 2
        centroid = (medioX + medioY + medioZ) / 3
        # print("Centroid : ", (medioX+medioY+medioZ)/3)
        prod1 = np.cross(vector1, vector2)
        prod2 = np.cross(vector1, vector3)
        prod3 = np.cross(vector2, vector3)
        normalVector = (prod1 + prod2 + prod3) / 3
        # print("Vetor normal: ",normalVector, "residuo: ", aminoacids.resid, aminoacids.resname)
        nome = str(aminoacids.resid) + str(aminoacids.resname)
        dicti2.update({'anel aromatico: ': aminoacids.resname, 'centroide ': centroid,
                       'normal vector': normalVector})
        dicti1.update({nome: dicti2})
        # print(dicti1)
        ficheiro.write(str(dicti1) + '\n')
        distances(centroid, nome, normalVector)
    ficheiro.close()
    return centroid_list, nome


def nine_points(aromatic_coords, aminoacids):
    '''
    Esta função irá tratar todos os aminoacidos que contem um eneágono adicionando a informação desejada a um ficheiro de outpur 'aromaticos.txt'
    :param aromatic_coords: coordenadas dos aminoacidos necessarios e pertencentes ao eneágono
    :param aminoacids: informação necessária sobre o aminoacido em questão
    '''
    with open('aromaticos.txt', 'a+') as ficheiro:
        dicti1 = {}
        dicti2 = {}
        vector1 = aromatic_coords[3] - aromatic_coords[7]
        vector2 = aromatic_coords[0] - aromatic_coords[8]
        vector3 = aromatic_coords[1] - ((aromatic_coords[7] + aromatic_coords[8]) * 2)
        medioX = (aromatic_coords[3] + aromatic_coords[7]) / 2
        medioY = (aromatic_coords[0] + aromatic_coords[8]) / 2
        medioZ = (aromatic_coords[1] + ((aromatic_coords[7] + aromatic_coords[8]) / 2)) / 2
        centroid = (medioX + medioY + medioZ) / 3
        # print("Centroid : ", (medioX + medioY + medioZ) / 3)
        prod1 = np.cross(vector1, vector2)
        prod2 = np.cross(vector1, vector3)
        prod3 = np.cross(vector2, vector3)
        normalVector = (prod1 + prod2 + prod3) / 3
        # print("Vetor normal: ",normalVector, "residuo: ", aminoacids.resid, aminoacids.resname)
        nome = str(aminoacids.resid) + str(aminoacids.resname)
        dicti2.update({'anel aromatico: ': aminoacids.resname, 'centroide ': centroid,
                       'normal vector': normalVector})
        dicti1.update({nome: dicti2})
        # print(dicti1)
        ficheiro.write(str(dicti1) + '\n')
        distances(centroid, nome, normalVector)
    ficheiro.close()


def create_universe(gro, xtc):
    '''
    Irá criar o universo para ficheiros com a extenção gro e xtc
    :param gro: ficheiro .gro de input
    :param xtc: ficheiro .xtc de input
    :return: o universo criado
    '''
    global universe
    universe = mda.Universe(gro, xtc)
    return universe


def tipo_universo(universe):
    '''
    irá criar o tipo de universo definindo o ferro e as outras ligações que nos interessão
    :param universe: universo que vem dos ficheiros gro e xtc
    :return: aminoacidos no universo com resnames escolhidos
    '''
    global interactions

    interactions = universe.select_atoms("(resname PHE) or (resname TRP)  or (resname TYR)"
                                         " or (resname TYO) or (resname HIS) or (resname HSD)"
                                         " or (resname HSE) or (resname HSP) or (resname HIT)").residues
    return interactions


def coordinates_phe_tyr_tyo(interactions):
    '''
    esta funcao irá vefificar todas as posiçoes(neste caso podendo conter 5, 6 ou 9 atomos) dos
    atomos que pertencem adicionando-as numa lista e envieando-a para outra função onde esta
    sera tratada
    :param interactions: aminoacidos no universo com resnames escolhidos
    '''
    resname = ["PHE", "TYR", "TYO"]
    for pos in resname:
        for aa in interactions:
            if pos in str(aa):
                coordenadas = []
                cg = aa.atoms.CG.position
                cd1 = aa.atoms.CD1.position
                cd2 = aa.atoms.CD2.position
                ce1 = aa.atoms.CE1.position
                ce2 = aa.atoms.CE2.position
                cz = aa.atoms.CZ.position
                coordenadas.append(cg)
                coordenadas.append(cd1)
                coordenadas.append(cd2)
                coordenadas.append(ce1)
                coordenadas.append(ce2)
                coordenadas.append(cz)
                six_points(coordenadas, aa)


def coordinates_his_hsd_hse_hsp_hit(interactions):
    '''
    esta funcao irá vefificar todas as posiçoes(neste caso contendo 5 atomos) dos atomos que
    pertencem adicionando-as numa lista e envieando-a para outra função onde esta sera tratada
    :param interactions: aminoacidos no universo com resnames escolhidos
    '''
    resname = ["HIS", "HSD", "HSE", "HSP", "HIT"]
    for pos in resname:
        for aa in interactions:
            if pos in str(aa):
                coordenadas = []
                cg = aa.atoms.CG.position
                nd1 = aa.atoms.ND1.position
                ce1 = aa.atoms.CE1.position
                ne2 = aa.atoms.NE2.position
                cd2 = aa.atoms.CD2.position
                coordenadas.append(cg)
                coordenadas.append(nd1)
                coordenadas.append(ce1)
                coordenadas.append(ne2)
                coordenadas.append(cd2)
                five_points(coordenadas, aa)


def coordinates_trp(interactions):
    '''
    esta funcao irá vefificar todas as posiçoes(neste caso podendo conter 5, 6 ou 9 atomos) dos
    atomos que pertencem adicionando-as numa lista e envieando-a para outra função onde esta
    sera tratada
    :param interactions: aminoacidos no universo com resnames escolhidos
    '''
    for aa in interactions:
        if "TRP" in str(aa):
            coordenadas = []
            cg = aa.atoms.CG.position
            cd1 = aa.atoms.CD1.position
            ne1 = aa.atoms.NE1.position
            ce2 = aa.atoms.CE2.position
            cd2 = aa.atoms.CD2.position
            coordenadas.append(cg)
            coordenadas.append(cd1)
            coordenadas.append(ne1)
            coordenadas.append(ce2)
            coordenadas.append(cd2)
            five_points(coordenadas, aa)
    for aa in interactions:
        if "TRP" in str(aa):
            coordenadas = []
            cd2 = aa.atoms.CD2.position
            ce3 = aa.atoms.CE3.position
            cz3 = aa.atoms.CZ3.position
            ch2 = aa.atoms.CH2.position
            cz2 = aa.atoms.CZ2.position
            ce2 = aa.atoms.CE2.position
            coordenadas.append(cd2)
            coordenadas.append(ce3)
            coordenadas.append(cz3)
            coordenadas.append(ch2)
            coordenadas.append(cz2)
            coordenadas.append(ce2)
            six_points(coordenadas, aa)
    for aa in interactions:
        if "TRP" in str(aa):
            coordenadas = []
            cg = aa.atoms.CG.position
            cd1 = aa.atoms.CD1.position
            ne1 = aa.atoms.NE1.position
            ce2 = aa.atoms.CE2.position
            cd2 = aa.atoms.CD2.position
            ce3 = aa.atoms.CE3.position
            cz3 = aa.atoms.CZ3.position
            ch2 = aa.atoms.CH2.position
            cz2 = aa.atoms.CZ2.position
            coordenadas.append(cg)
            coordenadas.append(cd1)
            coordenadas.append(ne1)
            coordenadas.append(ce2)
            coordenadas.append(cd2)
            coordenadas.append(ce3)
            coordenadas.append(cz3)
            coordenadas.append(ch2)
            coordenadas.append(cz2)
            nine_points(coordenadas, aa)


if __name__ == '__main__':
    UNIVERSO = create_universe("./md003/md3.gro",
                               "./md003/md3.xtc")

    INTERACTIONS = tipo_universo(UNIVERSO)

    eliminar_ficheiros_existentes(INTERACTIONS)

    lista_centroides = []
    lista_centroides.append(centroid_list)
