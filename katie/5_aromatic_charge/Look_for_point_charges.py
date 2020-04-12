import MDAnalysis as mda

universe = mda.Universe("./md003/md3.gro", "./md003/md3.xtc")
interactions = universe.select_atoms("(resname ARG and name NH1) or (resname ASP and name OD*) "
                                     "or (resname GLU and name OE*) or (resname HIS and name ND1)"
                                     "or (resname HIT and name NE2) or (resname LYS* and name NZ)"
                                     "or (resname TYO and name OH)")
ferro = universe.select_atoms("(resname FE3 and name FE) or (resname V3 and name V3)"
                              "or (resname VO and name V4) or (resname VO2 and name V5)"
                              "or (resname VO4 and name V5)")
lista = []  # criação de uma lista vazia
ids = []
dicti = {}  # criação de um dicionário

# Seleção dos NH1 das ARG

for line in interactions:  # correr por todas as interactions obtendo cada linha
    if 'NH1' in str(line):  # verificar de o atmo que pretendemos encontra-se
        frase = "resid " + str(line.resid) + " and (name HH11 or name HH12)"  ## criação de uma variavel para se poder criar uma seleção de um universo 'interativo'
        select = universe.select_atoms(frase)  ## criação de uma nova seleção no universo utilizando a variavel definida acima
        if str(select) != "<AtomGroup []>":  # esta funcao ira verificar se o HH12/HH11 existe, e se isso acontecer mante-lo
            juncao = str(line.resid) + str(line.resname)    ##irá concatenar os id's e os nomes dos residuos em questão
            ids.append(line.resid)
            lista.append(juncao)  ## Irá guardar todos os id's que  interessam numa lista
            dicti.update({'NH1': lista})  ## guardar com o indice desejado (neste caso trata-se do nome do atomo selecionado) com ids desejados

# Seleção dos OD* das ASP

lista = []  # necessario para limpar a lista
for line in interactions:
    if 'ASP' in str(line):
        if line.resid not in lista:
            juncao = str(line.resid) + str(line.resname)
            ids.append(line.resid)
            lista.append(juncao)
            dicti.update({'OD*': lista})

# Seleção dos OE* das GLU

lista = []
for line in interactions:
    if 'GLU' in str(line):
        if line.resid not in lista:
            juncao = str(line.resid) + str(line.resname)
            ids.append(line.resid)
            lista.append(juncao)
            dicti.update({'OE*': lista})

# Seleção dos ND1 das HIS

lista = []
for line in interactions:
    if 'ND1' in str(line):
        frase = "resid " + str(line.resid) + " and name HD1"
        select = universe.select_atoms(frase)
        for line1 in range(len(select)):  # apos a primeira seleção iremos necessitar de outra pois neste caso precisamos que ambos os requisitos existam
            frase1 = "resid " + str(select[line1].resid) + " and name HE2"
            select1 = universe.select_atoms(frase1)
            if str(select1) != "<AtomGroup []>":  # esta funcao ira verificar se existe algo em select1, e se isso acontecer mantem-no
                juncao = str(select1[line1].resid) + str(select1[line1].resname)
                ids.append(line.resid)
                lista.append(juncao)
                dicti.update({'ND1': lista})

# Seleção dos NE2 das HIT

lista = []
for line in interactions:
    if 'NE2' in str(line):
        if line.resid not in lista:
            juncao = str(line.resid) + str(line.resname)
            ids.append(line.resid)
            lista.append(juncao)
            dicti.update({'NE2': lista})

# Seleção dos HZ3 das LYS+

lista = []
lys = ['LYS', 'LYSX', 'LYSY']  # lista com lisinas existentes no gro
for pos in lys:  # correr na lista de lys e selecionar por cada volta um dos nomes
    for line in interactions:
        if pos in str(line):  # verificar se cada nome selecionado na pos existe no ficheiro
            frase = "resid " + str(line.resid) + " and name HZ3"
            select = universe.select_atoms(frase)
            if str(select) != "<AtomGroup []>":
                juncao = str(line.resid) + str(line.resname)
                ids.append(line.resid)
                lista.append(juncao)
                dicti.update({'HZ3': lista[:-1]})

lista = []
for line in interactions:
    if 'OH' in str(line):
        frase = "resid " + str(line.resid) + " and name OH"
        select = universe.select_atoms(frase)
        if str(select) != "<AtomGroup []>":
            juncao = str(line.resid) + str(line.resname)
            ids.append(line.resid)
            lista.append(juncao)
            dicti.update({'OH': lista})

lista = []
for line in ferro:
    if 'FE3' in str(line):
        frase = "resid " + str(line.resid) + " and name FE"
        select = universe.select_atoms(frase)
        if str(select) != "<AtomGroup []>":
            juncao = str(line.resid) + str(line.resname)
            ids.append(line.resid)
            lista.append(juncao)
            dicti.update({'FE': lista})


lista = []
for line in ferro:
    if 'V3' in str(line):
        frase = "resid " + str(line.resid) + " and name V3"
        select = universe.select_atoms(frase)
        if str(select) != "<AtomGroup []>":
            juncao = str(line.resid) + str(line.resname)
            ids.append(line.resid)
            lista.append(juncao)
            dicti.update({'V3': lista})

lista = []
for line in ferro:
    if 'VO' in str(line):
        frase = "resid " + str(line.resid) + " and name V4"
        select = universe.select_atoms(frase)
        if str(select) != "<AtomGroup []>":
            juncao = str(line.resid) + str(line.resname)
            ids.append(line.resid)
            lista.append(juncao)
            dicti.update({'V4': lista})

lista = []
for line in ferro:
    if 'VO2' in str(line):
        frase = "resid " + str(line.resid) + " and name V5"
        select = universe.select_atoms(frase)
        if str(select) != "<AtomGroup []>":
            juncao = str(line.resid) + str(line.resname)
            ids.append(line.resid)
            lista.append(juncao)
            dicti.update({'V5': lista})

lista = []
for line in ferro:
    if 'VO4' in str(line):
        frase = "resid " + str(line.resid) + " and name V5"
        select = universe.select_atoms(frase)
        if str(select) != "<AtomGroup []>":
            juncao = str(line.resid) + str(line.resname)
            ids.append(line.resid)
            lista.append(juncao)
            dicti.update({'V5': lista})

print(dicti)

'''
Esta parte do código irá criar uma lista temporaria que sera utilizada no ficheiro python 
'Aromaticos.py'
'''
temp = []
j = -1
while j < len(ids):
    j += 1
    for i in ids:
        if i not in temp:
            temp.append(i)

