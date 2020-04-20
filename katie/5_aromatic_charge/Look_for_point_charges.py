#IMPORTS
import MDAnalysis as mda

universe = mda.Universe("./md003/md3.gro", "./md003/md3.xtc")


interactions = universe.select_atoms("(resname ARG and name NH1) or (resname ASP and name OD*) "
                                     "or (resname GLU and name OE*) or (resname HIS and name ND1)"
                                     "or (resname HIT and name NE2) or (resname LYS* and name NZ)"
                                     "or (resname TYO and name OH)")

iron = universe.select_atoms("(resname FE3 and name FE) or (resname V3 and name V3)"
                              "or (resname VO and name V4) or (resname VO2 and name V5)"
                              "or (resname VO4 and name V5)")


my_list = []
ids = []
dicti = {}


"""
This function will create a dictionary that will contain as key all the atoms that are considered
essential, as they are characteristic due to their load, and as key the name and id of all
amino acids that meet the requirements, that is, it has a loaded atom in its composition
"""
# Selection of NH1 from ARG
for line in interactions:
    if 'NH1' in str(line):
        phrase = "resid " + str(line.resid) + " and (name HH11 or name HH12)"
        select = universe.select_atoms(phrase)
        if str(select) != "<AtomGroup []>":
            junction = str(line.resid) + str(line.resname)
            ids.append(line.resid)
            my_list.append(junction)
            dicti.update({'NH1': my_list})

# Selection of OD* from ASP
my_list = []
for line in interactions:
    if 'ASP' in str(line):
        if line.resid not in my_list:
            junction = str(line.resid) + str(line.resname)
            ids.append(line.resid)
            my_list.append(junction)
            dicti.update({'OD*': my_list})

# Selection of  OE* from GLU
my_list = []
for line in interactions:
    if 'GLU' in str(line):
        if line.resid not in my_list:
            junction = str(line.resid) + str(line.resname)
            ids.append(line.resid)
            my_list.append(junction)
            dicti.update({'OE*': my_list})

# Selection of  ND1 from HIS
my_list = []
for line in interactions:
    if 'ND1' in str(line):
        phrase = "resid " + str(line.resid) + " and name HD1"
        select = universe.select_atoms(phrase)
        for line1 in range(len(select)):
            phrase1 = "resid " + str(select[line1].resid) + " and name HE2"
            select1 = universe.select_atoms(phrase1)
            if str(select1) != "<AtomGroup []>":
                junction = str(select1[line1].resid) + str(select1[line1].resname)
                ids.append(line.resid)
                my_list.append(junction)
                dicti.update({'ND1': my_list})

# Selection of  NE2 from HIT
my_list = []
for line in interactions:
    if 'NE2' in str(line):
        if line.resid not in my_list:
            junction = str(line.resid) + str(line.resname)
            ids.append(line.resid)
            my_list.append(junction)
            dicti.update({'NE2': my_list})

# Selection of  HZ3 from LYS+
my_list = []
lys = ['LYS', 'LYSX', 'LYSY']
for pos in lys:
    for line in interactions:
        if pos in str(line):
            phrase = "resid " + str(line.resid) + " and name HZ3"
            select = universe.select_atoms(phrase)
            if str(select) != "<AtomGroup []>":
                junction = str(line.resid) + str(line.resname)
                ids.append(line.resid)
                my_list.append(junction)
                dicti.update({'HZ3': my_list[:-1]})

# Selection of  OH from TYO
my_list = []
for line in interactions:
    if 'OH' in str(line):
        phrase = "resid " + str(line.resid) + " and name OH"
        select = universe.select_atoms(phrase)
        if str(select) != "<AtomGroup []>":
            junction = str(line.resid) + str(line.resname)
            ids.append(line.resid)
            my_list.append(junction)
            dicti.update({'OH': my_list})

# Selection of FE from FE3
my_list = []
for line in iron:
    if 'FE3' in str(line):
        phrase = "resid " + str(line.resid) + " and name FE"
        select = universe.select_atoms(phrase)
        if str(select) != "<AtomGroup []>":
            junction = str(line.resid) + str(line.resname)
            ids.append(line.resid)
            my_list.append(junction)
            dicti.update({'FE': my_list})

# Selection of V3 from V3
my_list = []
for line in iron:
    if 'V3' in str(line):
        phrase = "resid " + str(line.resid) + " and name V3"
        select = universe.select_atoms(phrase)
        if str(select) != "<AtomGroup []>":
            junction = str(line.resid) + str(line.resname)
            ids.append(line.resid)
            my_list.append(junction)
            dicti.update({'V3': my_list})

# Selection of V4 from VO
my_list = []
for line in iron:
    if 'VO' in str(line):
        phrase = "resid " + str(line.resid) + " and name V4"
        select = universe.select_atoms(phrase)
        if str(select) != "<AtomGroup []>":
            junction = str(line.resid) + str(line.resname)
            ids.append(line.resid)
            my_list.append(junction)
            dicti.update({'V4': my_list})

# Selection of V5 from VO2
my_list = []
for line in iron:
    if 'VO2' in str(line):
        phrase = "resid " + str(line.resid) + " and name V5"
        select = universe.select_atoms(phrase)
        if str(select) != "<AtomGroup []>":
            junction = str(line.resid) + str(line.resname)
            ids.append(line.resid)
            my_list.append(junction)
            dicti.update({'V5': my_list})

# Selection of V5 from VO4
my_list = []
for line in iron:
    if 'VO4' in str(line):
        phrase = "resid " + str(line.resid) + " and name V5"
        select = universe.select_atoms(phrase)
        if str(select) != "<AtomGroup []>":
            junction = str(line.resid) + str(line.resname)
            ids.append(line.resid)
            my_list.append(junction)
            dicti.update({'V5': my_list})


"""
This code block will create a temporary list containing all the amino acid id's used in the 
dictionary and will be used in the python file 'Aromaticos.py'
"""
temp = []
j = -1
while j < len(ids):
    j += 1
    for i in ids:
        if i not in temp:
            temp.append(i)
