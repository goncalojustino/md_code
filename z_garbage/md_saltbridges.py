import numpy, MDAnalysis
from MDAnalysis.analysis import contacts
import os
os.chdir(/work3/protein_md001_A_TEMP/)

universe = MDAnalysis.Universe('md.gro','md_fitted.xtc')
acid = universe.select_atoms("(resname ASP and name OD*) or (resname GLU and name OE*) or (resname TYO and name OH) or (resname HIT and name ND*) or (resname HIT and name NE*)")
basic = universe.select_atoms("(resname ARG and name NH*) or (resname LYS* and name NZ*) or (resname HISH and name NE*) or (resname HISH and name ND*)")

distance = MDAnalysis.analysis.distances.distance_array(acid.positions,basic.positions)

numpy.savetxt('saltbridges_distance_matrix.dat',distance,delimiter=" ")
with open('saltbridges_acid_list.out','w') as f:
    print(list(acid),file=f)

with open('saltbridges_basic_list','w') as f:
    print(list(basic),file=f)


