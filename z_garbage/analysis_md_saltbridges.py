import numpy, MDAnalysis
from MDAnalysis.analysis import contacts

universe = MDAnalysis.Universe('md2.gro','md2_fitted.xtc')
acid = universe.select_atoms("(resname ASP and name OD*) or (resname GLU and name OE*) or (resname TYO and name OH) or (resname HIT and name ND*) or (resname HIT and name NE*)")
basic = universe.select_atoms("(resname ARG and name NH*) or (resname LYS* and name NZ*) or (resname HISH and name NE*) or (resname HISH and name ND*)")

distance = MDAnalysis.analysis.distances.distance_array(acid.positions,basic.positions)

numpy.savetxt('Saltbridges_distance_matrix.dat',distance,delimiter=" ")

with open('Saltbridges_acid_list.dat','w') as f:
    print(list(acid),file=f)

with open('Saltbridges_basic_list.dat','w') as f:
    print(list(basic),file=f)


