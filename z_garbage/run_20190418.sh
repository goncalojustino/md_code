*** we are running on /work7

*** CHANGE rcoulomb rvdw rlist from 1.0 to 1.4
*** CHECK groups on mdp files ALL

gmx pdb2gmx -water spce -ff zgj12 -o processed.gro -p topol2.top -f protein.pdb

gmx editconf -f processed.gro -o box.gro -bt cubic -d 1.6
gmx solvate -cp box.gro -cs spc216.gro -p topol2.top -o solv.gro
gmx grompp -f em.mdp -c solv.gro -p topol2.top -o ions.tpr -maxwarn 1
echo "16" | gmx genion -s ions.tpr -o solv_ions.gro -p topol2.top -pname NA -nname CL -neutral -conc 0.154
gmx grompp -f em.mdp -c solv_ions.gro -p topol2.top -o em.tpr
gmx mdrun -v -deffnm em -nb auto  -nt 12 -pin on

*** why is this here ? sed "s/SDMSO/SDmso/" *itp -i


printf "1 | 13 | 14 \n q \n" | gmx make_ndx -f em.gro  -o index_md.ndx

gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol2.top -o nvt.tpr -n index_md.ndx
gmx mdrun -v -deffnm nvt -nb auto -nt 12 -pin on
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol2.top -o npt.tpr -n index_md.ndx
gmx mdrun -v -deffnm npt -nb auto -nt 12 -pin on

*** change times in md2.mdp 

gmx grompp -f md2.mdp -c npt.gro -r npt.gro -t npt.cpt -p topol2.top -o md2.tpr -n index_md.ndx
gmx mdrun -v -deffnm md2 -nb auto -nt 12 -pin on
done