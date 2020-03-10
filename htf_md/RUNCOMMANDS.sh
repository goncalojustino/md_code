md001_[A-D] is doing
check that the differente cases show at least one coherent coordination

must change ggroups in run_me for others






cp ../mdp_templates/*mdp .
grep HET protein.pdb  | head -n 10
grep FE3 *mdp
sed "s/COA/NO3/" *mdp -i
sed "s/FE3/V3 /" *mdp -i


####
for f in md0[0-9][0-9]_[A-D]; do cd $f; bash run_me; cd ../; done

for f in md0[0-9][0-9]_[A-D]; do cp -fv md001_A/run_me $f; done

####


########################################
#clean up previous attemps
rm \#*

gmx-mod pdb2gmx -water spc -ff zgj6 -o processed.gro -p topol.top -f protein.pdb

#create index input file and create the actual index file
echo " \"Protein\" | \"COA\" | \"FE3\" " > index_creator_em
echo "q" >> index_creator_em
gmx-mod make_ndx -f processed.gro -o index_master.ndx < index_creator_em

gmx-mod editconf -f processed.gro -o box.gro -bt cubic -d 1.6
gmx-mod solvate -cp box.gro -cs spc216.gro -p topol.top -o solv.gro
gmx-mod grompp -f em.mdp -c solv.gro -p topol.top -o ions.tpr

#this replaces group 16 "SOL" molecules with ions
#CHECK THAT GROUP 16 IS SOL
echo 16 | gmx-mod genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15

gmx-mod grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr -n index_master.ndx

gmx-mod mdrun -v -deffnm em -nb gpu -ntomp 14 -pin on

#create index for dynamics with protein+prosthetics and with water+ions 
#for use in nvt, npt and md steps
echo " 1 | 13 | 14 " > index_creator_md
echo "\"Water_and_ions\"" >> index_creator_md
echo "q" >> index_creator_md
gmx-mod make_ndx -f em.gro -o index_md.ndx < index_creator_md




gmx-mod grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr -n index_md.ndx
gmx-mod mdrun -v -deffnm nvt -pin on -ntomp 14 -nb gpu
gmx-mod grompp -f npt.mdp -c nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -n index_md.ndx
gmx-mod mdrun -v -deffnm npt -pin on -ntomp 14 -nb gpu
gmx-mod grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr -n index_md.ndx
gmx-mod mdrun -v -deffnm md -nb gpu -pin on -ntomp 8



gmx editconf -f processed.gro -o box.gro -bt cubic -d 1.6
gmx solvate -cp box.gro -cs spc216.gro -p topol.top -o solv.gro
gmx grompp -f em.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 1
gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.154
gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em -nb gpu -nt 12 -pin on
sed "s/SDMSO/SDmso/" *itp -i
gmx grompp -f nvt.mdp -c em.gro -r rm.gro -p topol.top -o nvt.tpr
gmx mdrun -v -deffnm nvt -nb gpu -nt 12 -pin on
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -v -deffnm npt -nb gpu -nt 12 -pin on
gmx grompp -f md2.mdp -c npt.gro -r npt.gro -t npt.cpt -p topol.top -o md2.tpr
gmx mdrun -v -deffnm md2 -nb gpu -nt 12 -pin on


