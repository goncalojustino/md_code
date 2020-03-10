
source /workopt/gromacs/bin/GMXRC.bash


#first things first - re-order water molecules by distance to metal+ion
#create ANION + METAL, stored in 26
printf "13 | 14 \n q \n" | gmx  make_ndx -f md2.gro -o index_analysis.ndx
#order water (17) by distance to 26
printf " 26 \n 17 \n" | gmx  trjorder -f md2.xtc -s md2.gro -n index_analysis.ndx -nshell nshell.xvg -o md2_ordered.xtc -xvg xmgrace -r 0.5 -dt 1 
#create PROTEIN + METAL + ION, stored in 27
printf "1 | 13 | 14 \n q \n" | gmx make_ndx -f md2.gro -n index_analysis.ndx -o index_analysis.ndx

#remove roto+trans; analyse protein (1) and output system (0)
printf "1\n 0\n" | gmx trjconv -f md2_ordered.xtc -s md2.tpr -o  md2_fitted.xtc -fit rot+trans -n index_analysis.ndx -dt 1

#cluster by PROTEIN + METAL + ION (27) and output system (0) clusters
printf "27 \n 0\n " | gmx cluster -f md2_fitted.xtc -s md2.tpr -method gromos -g cluster_out.log -dt 10 -cl clusters.pdb -n index_analysis.ndx -wcl 5

#split multimodel pdb file into various single model pdb files
csplit -k -s -n 3 -b '%02d.pdb' -f cluster_ clusters.pdb '/^ENDMDL/+1' {*}
#correct an atom nomenclature issue from PDB format; only for cluster_00.pdb
for g in cluster_00.pdb
 do
 gmx editconf -f $g -o temp.pdb
 sed "s/0.00           F/0.00          FE/" temp.pdb -i
 mv temp.pdb $g
done

#hbond analysis from htf_md, FIRST SHELL
#prepare index files for generic 13 | 14 (ion and metal) and for everything else (protein + water + ions)

#this selects PROTEIN + WATER + IONS, stores as 28
printf "1 | 25 \n q \n" | gmx make_ndx -f md2.gro -n index_analysis.ndx -o index_analysis.ndx
#analyse hbonding between ANION + METAL (26) and PROTEIN + WATER + IONS (28)
printf "26 \n28 \nq \n" | gmx hbond  -n index_analysis.ndx -f md2_fitted.xtc -s md2.tpr -g hbond.log -dist -ang -hx -hbn -hbm -dan -nhbdist
#rename 1st shell results
for f in hb*; do rename 's/\./\_1./' $f; done

#hbond analysis from htf_md, SECOND SHELL
more hbond_1.log | sed 1d | awk '{print $1; print $2; print $3}' | sed '/SOL/d' | sed '/338/d' | sed '/339/d' | cut -c-6  | cut -c 4- | uniq | awk '{print "r "$1" |"}'  | tr -d [A-Z] | sed ':a;N;$!ba;s/\n/ /g'|  rev | cut -c 2- | rev > maker
more hbond_1.log | sed 1d | awk '{print $1; print $2; print $3}' | sed '/SOL/d' | sed '/338/d' | sed '/339/d' | cut -c-6  | cut -c 4- | uniq |  awk '{print " &! r "$1 }' | tr -d [A-Z] | sed ':a;N;$!ba;s/\n/ /g' | sed  's/^/1/' >> maker
echo "name 29 group29" >> maker
echo "name 30 group30" >> maker
echo "q" >> maker
#groups 29 and 30
gmx make_ndx -f md2.gro -n index_analysis.ndx -o index_analysis.ndx < maker
printf "29 \n 30 \n q \n" | gmx hbond  -n index_analysis.ndx -f md2_fitted.xtc -s md2.tpr -g hbond.log -dist -ang -hx -hbn -hbm -dan -nhbdist -dt 1 
#rename 2nd shell results
for f in hb*; do rename 's/\./\_2./' $f; done


#g_contacts must use older gromacs
source /workopt/gromacs457/bin/GMXRC.bash

# ANION + METAL group, stored as 34
printf " 17 | 18 \n q \n" | make_ndx -f md2.gro -o index_contacts.ndx
#analyse contacts from group 34 to system (water!) - no slowdown
printf " 34 \n 0 \n" | g_contacts -f md2_ordered.xtc -n index_contacts.ndx -dt 1 -resndx -s md2.gro -d 0.3 -on contacts_res_1.ndx -o contacts_1.dat

#get contacts to second coordination 
more contacts_1.dat | sed 1d | sed '/SOL/d' | awk '{print $7}' | sort | uniq  | awk '{print "r "$1" |"}' | sed ':a;N;$!ba;s/\n/ /g'  |  rev | cut -c 2- | rev > maker
more contacts_1.dat | sed 1d | sed '/SOL/d' | awk '{print $7}' | sort | uniq  | awk '{print " &! r "$1 }' | sed ':a;N;$!ba;s/\n/ /g' | sed  's/^/1/' >> maker
echo "name 35 group35" >> maker
echo "name 36 group36" >> maker
echo q >> maker
make_ndx -f md2.gro -n index_contacts.ndx -o index_contacts.ndx < maker

printf "35 \n 36\n" | g_contacts -f md2_ordered.xtc -n index_contacts.ndx -dt 1 -resndx -s md2.gro -d 0.35 -on contacts_res_2.ndx -o contacts_2.dat 

source /workopt/gromacs/bin/GMXRC.bash

