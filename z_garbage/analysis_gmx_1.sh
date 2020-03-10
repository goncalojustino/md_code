#processing of hTf + metal + ion MD runs

source /workopt/gromacs/bin/GMXRC.bash

#CLUSTER ANALYSIS START===========================================================================================
#1.RE-ORDER water molecules by distance to metal+ion
#1a. create anion+metal group, stored in 26
printf "13 | 14 \n q \n" | gmx  make_ndx -f md2.gro -o index_analysis.ndx
#1b. order water [17] by distance to anion+metal [26]
printf " 26 \n 17 \n" | gmx  trjorder -f md2.xtc -s md2.gro -n index_analysis.ndx -nshell nshell.xvg -o md_ordered.xtc -xvg xmgrace -r 0.5 -dt 1 
#1c. create protein+metal+ion, stored in 27
printf "1 | 13 | 14 \n q \n" | gmx make_ndx -f md2.gro -n index_analysis.ndx -o index_analysis.ndx

#1d. remove rotational and translational movement
#analyse protein [1] and output system [0]
printf "1\n 0\n" | gmx trjconv -f md_ordered.xtc -s md2.tpr -o  md_fitted.xtc -fit rot+trans -n index_analysis.ndx -dt 1

#1e. cluster by protein+metal+ion[27] and output system [0] in clusters
printf "27 \n 0\n " | gmx cluster -f md_fitted.xtc -s md2.tpr -method gromos -g cluster_out.log -dt 1 -cl clusters.pdb -n index_analysis.ndx -wcl 5

#1f. split multimodel pdb file into various single model pdb files
csplit -k -s -n 3 -b '%02d.pdb' -f cluster_ clusters.pdb '/^ENDMDL/+1' {*}
#correct an FE atom nomenclature issue from PDB format; only for cluster_00.pdb
for g in cluster_00.pdb
 do
 gmx editconf -f $g -o temp.pdb
 sed "s/0.00           F/0.00          FE/" temp.pdb -i
 mv temp.pdb $g
done
#CLUSTER ANALYSIS END===========================================================================================

#H BONDING ANALYSIS START===========================================================================================
#2. FIRST SHELL ANALYSIS of H bonding contacts
#2a. group metal+ion is alreadt created [26]
#2b. group everything else (protein+water+ions) and store to 28
printf "1 | 25 \n q \n" | gmx make_ndx -f md2.gro -n index_analysis.ndx -o index_analysis.ndx
#2c. analyse hbonding between ANION + METAL (26) and PROTEIN + WATER + IONS (28)
printf "26 \n28 \nq \n" | gmx hbond  -n index_analysis.ndx -f md_fitted.xtc -s md2.tpr -g hbond.log -dist -ang -hx -hbn -hbm -dan -nhbdist -dt 1
#2d rename 1st shell results
for f in hb*; do rename 's/\./\_1./' $f; done
#3. SECOND SHELL ANALYSIS of H bonding contacts
#3a. get the residues that are in the 1st shell into [29] and the others into [30]
more hbond_1.log | sed 1d | awk '{print $1; print $2; print $3}' | sed '/SOL/d' | sed '/338/d' | sed '/339/d' | cut -c-6  | cut -c 4- | uniq | awk '{print "r "$1" |"}'  | tr -d [A-Z] | sed ':a;N;$!ba;s/\n/ /g'|  rev | cut -c 2- | rev > maker
more hbond_1.log | sed 1d | awk '{print $1; print $2; print $3}' | sed '/SOL/d' | sed '/338/d' | sed '/339/d' | cut -c-6  | cut -c 4- | uniq |  awk '{print " &! r "$1 }' | tr -d [A-Z] | sed ':a;N;$!ba;s/\n/ /g' | sed  's/^/1/' >> maker
echo "name 29 group29" >> maker
echo "name 30 group30" >> maker
echo "q" >> maker
gmx make_ndx -f md2.gro -n index_analysis.ndx -o index_analysis.ndx < maker
#3b. Analyse H bonding from group 29 (1st shell) to the remainder of the system
printf "29 \n 30 \n q \n" | gmx hbond  -n index_analysis.ndx -f md_fitted.xtc -s md2.tpr -g hbond.log -dist -ang -hx -hbn -hbm -dan -nhbdist -dt 1
#3c. rename 2nd shell results
for f in hb*[a-z].*; do rename 's/\./\_2./' $f; done
#H BONDING ANALYSIS END===========================================================================================

# G_CONTACTS START==========================================================================================
# use 4.57
source /workopt/gromacs457/bin/GMXRC.bash
#4a. ANION + METAL group is present as [26]
#4b. analyse contacts from anion+metal [26] to system (water!) - no slowdown
printf " 26 \n 0 \n" | g_contacts -f md_ordered.xtc -n index_analysis.ndx -dt 1 -resndx -s md2.gro -d 0.3 -on contacts_res_1.ndx -o contacts_1.dat
#4c. get contacts from 1st sphere into [31] and eveything else into [32]
more contacts_1.dat | sed 1d | sed '/SOL/d' | awk '{print $7}' | sort | uniq  | awk '{print "r "$1" |"}' | sed ':a;N;$!ba;s/\n/ /g'  |  rev | cut -c 2- | rev > maker
more contacts_1.dat | sed 1d | sed '/SOL/d' | awk '{print $7}' | sort | uniq  | awk '{print " &! r "$1 }' | sed ':a;N;$!ba;s/\n/ /g' | sed  's/^/1/' >> maker
echo "name 31 group31" >> maker
echo "name 32 group32" >> maker
echo q >> maker
make_ndx -f md2.gro -n index_analysis.ndx -o index_analysis.ndx < maker
#4d. analyse contacts
printf "31 | 32" | g_contacts -f md_ordered.xtc -n index_analysis.ndx -dt 1 -resndx -s md2.gro -d 0.35 -on contacts_res_2.ndx -o contacts_2.dat 
# G_CONTACTS END==========================================================================================



CHANGE BACK TO NEW GROMACS
# G_MMPBSA for PROT+ANION VS MET==================================================================================
#4a. store Pr+An in 31
printf " 1 | 13 \n q \n" | gmx make_ndx -f md.gro -n index_analysis.ndx -o index_analysis.ndx 
#4b. calculation of potential energy in vaccuum
printf "31 \n 22 \n" | g_mmpbsa -f md_fitted.xtc -s md.tpr -n index_analysis.ndx -pdie 2 -decomp -mm PAvM_energy_mm.xvg -mmcon PAvM_contrib_mm.dat -dt 10
#4c. calculation of polar solvation energy
printf "31 \n 22 \n" | g_mmpbsa -f md_fitted.xtc -s md.tpr -n index_analysis.ndx -i /workopt/g_mmpbsa/parameter_files/mmpbsa_polar.mdp -nomme -pbsa -decomp -pol PAvM_polar.xvg -pcon PAvM_contrib_pol.dat -dt 10
#4d. calculation of non-polar solvation energy --> WCA only model
printf "31 \n 22 \n" | g_mmpbsa -f md_fitted.xtc -s md.tpr -n index_analysis.ndx -i /workopt/g_mmpbsa/parameter_files/mmpbsa_wca.mdp -nomme -pbsa -decomp -apol PAvM_wca.xvg -apcon  PAvM_contrib_wca.dat -dt 10
#4e. average free binding energy calculation
python3 /workopt/g_mmpbsa/scripts/MmPbSaStat.py -m PAvM_energy_mm.xvg  -p PAvM_polar.xvg -a PAvM_wca.xvg -bs -nbs 2000 -of PAvM_full_energy.dat -os PAvM_summary_energy.dat
#4f. contribution of residues to the binding energy
python3 /workopt/g_mmpbsa/scripts/MmPbSaDecomp.py -bs -nbs 2000 -m PAvM_contrib_mm.dat -p PAvM_contrib_pol.dat -a PAvM_contrib_wca.dat -o PAvM_final_contrib_energy.dat -om PAvM_energyMapIn.dat

# G_MMPBSA for PROT VS ANION+MET==================================================================================
# ANION + METAL is already defined in group 26
#5a. calculation of potential energy in vaccuum
printf "1 \n 26 \n" | g_mmpbsa -f md_fitted.xtc -s md.tpr -n index_analysis.ndx -pdie 2 -decomp -mm PAvM_energy_mm.xvg -mmcon PAvM_contrib_mm.dat -dt 10
#5b. calculation of polar solvation energy
printf "1 \n 26 \n" | g_mmpbsa -f md_fitted.xtc -s md.tpr -n index_analysis.ndx -i /workopt/g_mmpbsa/parameter_files/mmpbsa_polar.mdp -nomme -pbsa -decomp -pol PAvM_polar.xvg -pcon PAvM_contrib_pol.dat -dt 10
#5c. calculation of non-polar solvation energy --> WCA only model
printf "1 \n 26 \n" | g_mmpbsa -f md_fitted.xtc -s md.tpr -n index_analysis.ndx -i /workopt/g_mmpbsa/parameter_files/mmpbsa_wca.mdp -nomme -pbsa -decomp -apol PAvM_wca.xvg -apcon  PAvM_contrib_wca.dat -dt 10
#5d. average free binding energy calculation
python3 /workopt/g_mmpbsa/scripts/MmPbSaStat.py -m PAvM_energy_mm.xvg  -p PAvM_polar.xvg -a PAvM_wca.xvg -bs -nbs 2000 -of PAvM_full_energy.dat -os PAvM_summary_energy.dat
#5e.  contribution of residues to the binding energy
python3 /workopt/g_mmpbsa/scripts/MmPbSaDecomp.py -bs -nbs 2000 -m PAvM_contrib_mm.dat -p PAvM_contrib_pol.dat -a PAvM_contrib_wca.dat -o PAvM_final_contrib_energy.dat -om PAvM_energyMapIn.dat

#g_contacts must use older gromacs


source /workopt/gromacs/bin/GMXRC.bash

