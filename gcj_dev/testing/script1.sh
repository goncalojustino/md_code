for f in protein_md{001..020}_A
do

cd /work3/$f
#house-keeping
rm \#* cluster*

#first things first - re-order water molecules by distance to metal+ion
echo "13 | 14" > maker #this is anion+metal, stored in 26
echo "q" >> maker
gmx  make_ndx -f md.gro -o index_analysis.ndx < maker
echo "26" > maker
echo "17" >> maker
gmx  trjorder -f md.xtc -s md.gro -n index_analysis.ndx -nshell nshell.xvg -o md_ordered.xtc -xvg xmgrace -r 0.5 -dt 1 < maker

#group together Protein + Metal + Ion
echo "1 | 13 | 14 " > maker # this is protein+anion+metal, stored in 27
echo "q" >> maker
gmx make_ndx -f md.gro -n index_analysis.ndx -o index_analysis.ndx < maker

#remove roto+trans; this selects protein for least squares fit and system for output
echo "1" > maker
echo "0" >> maker
gmx trjconv -f md_ordered.xtc -s md.tpr -o  md_fitted.xtc -fit rot+trans -n index_analysis.ndx -dt 1 < maker

#cluster by Protein+metal+ion and output full system clusters
echo "27" > maker
echo "0" >> maker
gmx cluster -f md_fitted.xtc -s md.tpr -method gromos -g cluster_out.log -dt 1 -cl clusters.pdb -n index_analysis.ndx -wcl 5 < maker

#split multimodel pdb file into various single model pdb files
csplit -k -s -n 3 -b '%02d.pdb' -f cluster_ clusters.pdb '/^ENDMDL/+1' {*}

#correct an atom nomenclature issue from PDB format
for g in cluster_00.pdb
 do
 gmx editconf -f $g -o temp.pdb
 sed "s/0.00           F/0.00          FE/" temp.pdb -i
 mv temp.pdb $g
done

#hbond analysis from htf_md, FIRST SHELL
#prepare index files for generic 13 | 14 (ion and metal) and for everything else (protein + water + ions)

#this selects protein+water+ions, stores as 28
printf "1 | 25 \n q \n" | gmx make_ndx -f md.gro -n index_analysis.ndx -o index_analysis.ndx
printf "26 \n28 \nq \n" | gmx hbond  -n index_analysis.ndx -f md_fitted.xtc -s md.tpr -g hbond.log -dist -ang -hx -hbn -hbm -dan -nhbdist

#rename 1st shell results
mv hbond.ndx hbond1.ndx
mv hbond.log hbond1.log
mv hbnum.xvg hbnum1.xvg
mv hbdist.xvg hbdist1.xvg
mv hbang.xvg hbang1.xvg
mv hbhelix.xvg hbhelix1.xvg
mv hbmap.xpm hbmap1.xpm
mv danum.xvg danum1.xvg

#hbond analysis from htf_md, SECOND SHELL
more hbond1.log | sed 1d | awk '{print $1; print $2; print $3}' | sed '/SOL/d' | sed '/338/d' | sed '/339/d' | cut -c-6  | cut -c 4- | uniq | awk '{print "r "$1" |"}'  | tr -d [A-Z] | sed ':a;N;$!ba;s/\n/ /g'|  rev | cut -c 2- | rev > maker
more hbond1.log | sed 1d | awk '{print $1; print $2; print $3}' | sed '/SOL/d' | sed '/338/d' | sed '/339/d' | cut -c-6  | cut -c 4- | uniq |  awk '{print " &! r "$1 }' | tr -d [A-Z] | sed ':a;N;$!ba;s/\n/ /g' | sed  's/^/1/' >> maker
echo "q" >> maker
#groups 29 and 30

printf "29 \n 30 \n q \n" | gmx hbond  -n index_analysis.ndx -f md_fitted.xtc -s md.tpr -g hbond.log -dist -ang -hx -hbn -hbm -dan -nhbdist -dt 1 

#rename 2nd shell results
mv hbond.ndx hbond2.ndx
mv hbond.log hbond2.log
mv hbnum.xvg hbnum2.xvg
mv hbdist.xvg hbdist2.xvg
mv hbang.xvg hbang2.xvg
mv hbhelix.xvg hbhelix2.xvg
mv hbmap.xpm hbmap2.xpm
mv danum.xvg danum2.xvg


#g_mmpbsa for Pr+An vs Met
# store Pr+An in 31
printf " 1 | 13 \n q \n" | gmx make_ndx -f md.gro -n index_analysis.ndx -o index_analysis.ndx 

#1 calculation of potential energy in vaccuum
printf "31 \n 14 \n" | g_mmpbsa -f md_fitted.xtc -s md.tpr -n index_analysis.ndx -pdie 2 -decomp -mm PAvM_energy_mm.xvg -mmcon PAvM_contrib_mm.dat -dt 10
#2 calculation of polar solvation energy
printf "31 \n 14 \n" | g_mmpbsa -f md_fitted.xtc -s md.tpr -n index_analysis.ndx -i /workopt/g_mmpbsa/parameter_files/mmpbsa_polar.mdp -nomme -pbsa -decomp -pol PAvM_polar.xvg -pcon PAvM_contrib_pol.dat -dt 10 < maker
#3 calculation of non-polar solvation energy --> WCA only model
printf "31 \n 14 \n" | g_mmpbsa -f md_fitted.xtc -s md.tpr -n index_analysis.ndx -i /workopt/g_mmpbsa/parameter_files/mmpbsa_wca.mdp -nomme -pbsa -decomp -apol PAvM_wca.xvg -apcon  PAvM_contrib_wca.dat -dt 10 < maker
#average free binding energy calculation
python3 /workopt/g_mmpbsa/scripts/MmPbSaStat.py -m PAvM_energy_mm.xvg  -p PAvM_polar.xvg -a PAvM_wca.xvg -bs -nbs 2000 -of PAvM_full_energy.dat -os PAvM_summary_energy.dat
# contribution of residues to the binding energy
python3 /workopt/g_mmpbsa/scripts/MmPbSaDecomp.py -bs -nbs 2000 -m PAvM_contrib_mm.dat -p PAvM_contrib_pol.dat -a PAvM_contrib_wca.dat -o PAvM_final_contrib_energy.dat -om PAvM_energyMapIn.dat

#g_mmpbsa for Pr vs An+Met
# stored in 26
#1 calculation of potential energy in vaccuum
printf "1 \n 26 \n" |  g_mmpbsa -f md_fitted.xtc -s md.tpr -n index_analysis.ndx -pdie 2 -decomp -mm PvAM_energy_mm.xvg -mmcon PvAM_contrib_mm.dat -dt 5 < maker
 #2 calculation of polar solvation energy
printf "1 \n 26 \n" | g_mmpbsa -f md_fitted.xtc -s md.tpr -n index_analysis.ndx -i /workopt/g_mmpbsa/parameter_files/mmpbsa_polar.mdp -nomme -pbsa -decomp -pol PvAM_polar.xvg -pcon PvAM_contrib_pol.dat -dt 5 < maker
 #3 calculation of non-polar solvation energy --> WCA only model
printf "1 \n 26 \n" | g_mmpbsa -f md_fitted.xtc -s md.tpr -n index_analysis.ndx -i /workopt/g_mmpbsa/parameter_files/mmpbsa_wca_ap.mdp -nomme -pbsa -decomp -apol PvAM_wca.xvg -apcon  PvAM_contrib_wca.dat -dt 5 < maker
 #average free binding energy calculation
 python3 /work_opt/g_mmpbsa/scripts/MmPbSaStat.py -m PvAM_energy_mm.xvg  -p PvAM_polar.xvg -a PvAM_wca.xvg -bs -nbs 2000 -of PvAM_full_energy.dat -os PvAM_summary_energy.dat
 # contribution of residues to the binding energy
 python3 /work_opt/g_mmpbsa/scripts/MmPbSaDecomp.py -bs -nbs 2000 -m PvAM_contrib_mm.dat -p PvAM_contrib_pol.dat -a PvAM_contrib_wca.dat -o PvAM_final_contrib_energy.dat -om PvAM_energyMapIn.dat

done
