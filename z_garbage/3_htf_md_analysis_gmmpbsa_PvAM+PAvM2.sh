!#\bin\bash

# g_mmpbsa analysis on protein_md runs for htf @marycooper

for f in protein_md{050..056}_A    # CORRECT FOLDERS
 do
  cd /work6/$f              # CORRECT LOCATION
  rm \#*                    # house cleaning never hurt a soul

 #PROTEIN+ANION vs. METAL
 echo "1 | 13" > maker # this is protein + anion
 echo "q" >> maker
 gmx make_ndx -f md.gro -o index_mmpbsa.ndx < maker

 #this handles protein+anion vs. metal
 echo "14" > maker
 echo "26" >> maker
 #1 calculation of potential energy in vaccuum
 g_mmpbsa -f md_fitted.xtc -s md.tpr -n index_mmpbsa.ndx -pdie 2 -decomp -mm PAvM_energy_mm.xvg -mmcon PAvM_contrib_mm.dat -dt 1000 < maker
 #2 calculation of polar solvation energy
 g_mmpbsa -f md_fitted.xtc -s md.tpr -n index_mmpbsa.ndx -i /work_opt/g_mmpbsa/mdp_samples/mmpbsa_polar.mdp -nomme -pbsa -decomp -pol PAvM_polar.xvg -pcon PAvM_contrib_pol.dat -dt 1000 < maker
 #3 calculation of non-polar solvation energy --> WCA only model
 g_mmpbsa -f md_fitted.xtc -s md.tpr -n index_mmpbsa.ndx -i /work_opt/g_mmpbsa/mdp_samples/mmpbsa_wca_ap.mdp -nomme -pbsa -decomp -apol PAvM_wca.xvg -apcon  PAvM_contrib_wca.dat -dt 1000 < maker
 #average free binding energy calculation
 python3 /work_opt/g_mmpbsa/scripts/MmPbSaStat.py -m PAvM_energy_mm.xvg  -p PAvM_polar.xvg -a PAvM_wca.xvg -bs -nbs 2000 -of PAvM_full_energy.dat -os PAvM_summary_energy.dat
 # contribution of residues to the binding energy
 python3 /work_opt/g_mmpbsa/scripts/MmPbSaDecomp.py -bs -nbs 2000 -m PAvM_contrib_mm.dat -p PAvM_contrib_pol.dat -a PAvM_contrib_wca.dat -o PAvM_final_contrib_energy.dat -om PAvM_energyMapIn.dat

 
 #this handles protein vs. anion+metal
 echo "26" > maker
 echo "1" >> maker
 #1 calculation of potential energy in vaccuum
 g_mmpbsa -f md_fitted.xtc -s md.tpr -n index_mmpbsa.ndx -pdie 2 -decomp -mm PvAM_energy_mm.xvg -mmcon PvAM_contrib_mm.dat -dt 20 < maker
 #2 calculation of polar solvation energy
 g_mmpbsa -f md_fitted.xtc -s md.tpr -n index_mmpbsa.ndx -i /work_opt/g_mmpbsa/mdp_samples/mmpbsa_polar.mdp -nomme -pbsa -decomp -pol PvAM_polar.xvg -pcon PvAM_contrib_pol.dat -dt 20 < maker
 #3 calculation of non-polar solvation energy --> WCA only model
 g_mmpbsa -f md_fitted.xtc -s md.tpr -n index_mmpbsa.ndx -i /work_opt/g_mmpbsa/mdp_samples/mmpbsa_wca_ap.mdp -nomme -pbsa -decomp -apol PvAM_wca.xvg -apcon  PvAM_contrib_wca.dat -dt 20 < maker
 #average free binding energy calculation
 python3 /work_opt/g_mmpbsa/scripts/MmPbSaStat.py -m PvAM_energy_mm.xvg  -p PvAM_polar.xvg -a PvAM_wca.xvg -bs -nbs 2000 -of PvAM_full_energy.dat -os PvAM_summary_energy.dat
 # contribution of residues to the binding energy
 python3 /work_opt/g_mmpbsa/scripts/MmPbSaDecomp.py -bs -nbs 2000 -m PvAM_contrib_mm.dat -p PvAM_contrib_pol.dat -a PvAM_contrib_wca.dat -o PvAM_final_contrib_energy.dat -om PvAM_energyMapIn.dat


 cd ..
 done
