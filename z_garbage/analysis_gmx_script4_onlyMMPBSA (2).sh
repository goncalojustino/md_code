source /workopt/gromacs/bin/GMXRC.bash


# store Pr+An in 26
printf " 1 | 13 \n q \n" | gmx make_ndx -f md2.gro -o index_mmpbsa.ndx 
# store An+Met in 27
printf " 13 | 14 \n q \n" | gmx make_ndx -f md2.gro -n index_mmpbsa.ndx -o index_mmpbsa.ndx 

#g_mmpbsa for Pr+An vs Met
#1 calculation of potential energy in vaccuum
printf "26 \n 22 \n" | g_mmpbsa -f md2_fitted.xtc -s md2.tpr -n index_mmpbsa.ndx -pdie 2 -decomp -mm PAvM_energy_mm.xvg -mmcon PAvM_contrib_mm.dat -dt 1000
#2 calculation of polar solvation energy
printf "26 \n 22 \n" | g_mmpbsa -f md2_fitted.xtc -s md2.tpr -n index_mmpbsa.ndx -i /workopt/g_mmpbsa/parameter_files/mmpbsa_polar.mdp -nomme -pbsa -decomp -pol PAvM_polar.xvg -pcon PAvM_contrib_pol.dat -dt 1000
#3 calculation of non-polar solvation energy --> WCA only model
printf "26 \n 22 \n" | g_mmpbsa -f md2_fitted.xtc -s md2.tpr -n index_mmpbsa.ndx -i /workopt/g_mmpbsa/parameter_files/mmpbsa_wca.mdp -nomme -pbsa -decomp -apol PAvM_wca.xvg -apcon  PAvM_contrib_wca.dat -dt 1000
#average free binding energy calculation
python3 /workopt/g_mmpbsa/scripts/MmPbSaStat.py -m PAvM_energy_mm.xvg  -p PAvM_polar.xvg -a PAvM_wca.xvg -bs -nbs 2000 -of PAvM_full_energy.dat -os PAvM_summary_energy.dat
# contribution of residues to the binding energyls 
python3 /workopt/g_mmpbsa/scripts/MmPbSaDecomp.py -bs -nbs 2000 -m PAvM_contrib_mm.dat -p PAvM_contrib_pol.dat -a PAvM_contrib_wca.dat -o PAvM_final_contrib_energy.dat -om PAvM_energyMapIn.dat

#g_mmpbsa for Pr vs An+Met
# ANION + METAL is already defined in group 27
#1 calculation of potential energy in vaccuum
printf "1 \n 27 \n" | g_mmpbsa -f md2_fitted.xtc -s md2.tpr -n index_mmpbsa.ndx -pdie 2 -decomp -mm PvAM_energy_mm.xvg -mmcon PvAM_contrib_mm.dat -dt 1000
#2 calculation of polar solvation energy
printf "1 \n 27 \n" | g_mmpbsa -f md2_fitted.xtc -s md2.tpr -n index_mmpbsa.ndx -i /workopt/g_mmpbsa/parameter_files/mmpbsa_polar.mdp -nomme -pbsa -decomp -pol PvAM_polar.xvg -pcon PvAM_contrib_pol.dat -dt 1000
#3 calculation of non-polar solvation energy --> WCA only model
printf "1 \n 27 \n" | g_mmpbsa -f md2_fitted.xtc -s md2.tpr -n index_mmpbsa.ndx -i /workopt/g_mmpbsa/parameter_files/mmpbsa_wca.mdp -nomme -pbsa -decomp -apol PvAM_wca.xvg -apcon  PvAM_contrib_wca.dat -dt 1000
#average free binding energy calculation
python3 /workopt/g_mmpbsa/scripts/MmPbSaStat.py -m PvAM_energy_mm.xvg  -p PvAM_polar.xvg -a PvAM_wca.xvg -bs -nbs 2000 -of PvAM_full_energy.dat -os PvAM_summary_energy.dat
# contribution of residues to the binding energy
python3 /workopt/g_mmpbsa/scripts/MmPbSaDecomp.py -bs -nbs 2000 -m PvAM_contrib_mm.dat -p PvAM_contrib_pol.dat -a PvAM_contrib_wca.dat -o PvAM_final_contrib_energy.dat -om PvAM_energyMapIn.dat


