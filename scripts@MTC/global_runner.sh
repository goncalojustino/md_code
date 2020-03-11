for f in protein_md0{02..28}_A protein_md0{43..56}_A
#for f in protein_md001_A
do 

cd /work9/$f
echo $f
#/workcode/analysis/analysis_gmx_script2_complete.sh
#/workcode/analysis/analysis_gmx_script3_noMMPBSA.sh
#/workcode/analysis/analysis_gmx_script3_onlyHbond.sh
/workcode/analysis/analysis_gmx_script4_onlyMMPBSA.sh

#the following calls conda md env and requires md.gro amd md_fitted.xtc in each folder
#/workopt/miniconda3/envs/md/bin/ipython "-c %run /workcode/analysis/analysis_md_saltbridges.py"
#/workopt/miniconda3/envs/md/bin/ipython "-c %run /workcode/analysis/analysis_md_pipi.py"
#/workopt/miniconda3/envs/md/bin/ipython "-c %run /workcode/analysis/analysis_md_cationpi.py"
#/workopt/miniconda3/envs/md/bin/ipython "-c %run /workcode/analysis/analysis_md_anionpi.py"
#/workopt/miniconda3/envs/md/bin/ipython "-c %run /workcode/analysis/analysis_md_argarg.py"
done
