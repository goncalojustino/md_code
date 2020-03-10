for f in protein_md{001..063}_A
do

cd /work6/$f

#hbond analysis from htf_md

#prepare index files for generic 13 | 14 (ion and metal) and for everything else (protein + water + ions)
echo " 13 | 14 " > maker
echo " 1 | 25 " >> maker
echo "q" >> maker
gmx make_ndx -f md.gro -o index_hbond.ndx < maker
#run bond analysis
echo " 26 " > maker
echo " 27 " >> maker
gmx hbond  -n index_hbond.ndx -f md_fitted.xtc -s md.tpr -dt 1 -g hbond.log -dist -ang -hx -hbn -hbm -dan -nhbdist < maker

done
