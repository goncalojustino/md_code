
source /workopt/gromacs/bin/GMXRC.bash


rm *hb*
#hbond analysis from htf_md, FIRST SHELL
#prepare index files for generic 13 | 14 (ion and metal) and for everything else (protein + water + ions)

# this selects ANION + METAL, stores as 26
printf "13 | 14 \n q \n" | gmx make_ndx -f md2.gro  -o index_hbond.ndx

#this selects PROTEIN + WATER + IONS, stores as 27
printf "1 | 25 \n q \n" | gmx make_ndx -f md2.gro  -n index_hbond.ndx -o index_hbond.ndx

#analyse hbonding between ANION + METAL (26) and PROTEIN + WATER + IONS (27)
printf "26 \n27 \nq \n" | gmx hbond  -n index_hbond.ndx -f md2_fitted.xtc -s md2.tpr -g hbond.log -dist -ang -hx -hbn -hbm -dan hbdanum.xvg -nhbdist -dt 20000
#rename 1st shell results
for f in hb*; do mv $f FIRST_$f; done

#hbond analysis from htf_md, SECOND SHELL
more FIRST_hbond.log | sed 1d | awk '{print $1; print $2; print $3}' | sed '/SOL/d' | sed '/338/d' | sed '/339/d' | cut -c-6  | cut -c 4- | uniq | awk '{print "r "$1" |"}'  | tr -d [A-Z] | sed ':a;N;$!ba;s/\n/ /g'|  rev | cut -c 2- | rev > maker
more FIRST_hbond.log | sed 1d | awk '{print $1; print $2; print $3}' | sed '/SOL/d' | sed '/338/d' | sed '/339/d' | cut -c-6  | cut -c 4- | uniq |  awk '{print " &! r "$1 }' | tr -d [A-Z] | sed ':a;N;$!ba;s/\n/ /g' | sed  's/^/1/' >> maker
echo "name 28 group28" >> maker
echo "name 29 group29" >> maker
echo "q" >> maker
#groups 29 and 30
gmx make_ndx -f md2.gro -n index_hbond.ndx -o index_hbond.ndx < maker
printf "28 \n 29 \n q \n" | gmx hbond  -n index_hbond.ndx -f md2_fitted.xtc -s md2.tpr -g hbond.log -dist -ang -hx -hbn -hbm -dan hbdanum.xvg -nhbdist -dt 20000
#rename 2nd shell results
# for f in hb*1.*; do rename 's/\./\_2./' $f; done
for f in hb*; do mv $f SECOND_$f; done

for f in FIRST*; do rename 's/\./\_1./' $f; done
for f in FIRST*; do rename 's/FIRST_//' $f; done
for f in SECOND*; do rename 's/\./\_2./' $f; done
for f in SECOND*; do rename 's/SECOND_//' $f; done


source /workopt/gromacs/bin/GMXRC.bash

