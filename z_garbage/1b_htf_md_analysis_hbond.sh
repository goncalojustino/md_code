for f in protein_md{001..063}_A # CORRECT FOLDERS
 do
  cd /work6/$f    # CORRECT LOCATION
  rm \#*          # house cleaning never hurt a soul

#hbond analysis from htf_md
#prepare index files for generic 13 | 14 (ion and metal) and for everything else (protein + water + ions)
echo "13 | 14 " > maker
echo "1 | 25 " >> maker
echo "q" >> maker
gmx make_ndx -f md.gro -o index_hbond.ndx < maker
#run bond analysis
echo " 26 " > maker
echo " 27 " >> maker
gmx hbond  -n index_hbond.ndx -f md_fitted.xtc -s md.tpr -g hbond.log -dist -ang -hx -hbn -hbm -dan -nhbdist < maker

#2nd shell
mv hbond.ndx hbond1.ndx
mv hbond.log hbond1.log
mv hbnum.xvg hbnum1.xvg
mv hbdist.xvg hbdist1.xvg
mv hbang.xvg hbang1.xvg
mv hbhelix.xvg hbhelix1.xvg
mv hbmap.xpm hbmap1.xpm
mv danum.xvg danum1.xvg

more hbond1.log | sed 1d | awk '{print $1; print $2; print $3}' | sed '/SOL/d' | sed '/338/d' | sed '/339/d' | cut -c-6  | cut -c 4- | uniq | awk '{print "r "$1" |"}'  | tr -d [A-Z] | sed ':a;N;$!ba;s/\n/ /g'|  rev | cut -c 2- | rev > maker
more hbond1.log | sed 1d | awk '{print $1; print $2; print $3}' | sed '/SOL/d' | sed '/338/d' | sed '/339/d' | cut -c-6  | cut -c 4- | uniq |  awk '{print " &! r "$1 }' | tr -d [A-Z] | sed ':a;N;$!ba;s/\n/ /g' | sed  's/^/1/' >> maker
echo "q" >> maker
gmx make_ndx -f md.gro -o index_hbond2.ndx < maker
echo "26 " > maker
echo "27 " >> maker
gmx hbond  -n index_hbond2.ndx -f md_fitted.xtc -s md.tpr -g hbond.log -dist -ang -hx -hbn -hbm -dan -nhbdist -dt 1 < maker

mv hbond.ndx hbond2.ndx
mv hbond.log hbond2.log
mv hbnum.xvg hbnum2.xvg
mv hbdist.xvg hbdist2.xvg
mv hbang.xvg hbang2.xvg
mv hbhelix.xvg hbhelix2.xvg
mv hbmap.xpm hbmap2.xpm
mv danum.xvg danum2.xvg

done
