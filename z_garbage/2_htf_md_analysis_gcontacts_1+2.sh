#!\bin\bash
#on md_ordered.xtc in all folders

#for f in protein_md{001..063}_A                    # CORRECT FOLDERS
for f in protein_md036_A protein_md038_A protein_md{040..042}_A protein_md063_A
 do
  cd /work6/$f       # CORRECT LOCATION
  rm \#*             # house cleaning never hurt a soul
 
 source /workopt/gromacs46/bin/GMXRC.bash

# echo "17|18" > maker
# echo "q" >> maker
# make_ndx_46 -f md.gro -o temp.ndx < maker

 #this compares metal+ion vs. all system (water!)
 #and does not slow down the child
# echo "34" > maker # this is metal + ion
# echo "0" >> maker

# g_contacts -f md_ordered.xtc -n temp.ndx -dt 1 -resndx -s md.gro -d 0.3 -on g_cont_res.ndx -o g_contacts.dat < maker
 
 #get contacts to second coordination 
# more g_contacts.dat | sed 1d | sed '/SOL/d' | awk '{print $7}' | sort | uniq  | awk '{print "r "$1" |"}' | sed ':a;N;$!ba;s/\n/ /g'  |  rev | cut -c 2- | rev > maker
# more g_contacts.dat | sed 1d | sed '/SOL/d' | awk '{print $7}' | sort | uniq  | awk '{print " &! r "$1 }' | sed ':a;N;$!ba;s/\n/ /g' | sed  's/^/1/' >> maker
# echo q >> maker
 make_ndx_46 -f md.gro -o temp.ndx < maker

 echo "34" > maker 
 echo "35" >> maker
 g_contacts -f md_ordered.xtc -n temp.ndx -dt 1 -resndx -s md.gro -d 0.35 -on g_cont_res2.ndx -o g_contacts2.dat < maker


 source /workopt/gromacs/bin/GMXRC.bash

 done

