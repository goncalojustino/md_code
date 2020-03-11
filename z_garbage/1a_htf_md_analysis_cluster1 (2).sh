for f in protein_md{001..020}_A
do

cd /work6/$f
#house-keeping
rm \#* cluster*

#first things first - re-order water molecules by distance to metal+ion
echo "13 | 14" > maker #this is metal+ion
echo "q" >> maker
gmx  make_ndx -f md.gro -o order.ndx < maker
echo "26" > maker
echo "17" >> maker
gmx  trjorder -f md.xtc -s md.gro -n order.ndx -nshell nshell.xvg -o md_ordered.xtc -xvg xmgrace -r 0.5 -dt 1 < maker

#group together Protein + Metal + Ion
echo "1 | 13 | 14 " > maker
echo "q" >> maker
gmx make_ndx -f md.gro -o index_cluster.ndx < maker

#remove roto+trans; this selects protein for least squares fit and system for output
echo "1" > maker
echo "0" >> maker
gmx trjconv -f md_ordered.xtc -s md.tpr -o  md_fitted.xtc -fit rot+trans -n index_cluster.ndx -dt 1 < maker

#cluster by Protein+metal+ion and output full system clusters
echo "26" > maker
echo "0" >> maker
gmx cluster -f md_fitted.xtc -s md.tpr -method gromos -g cluster_out.log -dt 1 -cl clusters.pdb -n index_cluster.ndx -wcl 5 < maker

#split multimodel pdb file into various single model pdb files
#csplit -k -s -n 3 -b '%02d.pdb' -f cluster_ clusters.pdb '/^ENDMDL/+1' {*}

#correct an atom nomenclature issue from PDB format
for g in cluster*pdb
 do
 gmx editconf -f $g -o temp.pdb
 sed "s/0.00           F/0.00          FE/" temp.pdb -i
 mv temp.pdb $g
 done

done
