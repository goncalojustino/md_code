for f in protein_md259* protein_md26[0-4]*; do 
 echo " \"Protein\" | \"OXL\" | \"V3\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done
for f in protein_md25[3-8]*; do 
 echo " \"Protein\" | \"OXL\" | \"VO\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done
for f in protein_md24[7-9]* protein_md25[0-2]*; do 
 echo " \"Protein\" | \"OXL\" | \"FE3\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done
for f in protein_md23[5-9]* protein_md24[0-6]*; do 
 echo " \"Protein\" | \"ODA\" | \"VO\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done
for f in protein_md22[3-9]* protein_md23[0-4]*; do 
 echo " \"Protein\" | \"ODA\" | \"V3\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done
 for f in protein_md21[1-9]* protein_md22[0-2]*; do 
 echo " \"Protein\" | \"ODA\" | \"FE3\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done           
for f in protein_md20[5-9]* protein_md210*; do 
 echo " \"Protein\" | \"LAC\" | \"VO\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done   
for f in protein_md199* protein_md20[0-4]*; do 
 echo " \"Protein\" | \"LAC\" | \"V3\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done   
for f in protein_md19[3-8]*; do 
 echo " \"Protein\" | \"LAC\" | \"FE3\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done  
for f in protein_md18[7-9]* protein_md19[0-2]*; do 
 echo " \"Protein\" | \"H1O\" | \"VO\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done   
for f in protein_md17[8-9]* protein_md18[0-3]*; do 
 echo " \"Protein\" | \"H1O\" | \"V3\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done  
for f in protein_md17[2-7]*; do 
 echo " \"Protein\" | \"H1O\" | \"FE3\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done   
for f in protein_md16[6-9]* protein_md17[0-1]* protein_md15[4-6]*; do 
 echo " \"Protein\" | \"AC2\" | \"VO\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done   
for f in protein_md16[3-5]* protein_md14[8-9]* protein_md15[0-3]*; do 
 echo " \"Protein\" | \"AC2\" | \"V3\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done     
for f in protein_md15[7-9]* protein_md16[0-2]* protein_md14[5-7]*; do 
 echo " \"Protein\" | \"AC2\" | \"FE3\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done      
for f in protein_md139* protein_md14[0-4]*; do 
 echo " \"Protein\" | \"HHO\" | \"VO\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done   
for f in protein_md13[3-8]*; do 
 echo " \"Protein\" | \"HHO\" | \"V3\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done     
for f in protein_md12[7-9]* protein_md13[0-2]*; do 
 echo " \"Protein\" | \"HHO\" | \"FE3\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done       
for f in protein_md12[1-6]*; do 
 echo " \"Protein\" | \"CTT\" | \"VO\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done   
for f in protein_md11[5-9]* protein_md120*; do 
 echo " \"Protein\" | \"CTT\" | \"V3\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done     
for f in protein_md109* protein_md11[0-4]*; do 
 echo " \"Protein\" | \"CTT\" | \"FE3\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done   
for f in protein_md10[3-8]*; do 
 echo " \"Protein\" | \"CT2\" | \"VO\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done   
for f in protein_md09[7-9]* protein_md10[0-2]*; do 
 echo " \"Protein\" | \"CT2\" | \"V3\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done     
for f in protein_md09[1-6]*; do 
 echo " \"Protein\" | \"CT2\" | \"FE3\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done 
for f in protein_md08[2-9]* protein_md090*; do 
 echo " \"Protein\" | \"AC1\" | \"VO\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done   
for f in protein_md07[3-9]* protein_md08[0-1]*; do 
 echo " \"Protein\" | \"AC1\" | \"V3\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done     
for f in protein_md06[4-9]* protein_md07[0-2]*; do 
 echo " \"Protein\" | \"AC1\" | \"FE3\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done  
 
for f in protein_md001* protein_md007* protein_md029*; do
 echo " \"Protein\" | \"COA\" | \"FE3\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done  
for f in protein_md004* protein_md010* protein_md032*; do
 echo " \"Protein\" | \"COA\" | \"V3\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done 
for f in protein_md043* protein_md046* protein_md057*; do
 echo " \"Protein\" | \"COA\" | \"VO\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done  
for f in protein_md002* protein_md008* protein_md030*; do
 echo " \"Protein\" | \"COB\" | \"FE3\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done  
for f in protein_md005* protein_md011* protein_md033*; do
 echo " \"Protein\" | \"COB\" | \"V3\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done 
for f in protein_md044* protein_md047* protein_md058*; do
 echo " \"Protein\" | \"COB\" | \"VO\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done   
for f in protein_md003* protein_md009* protein_md031*; do
 echo " \"Protein\" | \"COC\" | \"FE3\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done  
for f in protein_md006* protein_md012* protein_md034*; do
 echo " \"Protein\" | \"COC\" | \"V3\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done 
for f in protein_md045* protein_md048* protein_md059*; do
 echo " \"Protein\" | \"COC\" | \"VO\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done    
for f in protein_md018* protein_md022* protein_md038*; do
 echo " \"Protein\" | \"H2P\" | \"FE3\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done  
for f in protein_md020* protein_md024* protein_md040*; do
 echo " \"Protein\" | \"H2P\" | \"V3\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done 
for f in protein_md052* protein_md054* protein_md062*; do
 echo " \"Protein\" | \"H2P\" | \"VO\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done    
for f in protein_md017* protein_md021* protein_md037*; do
 echo " \"Protein\" | \"HPO\" | \"FE3\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done  
for f in protein_md019* protein_md023* protein_md039*; do
 echo " \"Protein\" | \"HPO\" | \"V3\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done 
for f in protein_md051* protein_md053* protein_md061*; do
 echo " \"Protein\" | \"HPO\" | \"VO\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done  
for f in protein_md025* protein_md027* protein_md041*; do
 echo " \"Protein\" | \"NO3\" | \"FE3\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done  
for f in protein_md026* protein_md028* protein_md042*; do
 echo " \"Protein\" | \"NO3\" | \"V3\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done 
for f in protein_md055* protein_md056* protein_md063*; do
 echo " \"Protein\" | \"NO3\" | \"VO\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done  
for f in protein_md013* protein_md015* protein_md035*; do
 echo " \"Protein\" | \"SO4\" | \"FE3\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done  
for f in protein_md014* protein_md016* protein_md036*; do
 echo " \"Protein\" | \"SO4\" | \"V3\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done 
for f in protein_md049* protein_md050* protein_md060*; do
 echo " \"Protein\" | \"SO4\" | \"VO\" " > $f/index_creator_em
 echo "q" >> $f/index_creator_em
 done    
 
