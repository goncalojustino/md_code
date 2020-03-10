for f in protein_md*_[A-D]; do
 cp npt_template.mdp $f/npt.mdp
 done
 
for f in protein_md259* protein_md26[0-4]*; do 
 sed "s/XGRX/Protein_OXL_V3/" $f/npt.mdp -i
 done
for f in protein_md25[3-8]*; do 
 sed "s/XGRX/Protein_OXL_VO/" $f/npt.mdp -i
 done
for f in protein_md24[7-9]* protein_md25[0-2]*; do 
 sed "s/XGRX/Protein_OXL_FE3/" $f/npt.mdp -i
 done
for f in protein_md23[5-9]* protein_md24[0-6]*; do 
 sed "s/XGRX/Protein_ODA_VO/" $f/npt.mdp -i
 done
for f in protein_md22[3-9]* protein_md23[0-4]*; do 
 sed "s/XGRX/Protein_ODA_V3/" $f/npt.mdp -i
 done
for f in protein_md21[1-9]* protein_md22[0-2]*; do 
 sed "s/XGRX/Protein_ODA_FE3/" $f/npt.mdp -i
 done           
for f in protein_md20[5-9]* protein_md210*; do 
 sed "s/XGRX/Protein_LAC_VO/" $f/npt.mdp -i
 done   
for f in protein_md199* protein_md20[0-4]*; do 
 sed "s/XGRX/Protein_LAC_V3/" $f/npt.mdp -i
 done   
for f in protein_md19[3-8]*; do 
 sed "s/XGRX/Protein_LAC_FE3/" $f/npt.mdp -i
 done  
for f in protein_md18[7-9]* protein_md19[0-2]*; do 
 sed "s/XGRX/Protein_H1O_VO/" $f/npt.mdp -i
 done   
for f in protein_md17[8-9]* protein_md18[0-3]*; do 
 sed "s/XGRX/Protein_H1O_V3/" $f/npt.mdp -i
 done  
for f in protein_md17[2-7]*; do 
 sed "s/XGRX/Protein_H1O_FE3/" $f/npt.mdp -i
 done   
for f in protein_md16[6-9]* protein_md17[0-1]* protein_md15[4-6]*; do 
 sed "s/XGRX/Protein_AC2_VO/" $f/npt.mdp -i
 done   
for f in protein_md16[3-5]* protein_md14[8-9]* protein_md15[0-3]*; do 
 sed "s/XGRX/Protein_AC2_V3/" $f/npt.mdp -i
 done     
for f in protein_md15[7-9]* protein_md16[0-2]* protein_md14[5-7]*; do 
 sed "s/XGRX/Protein_AC2_FE3/" $f/npt.mdp -i
 done      
for f in protein_md139* protein_md14[0-4]*; do 
 sed "s/XGRX/Protein_HHO_VO/" $f/npt.mdp -i
 done   
for f in protein_md13[3-8]*; do 
 sed "s/XGRX/Protein_HHO_V3/" $f/npt.mdp -i
 done     
for f in protein_md12[7-9]* protein_md13[0-2]*; do 
 sed "s/XGRX/Protein_HHO_FE3/" $f/npt.mdp -i
 done       
for f in protein_md12[1-6]*; do 
 sed "s/XGRX/Protein_CTT_VO/" $f/npt.mdp -i
 done   
for f in protein_md11[5-9]* protein_md120*; do 
 sed "s/XGRX/Protein_CTT_V3/" $f/npt.mdp -i
 done     
for f in protein_md109* protein_md11[0-4]*; do 
 sed "s/XGRX/Protein_CTT_FE3/" $f/npt.mdp -i
 done   
for f in protein_md10[3-8]*; do 
 sed "s/XGRX/Protein_CT2_VO/" $f/npt.mdp -i
 done   
for f in protein_md09[7-9]* protein_md10[0-2]*; do 
 sed "s/XGRX/Protein_CT2_V3/" $f/npt.mdp -i
 done     
for f in protein_md09[1-6]*; do 
 sed "s/XGRX/Protein_CT2_FE3/" $f/npt.mdp -i
 done 
for f in protein_md08[2-9]* protein_md090*; do 
 sed "s/XGRX/Protein_AC1_VO/" $f/npt.mdp -i
 done   
for f in protein_md07[3-9]* protein_md08[0-1]*; do 
 sed "s/XGRX/Protein_AC1_V3/" $f/npt.mdp -i
 done     
for f in protein_md06[4-9]* protein_md07[0-2]*; do 
 sed "s/XGRX/Protein_AC1_FE3/" $f/npt.mdp -i
 done  
 
for f in protein_md001* protein_md007* protein_md029*; do
 sed "s/XGRX/Protein_COA_FE3/" $f/npt.mdp -i
 done  
for f in protein_md004* protein_md010* protein_md032*; do
 sed "s/XGRX/Protein_COA_V3/" $f/npt.mdp -i
 done 
for f in protein_md043* protein_md046* protein_md057*; do
 sed "s/XGRX/Protein_COA_VO/" $f/npt.mdp -i
 done  
for f in protein_md002* protein_md008* protein_md030*; do
 sed "s/XGRX/Protein_COA_FE3/" $f/npt.mdp -i
 done  
for f in protein_md005* protein_md011* protein_md033*; do
 sed "s/XGRX/Protein_COB_V3/" $f/npt.mdp -i
 done 
for f in protein_md044* protein_md047* protein_md058*; do
 sed "s/XGRX/Protein_COB_VO/" $f/npt.mdp -i
 done   
for f in protein_md003* protein_md009* protein_md031*; do
 sed "s/XGRX/Protein_COC_FE3/" $f/npt.mdp -i
 done  
for f in protein_md006* protein_md012* protein_md034*; do
 sed "s/XGRX/Protein_COC_V3/" $f/npt.mdp -i
 done 
for f in protein_md045* protein_md048* protein_md059*; do
 sed "s/XGRX/Protein_COC_VO/" $f/npt.mdp -i
 done    
for f in protein_md018* protein_md022* protein_md038*; do
 sed "s/XGRX/Protein_H2P_FE3/" $f/npt.mdp -i
 done  
for f in protein_md020* protein_md024* protein_md040*; do
 sed "s/XGRX/Protein_H2P_V3/" $f/npt.mdp -i
 done 
for f in protein_md052* protein_md054* protein_md062*; do
 sed "s/XGRX/Protein_H2P_VO/" $f/npt.mdp -i
 done    
for f in protein_md017* protein_md021* protein_md037*; do
 sed "s/XGRX/Protein_HPO_FE3/" $f/npt.mdp -i
 done  
for f in protein_md019* protein_md023* protein_md039*; do
 sed "s/XGRX/Protein_HPO_V3/" $f/npt.mdp -i
 done 
for f in protein_md051* protein_md053* protein_md061*; do
 sed "s/XGRX/Protein_HPO_VO/" $f/npt.mdp -i
 done  
for f in protein_md025* protein_md027* protein_md041*; do
 sed "s/XGRX/Protein_NO3_FE3/" $f/npt.mdp -i
 done  
for f in protein_md026* protein_md028* protein_md042*; do
 sed "s/XGRX/Protein_NO3_V3/" $f/npt.mdp -i
 done 
for f in protein_md055* protein_md056* protein_md063*; do
 sed "s/XGRX/Protein_NO3_VO/" $f/npt.mdp -i
 done  
for f in protein_md013* protein_md015* protein_md035*; do
 sed "s/XGRX/Protein_SO4_FE3/" $f/npt.mdp -i
 done  
for f in protein_md014* protein_md016* protein_md036*; do
 sed "s/XGRX/Protein_SO4_V3/" $f/npt.mdp -i
 done 
for f in protein_md049* protein_md050* protein_md060*; do
 sed "s/XGRX/Protein_SO4_VO/" $f/npt.mdp -i
 done    
 
