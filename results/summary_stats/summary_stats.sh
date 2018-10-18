#how many PE reads were assembled
grep -F 'Assembled reads .' /home/renaut/scratch/ATRAPP_Champlain_2016_Metat/samsa2_run_St1/step_1_merging/pear_log >temp1
grep -F 'Assembled reads file.' /home/renaut/scratch/ATRAPP_Champlain_2016_Metat/samsa2_run_St1/step_1_merging/pear_log >temp1_n
grep -F 'Assembled reads .' /home/renaut/scratch/ATRAPP_Champlain_2016_Metat/samsa2_run_St2/step_1_merging/pear_log >temp2
grep -F 'Assembled reads file.' /home/renaut/scratch/ATRAPP_Champlain_2016_Metat/samsa2_run_St2/step_1_merging/pear_log >temp2_n
grep -F 'Assembled reads .' /home/renaut/scratch/ATRAPP_Champlain_2016_Metat/samsa2_run_PRM/step_1_merging/pear_log >temp3
grep -F 'Assembled reads file.' /home/renaut/scratch/ATRAPP_Champlain_2016_Metat/samsa2_run_PRM/step_1_merging/pear_log >temp3_n


paste temp1 temp1_n >temp1_nn
paste temp2 temp2_n >temp2_nn
paste temp3 temp3_n >temp3_nn

cat temp1_nn temp2_nn temp3_nn >pe_assembly

rm temp*

#how many merged reads were dropped
grep 'Dropped' ../samsa2_run_St1/step_2_trimming/trimmomatic_log >temp1
grep 'Dropped' ../samsa2_run_St2/step_2_trimming/trimmomatic_log >temp2
grep 'Dropped' ../samsa2_run_PRM/step_2_trimming/trimmomatic_log >temp3
grep 'phred33' ../samsa2_run_St1/step_2_trimming/trimmomatic_log >temp1_n 
grep 'phred33' ../samsa2_run_St2/step_2_trimming/trimmomatic_log >temp2_n
grep 'phred33' ../samsa2_run_PRM/step_2_trimming/trimmomatic_log >temp3_n

paste temp1 temp1_n >temp1_nn
paste temp2 temp2_n >temp2_nn
paste temp3 temp3_n >temp3_nn

cat temp1_nn temp2_nn temp3_nn >trimmed

rm temp*

#how many reads were IDed as ribosomes
grep 'passing' /home/renaut/scratch/ATRAPP_Champlain_2016_Metat/samsa2_run_St1/step_2_trimming/*ribosomes.log >ribosomes
grep 'passing' /home/renaut/scratch/ATRAPP_Champlain_2016_Metat/samsa2_run_St2/step_2_trimming/*ribosomes.log >>ribosomes
grep 'passing' /home/renaut/scratch/ATRAPP_Champlain_2016_Metat/samsa2_run_PRM/step_2_trimming/*ribosomes.log >>ribosomes

#how many reads were annotated Refseq
wc -l /home/renaut/scratch/ATRAPP_Champlain_2016_Metat/samsa2_run_St1/step_4_annotation/*RefSeq_annotated >refseq
wc -l /home/renaut/scratch/ATRAPP_Champlain_2016_Metat/samsa2_run_St2/step_4_annotation/*RefSeq_annotated >>refseq
wc -l /home/renaut/scratch/ATRAPP_Champlain_2016_Metat/samsa2_run_PRM/step_4_annotation/*RefSeq_annotated >>refseq

#how many reads were annotated subsys
wc -l /home/renaut/scratch/ATRAPP_Champlain_2016_Metat/samsa2_run_St1/step_4_annotation/*subsys_annotated >subsys
wc -l /home/renaut/scratch/ATRAPP_Champlain_2016_Metat/samsa2_run_St2/step_4_annotation/*subsys_annotated >>subsys
wc -l /home/renaut/scratch/ATRAPP_Champlain_2016_Metat/samsa2_run_PRM/step_4_annotation/*subsys_annotated >>subsys


