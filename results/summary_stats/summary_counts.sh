#!/bin/bash

#SBATCH --time=0-00:59
#SBATCH --mem=4000
#SBATCH --account=def-breton23
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1


cd /home/renaut/scratch/ATRAPP_Champlain_2016_Metat/

touch temp_counts

for file in sequences/hiseq_raw/*RNA?_R2.fastq.gz
	do
		zcat $file | wc -l >>temp_counts
	done

ls -1 sequences/hiseq_raw/*RNA?_R2.fastq.gz >temp_seq

paste temp_counts temp_seq >summary_stats/summary_counts

rm temp_counts temp_seq
