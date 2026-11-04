#!/bin/bash
#SBATCH --partition=ada
#SBATCH --job-name="Human_vs_Mouse_BLAST"
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=human_vs_mouse_123456.out
#SBATCH --error=human_vs_mouse_123456.err

#only if blastn is installed in blast_env
module load python/miniconda3-py3.12
source /opt/exp_soft/miniconda3-py3.12/etc/profile.d/conda.sh
conda activate bedtools_env    #this is the environmnet in which I installed BLASTN

cd /scratch/irlgem001/Step2-blast

blastn \
-query Human_PCG_Blocks.fa \
-db Mouse_db \
-out Human_PCG_vs_Mouse.tsv \
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
-max_target_seqs 5 \
-num_threads 8
