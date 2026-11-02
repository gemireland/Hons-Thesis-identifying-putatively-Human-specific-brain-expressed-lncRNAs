#!/bin/bash
#SBATCH --partition=ada
#SBATCH --job-name="Human_vs_Chimp_BLAST"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH --output=human_vs_chimp_123456.out
#SBATCH --error=human_vs_chimp_123456.err

module load python/miniconda3-py3.12

#only if blastn is installed in blast_env
source /opt/exp_soft/miniconda3-py3.12/etc/profile.d/conda.sh
conda activate bedtools_env #ths is the environment in which I installed BLAST

cd /scratch/irlgem001/Step2-blast

blastn \
-query Human_PCG_Blocks.fa \
-db Chimp_db \
-out Human_PCG_vs_Chimp.tsv \
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
-max_target_seqs 5 \
-num_threads 8
