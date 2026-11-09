Python3 scripts for using BLASTN (v2.17.0) analysis to identify conserved PCGs 

This folder has the Python3 and Shell scripts used to identify conserved protein-coding genes (PCGs) anchors between Human and Mouse, Chimpanzee and Gorilla using BLASTN (v2.17.0). These results were incorporated into the synteny-based conservation analysis.

The scripts:

1. blaststep2MvH.sh: executes BLASTN analysis comparison (on the UCT HPC) between Human and Mouse PCG sequences using a species-specific PCG BLAST database

2. blaststep2CvH.sh: executes BLASTN analysis comparison (on the UCT HPC) between Human and Chimpanzee PCG sequences using a species-specific PCG BLAST database

3. blaststep2GvH.sh: executes BLASTN analysis comparison (on the UCT HPC) between Human and Gorilla PCG sequences using a species-specific PCG BLAST database

4. BLASTTop_hits.py: Python3 script used to process BLASTN alignment results for each species, combine their alignments and calculate statistics for the combined alignments for each Human-Comparison species PCG pair and ranks the matches (as stipulated in the ranking section of the script). This script identifies the highest-ranking PCG match (Top1) for each human PCG. 
   
Input files: PCG nucleotide sequences extracted from the corresponding reference genome assemblies using BEDTools getfasta

Important to note: Species-specific nucleotide BLAST databases were made from extracted PCG sequences. The highest-ranking PCG matches were further integrated with whole genome alignment data (WGA) and then integrated into the synteny-based conservation analysis 
