Python3 scripts for processing Chimpanzee's reference annotation GTF file (Ensembl release 116)

This folder has the Python3 scripts used to process Chimpanzee's reference genome annotation, to identify protein-coding genes (PCGs), long non-coding RNAs (lncRNAs) and PCG flanking anchors for each human lncRNA

The scripts:

1. finalchimpensstep1.py: Extracts the annotated human lncRNAs and PCGs from Chimpanzee's reference GTF file and groups transcript isoforms into gene-level genomic blocks. Identifies neighbouring upstream and downstream PCGs.

2. finalchimp_PCG_anchors.py: Extracts PCG anchor information from the output of finalchimpensstep1.py and generates Chimpanzee PCG genomic blocks in BED format

3. make_fullchimp_lncRNA_anchors.py: Combines the Chimpanzee lncRNA genomic coordinates with the upstream and downstream PCG anchor information to generate a full lncRNA anchor table
   
4. allChimpPCGs.py: Refined script encompassing scripts 1, 2 and 3. Extracts all annotated Chimpanzee PCGs from the reference Ensembl GTF file and generates a BED file containing their genomic coordinates for subsequent BLAST analysis
   
5. alllncRNAsChimp.py: extracts all annotated Chimpanzee lncRNAs from reference Ensembl GTF file and groups their transcript isoforms into gene-level blocks to create tab-delimited catalogues of lncRNAs without PCG anchor information. Output is used in synteny-based conservation

Input files (reference annotation): Pan_troglodytes.Pan_tro_3.0.116.gtf

Important to note: These scripts are the main processors of the Chimpanzee annotation processing stage of the pipeline, some intermediate scripts and details will need to be confirmed for full reproducibility
