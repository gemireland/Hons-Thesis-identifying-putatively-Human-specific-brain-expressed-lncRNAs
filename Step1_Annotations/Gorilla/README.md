Python3 scripts for processing Gorilla's reference annotation GTF file (Ensembl release 116)

This folder has the Python3 scripts used to process Gorilla's reference genome annotation, to identify protein-coding genes (PCGs), long non-coding RNAs (lncRNAs) and PCG flanking anchors for each human lncRNA

The scripts:

1. gorillastep1.py: Extracts the annotated Gorilla lncRNAs and PCGs from Gorilla's reference GTF file and groups transcript isoforms into gene-level genomic blocks. Identifies neighbouring upstream and downstream PCGs.

2. Gorilla_PCG_anchors.py: Extracts PCG anchor information from the output of gorillastep1.py and generates Gorilla PCG genomic blocks in BED format

3. make_fullgor_lncRNA_anchors.py: Combines the Gorilla lncRNA genomic coordinates with the upstream and downstream PCG anchor information to generate a full lncRNA anchor table

4. AllGorillaPCGs.py: Refined script encompassing scripts 1, 2 and 3. Extracts all annotated Gorilla PCGs from the reference Ensembl GTF file and generates a BED file containing their genomic coordinates for subsequent BLAST analysis

alllncRNAsGorilla.py: extracts all annotated Gorilla lncRNAs from reference Ensembl GTF file and groups their transcript isoforms into gene-level blocks to create tab-delimited catalogues of lncRNAs without PCG anchor information. Output is used in synteny-based conservation

Input files (reference annotation): Gorilla_gorilla.gorGor4.116.gtf

Important to note: These scripts are the main processors of the Gorilla annotation processing stage of the pipeline, some intermediate scripts and details will need to be confirmed for full reproducibility
