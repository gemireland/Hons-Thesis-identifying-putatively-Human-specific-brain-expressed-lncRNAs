Python3 scripts for processing Human reference annotation GTF file (GENCODE v49)

This folder has the Python3 scripts used to process GENCODE v49, the reference genome annotation, to identify protein-coding genes (PCGs), long non-coding RNAs (lncRNAs) and PCG flanking anchors for each human lncRNA

The scripts:

1. humangencodestep1.py: Extracts the annotated human lncRNAs and PCGs from the GENCODE v49 GTF file and groups transcript isoforms into gene-level genomic blocks. Identifies neighbouring upstream and downstream PCGs.

2. Human_PCG_anchorsscript.py: Used to extract PCG blocks from the output of humangencodestep1.py and creates Human_PCG_Blocks.bed files used in BLASTN analysis

3. Merged_Human_lncRNA_andAnchors.py: Combines the human lncRNA genomic coordinates with PCGs upstream and downstream (its anchors) to generate Human_lncRNA_Full_Anchor_Table.tsv used in synteny-based conservation analysis.

4. alllncRNAsHum.py: Refined script encompassing scripts 1, 2 and 3. Extracts annoatted human lncRNA transcripts from GENCODE v49 and groups transcript isoforms into gene-level blocks to create tab-delimited catalogues of lncRNAs without PCG anchor information

Input files (reference annotation): gencode.v49.annotation.gtf

Output files: 
- Human_PCG_Blcoks.bed
- Human_lncRNA_Full_Anchor_Table.tsv
- Human_lncRNA_blocks.tsv

Important to note: These scripts are the main processors of the human annotation processing stage of the pipeline, some intermediate scripts and details will need to be confirmed for full reproducibility
