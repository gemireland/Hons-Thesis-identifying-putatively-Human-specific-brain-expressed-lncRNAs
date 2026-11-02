Python3 scripts used to process the reference genome annotation for Mouse (GENCODE vM38) 

These scripts in this folder are used to process the reference annotation (GENCODE vM38) to idnetify annotated lncRNAs, PCGs and flanking/neighbouring PCG anchors 
for use in synteny-based conservation analysis 

The Python3 scripts:

1. mousegencodestep1.py: Initial script used to extract annotated mouse lncRNAs and PCGs from GENCODE vM38 reference annotation, groups all transcript isoforms into gene-level genomic blocks and identifies upstream and downstream PCG blocks

2. AllMousePCGs.py: Used to extract annotated mouse PCGs from the reference GTF file for BLASTN analysis 

3. alllncRNAsM.py: Refined script encompassing scripts 1 and 2 to extract annotated mouse lncRNAs and group the individual transcript isoforms into gene-level genomic blocks 

4. Make_Mouse_lncRNA_Full_Anchor_Table.py: Combines mouse lncRNA genomic coordinates with their upstream and downstream PCG anchor information to create a lncRNA and PCG anchor table

Input file: 
- GENCODE vM38: gencode.vM38.annotation.gtf

Note: These scripts form the bulk of the Mouse reference genome processing step, intermediary files and scripts may need to be confirmed for exact reproducibility of the pipeline
