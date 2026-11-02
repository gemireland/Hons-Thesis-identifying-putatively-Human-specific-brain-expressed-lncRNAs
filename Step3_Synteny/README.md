Python3 scripts for synteny-based lncRNA conservation analysis 

This folder has the Python3 scripts used to identify syntenically conserved lncRNAs between Human and Mouse, Chimpanzee and Gorilla by integrating conserved PCG matches and WGA regions in addition to previously identified lncRNA genomic coordinates

The Python3 scripts:

1. WGASynNet.py: Executes synteny-based conservation analysis on Human-Mouse comparison by integrating PCG matches, WGA regions and lncRNA genomic blocks to identify conserved lncRNA candidates 

2. WGASynNetChimp.py: Executes synteny-based conservation analysis on Human-Chimpanzee comparison by integrating PCG matches, WGA regions and lncRNA genomic blocks to identify conserved lncRNA candidates
   
3. WGASynNetGorilla.py: Executes synteny-based conservation analysis on Human-Gorilla comparison by integrating PCG matches, WGA regions and lncRNA genomic blocks to identify conserved lncRNA candidates

4. cleanupConservedlncs.py: applies an additional filtering step to the conserved lncRNA candidates (described in the scripts comments) to create a final set of syntenically conserved lncRNAs
   
Input files: Human and comparison species lncRNA genomic blocks, BLASTN Top1 PCG matches and UCSC pairwise Syntenic-Network MAF (WGA) files

Important to note: the conserved lncRNAs idnetified were candidates and later on in the study were investigated using brain-expression data to distinguish conserved from putatively human-specific lncRNAs
