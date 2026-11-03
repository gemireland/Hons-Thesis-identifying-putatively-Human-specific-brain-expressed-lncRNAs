Python3 and shell scripts for RNA-seq based brain transcriptome assembly and expansion of Chimpanzee and Gorilla lncRNA catalogues. 

This folder has the Python3 and shell scripts used to process the StringTie assembled Chimpanzee and Gorilla transcriptomes, identify previously unannotated lncRNAs and to expand the existing lncRNA catalogues for both Chimpanzee and Gorilla.

The scripts:

1. Chimp_StringTiegff_tobed.py: used to convert the newly assembled StringTie transcript annotations for Chimpanzee into BED format from GTF format

2. Gorilla_StringTiegff_tobed.py: used to convert the newly assembled StringTie transcript annotations for Gorilla into BED format from GTF format

3. createlncblocksfromStringTie.py: used to group all transcript isoform coordinates for each gene locus into genomic blocks

4. filterfeelncandcpc2: takes the results from coding potential analysis (CPC2 and FEELnc) and analyses them to identify non-coding transcripts

5. makenovelintergenicstogtf.py: used to convert the identified intergenic transcripts (from gffcompare) into GTF format

6. finalannotatedlncblocksChimp.py: used to process the expanded Chimpanzee lncRNA annotations and create gene-level lncRNA blocks (reference annotated and StringTie identified) for synteny-based reassessment

7. finalannotatedlncblocksGorilla.py: used to process the expanded Gorilla lncRNA annotations and create gene-level lncRNA blocks (reference annotated and StringTie identified) for synteny-based reassessment

Input files: Chimpanzee and Gorilla reference genome annotations (Ensembl release 116), newly assembled StringTie transcriptomes for Chimpanzee and Gorilla, results from CPC2 and FEELnc coding-potential analysis.

Important to note: These scripts represent only a fraction of the processing of the newly assembled transcriptomes. The initial RNA-seq alignment and transcriptome assembly were performed using HISAT2 and StringTie, while gffcomapre was used to compare reconstructed transcripts against the corresponding reference annotations. Full methodology can be followed in Section 2.5 of the Materials and Methodology section.  
