Python3 scripts for human-brain expression analysis and integration with synteny-based conservation results

This folder contains the Python3 scripts used to identify human lncRNAs expressed in the brain using BrainSpan, developmental Genotype-Tissue Expression (dGTEx) and Genotype-Tissue Expression (GTEx) datasets. The resulting non-redundant catalogue was combined with the synteny-based conservation results to identify putatively human-specific lncRNAs.

The Python3 scripts:

1. humlnc_brainspan_expression.py: Identifies annotated human lncRNAs expressed in BrainSpan samples using an expression threshold of RPKM > 1

2. extractingbraindGTEx.py: used to extract gene-level data from brain tissue samples in the dGTEx dataset

3. dGTEx_lncRNA_geneexpression.py: used to identify annotated human lncRNAs expressed in dGTEx brain samples using an expression threshold of TPM>1

4. GTEx_lncRNA_expression.py: used to identify  annotated human lncRNAs expressed in adult GTEx brain samples using gene-level expression data

5. combined_humbrain_expressed.py: combines BrainSpan, dGTEx and GTEx expression results to generate a non-redundant catalogue of putatively human-specific lncRNAs with evidence of expression in the brain

6. classify_human_brain_lncRNAs.py: integrates the putatively human-specific lncRNAs with evidence of expression in the brain catalogue with the expanded Chimpanzee and Gorilla synteny-based conservation results and the original Mouse conservation results to classify human lncRNAs as mammalian-conserved, primate-conserved or putatively human-specific.

Input files: Human reference genome (GENCODE v49) lncRNA genomic blocks, BrainSpan expression data, dGTEx and GTEx gene-level expression matrices and sample metadata, and synteny-based conservation results for Chimpanzee, Gorilla and Mouse

Important to note: Human lncRNAs classified as putatively human-specific lacked detectable syntenic counterparts in the three comparison species using our developed computational pipeline. This classification does not definitively establish their absence for all other species.
