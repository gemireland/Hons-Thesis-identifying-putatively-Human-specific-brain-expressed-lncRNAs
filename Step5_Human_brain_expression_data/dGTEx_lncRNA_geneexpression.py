#!/usr/bin/env python3

#now we want to compare the gene_ids from our lncRNA blcoks to the gene_ids in teh dGTEx brain expression file we created from extracting teh sample IDs for the huge dGTEx expression file

#firstly we are using pandas again to help with ease of working in a table
import pandas as pd

#to read in the human lncRNA blocks file
human_lncRNA_blocks = pd.read_csv("Human_lncRNA_blocks.tsv", sep="\t")

#reading in our extracted brain expression file
dGTEx_brain_expression = pd.read_csv("dGTEx_onlybrain_genes_TPM.tsv", sep="\t")

#remember as it is the same as the BrainSpan expression file, we do not need the .2 suffix on the ids so remove those

human_lncRNA_blocks["gene_id"] = human_lncRNA_blocks["lncrna_gene_id"].str.split(".").str[0]

print(dGTEx_brain_expression.columns.tolist())
dGTEx_brain_expression["gene_id_refined"] = dGTEx_brain_expression["gene_id"].str.split(".").str[0]

#these are the names of the columns of our 4 brain samples from the expression file
brain_sample_columns = [ 
    "DGTEX001-2BZ-01-SM-OLAQB",
    "DGTEX002-2BZ-01-SM-OLMRP",
    "DGTEX031-2BZ-01-SM-OLASE",
    "DGTEX047-2BZ-01-SM-OLATP"
]

#find the maximum expression value for each gene in each of the samples
dGTEx_brain_expression["max_expression"] = dGTEx_brain_expression[brain_sample_columns].max(axis=1)

#retain only the transcripts with TPM > 1
expressed_transcripts = dGTEx_brain_expression[dGTEx_brain_expression["max_expression"]>1]

#get the refined gene_ids that match the expressed transcripts
expressed_genes = set(expressed_transcripts["gene_id_refined"])

#now we want to see which of our human lncRNA blocks are expressed in the dGTEx brain expression file
human_lncRNA_blocks["dGTEx_brain_expression"] = human_lncRNA_blocks["gene_id"].isin(expressed_genes)

#output of all lncs with expression status 
human_lncRNA_blocks.to_csv("Human_lncRNAs_dGTEx_brain_expression_status.tsv", sep="\t", index=False)

#now only saving the expressed ones
expressed_human_lncRNAs = human_lncRNA_blocks[human_lncRNA_blocks["dGTEx_brain_expression"]== True]

expressed_human_lncRNAs.to_csv("Human_lncRNAs_expressed_in_dGTEx_brain.tsv", sep="\t", index=False)

print("total human lncRNAs:", len(human_lncRNA_blocks))
print("Human lncRNAs expressed in dGTEx brain:", len(expressed_human_lncRNAs))





