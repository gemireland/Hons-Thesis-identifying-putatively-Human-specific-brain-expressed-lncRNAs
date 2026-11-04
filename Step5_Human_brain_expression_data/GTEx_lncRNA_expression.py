#!/usr/bin/env python3

#using pandas to easily work with table data (.tsv file)
import pandas as pd

#read human lncRNA blocks
human_lncRNA_blocks = pd.read_csv("Human_lncRNA_blocks.tsv", sep="\t")

#read Adult GTEx brain expression data
GTEx_brain_expression = pd.read_csv("GTEx_onlybrain_genes_TPM.tsv", sep="\t")

#remove Ensembl version suffix from IDs so that we can compare across the datasets 
human_lncRNA_blocks["gene_id"] = (human_lncRNA_blocks["lncrna_gene_id"].str.split(".").str[0])

GTEx_brain_expression["gene_id_refined"] = (GTEx_brain_expression["Name"].str.split(".").str[0])

#identify all brain sample columns to remove the non-brain sample tissues
brain_sample_columns = [
    column for column in GTEx_brain_expression.columns
    if column not in [
        "Name",
        "Description",
        "gene_id_refined"
    ]
]


#find maximum expression across all brain samples
GTEx_brain_expression["max_expression"] = (
    GTEx_brain_expression[brain_sample_columns]
    .max(axis=1)
)

#retain only genes with TPM > 0.1 (did this the first time but ended up with too many which didn't make biological sense)
#trying again with a TPM > 1 as 0.1 resulted in over 33 000 human lncRNAs being expressed in the brain (so refining our threshold)
#we ended up going with threshold of TPM > 1
expressed_genes_GTEx = GTEx_brain_expression[
    GTEx_brain_expression["max_expression"] > 1
]

#get gene IDs of expressed genes (TMP>1)
expressed_genes = set(expressed_genes_GTEx["gene_id_refined"])

#see which expressed genes are in our human lncRNA blocks
human_lncRNA_blocks["GTEx_brain_expression"] = (
    human_lncRNA_blocks["gene_id"]
    .isin(expressed_genes)
)

#save all lncRNAs with evidence of expression
human_lncRNA_blocks.to_csv(
    "Human_lncRNAs_GTEx_brain_expression_status.tsv",
    sep="\t",
    index=False
)

#save only lncRNAs that are expressed within the brain
expressed_human_lncRNAs = human_lncRNA_blocks[
    human_lncRNA_blocks["GTEx_brain_expression"] == True]

expressed_human_lncRNAs.to_csv(
    "Human_lncRNAs_expressed_in_GTEx_brain.tsv",
    sep="\t",
    index=False)


print("total no. human lncRNAs",len(human_lncRNA_blocks))

print("no. human lncRNAs expressed in Adult GTEx brain",len(expressed_human_lncRNAs))
