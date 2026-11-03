#!/usr/bin/env python3

# This script identifies human lncRNA genes expressed in brain samples
# from developmental GTEx (dGTEx) and Adult GTEx using gene-level TPM data.

import pandas as pd

#defining file paths 

#human GENCODE lncRNA gene blocks
lncrna_file = "Human_lncRNA_blocks.tsv"

#dGTEx

dgtex_expression_file = (
    "dGTEx_Analysis_2025-07-21_v1_RNASeQCv2.4.3_gene_tpm.gct"
)

dgtex_metadata_file = (
    "dGTEx_Analysis_2026-01-30_v1_Sample_Attributes_Dataset_DS_open.xlsx"
)

#adult GTEx

gtex_expression_file = (
    "GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_tpm_non_lcm.gct"
)

gtex_metadata_file = (
    "GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt"
)

#loading the human lncRNA catalogue

print("\nLoading human lncRNA catalogue...")

human_lncs = pd.read_csv(
    lncrna_file,
    sep="\t"
)

#remove Ensembl version suffix so they can be comparable
#example: ENSG00000123456.4 -> ENSG00000123456

human_lncs["gene_id"] = (
    human_lncs["lncrna_gene_id"]
    .astype(str)
    .str.split(".")
    .str[0]
)

print("Human lncRNA genes:", len(human_lncs))



#analysing dGTEx


print("\ndGTEx")

# Read dGTEx metadata
dgtex_metadata = pd.read_excel(
    dgtex_metadata_file
)

# Select brain samples
dgtex_brain = dgtex_metadata[
    dgtex_metadata["SMTS"] == "Brain"
]

dgtex_brain_ids = dgtex_brain["SAMPLE_ID"].astype(str).tolist()

print("Brain samples in dGTEx metadata:", len(dgtex_brain_ids))


#read gene-level TPM GCT file
#GCT files contain two header lines before the table

dgtex_expression = pd.read_csv(
    dgtex_expression_file,
    sep="\t",
    skiprows=2
)

print("Genes in dGTEx expression matrix:", len(dgtex_expression))


#find which metadata brain samples are present in expression matrix as it is inclusive of many other tissue types

dgtex_samples_present = [
    sample
    for sample in dgtex_brain_ids
    if sample in dgtex_expression.columns
]

print(
    "dGTEx brain samples present in expression matrix:",
    len(dgtex_samples_present)
)

print("dGTEx brain sample IDs:")
for sample in dgtex_samples_present:
    print(sample)


#cleaning Ensembl IDs
#GCT gene identifier should be in the Name column

dgtex_expression["gene_id"] = (
    dgtex_expression["Name"]
    .astype(str)
    .str.split(".")
    .str[0]
)


#to calculate maximum TPM across brain samples

dgtex_expression["max_dGTEx_brain_TPM"] = (
    dgtex_expression[dgtex_samples_present]
    .max(axis=1)
)


#defining expression threshold
#using the same criterion as previous analysis: TPM > 0.1

dgtex_expression["dGTEx_brain_expressed"] = (
    dgtex_expression["max_dGTEx_brain_TPM"] > 0.1
)


#set of genes expressed in dGTEx brain

dgtex_expressed_genes = set(
    dgtex_expression.loc[
        dgtex_expression["dGTEx_brain_expressed"],
        "gene_id"
    ]
)

print(
    "All genes expressed in dGTEx brain:",
    len(dgtex_expressed_genes)
)



# ADULT GTEx


print("\nAdult GTEx")

# Read Adult GTEx sample metadata

gtex_metadata = pd.read_csv(
    gtex_metadata_file,
    sep="\t",
    low_memory=False
)

print("GTEx metadata columns:")
print(gtex_metadata.columns.tolist())


#GTEx uses SAMPID for sample ID and SMTS for broad tissue
#Select all samples classified as Brain

gtex_brain = gtex_metadata[
    gtex_metadata["SMTS"] == "Brain"
]

gtex_brain_ids = gtex_brain["SAMPID"].astype(str).tolist()

print(
    "Brain samples in Adult GTEx metadata:",
    len(gtex_brain_ids)
)


#looking at the adult GTEx
#read only the header first so the entire huge GTEx file is not loaded
gtex_header = pd.read_csv(
    gtex_expression_file,
    sep="\t",
    skiprows=2,
    nrows=0
)

#find which brain samples are present in the expression matrix
gtex_samples_present = [
    sample for sample in gtex_brain_ids
    if sample in gtex_header.columns
]

print(
    "Adult GTEx brain samples present in expression matrix:",
    len(gtex_samples_present)
)

#read the expression matrix in chunks to reduce memory use
gtex_expressed_genes = set()

for chunk in pd.read_csv(
    gtex_expression_file,
    sep="\t",
    skiprows=2,
    usecols=["Name"] + gtex_samples_present,
    chunksize=1000
):

    #remove Ensembl version suffix
    chunk["gene_id"] = (
        chunk["Name"]
        .astype(str)
        .str.split(".")
        .str[0]
    )

    # For each chunk
    brain_tpm = chunk[gtex_samples_present]

    chunk["samples_above_0.1"] = (
    brain_tpm > 0.1
    ).sum(axis=1)

    chunk["proportion_above_0.1"] = (
    chunk["samples_above_0.1"] / len(gtex_samples_present)
    )

    #find genes with TPM > 0.1 in at least one brain sample
    expressed = chunk[gtex_samples_present].max(axis=1) > 0.1

    gtex_expressed_genes.update(
        chunk.loc[expressed, "gene_id"]
    )

print(
    "All genes expressed in Adult GTEx brain:",
    len(gtex_expressed_genes)
)

#we are mapping expression back to the human lncRNA catalogue


print("\nlncRNA results")

human_lncs["dGTEx_brain_expression"] = (
    human_lncs["gene_id"].isin(dgtex_expressed_genes)
)

human_lncs["GTEx_brain_expression"] = (
    human_lncs["gene_id"].isin(gtex_expressed_genes)
)


#showing expression in either developmental OR adult GTEx

human_lncs["GTEx_any_brain_expression"] = (
    human_lncs["dGTEx_brain_expression"] |
    human_lncs["GTEx_brain_expression"]
)



dgtex_lnc_count = human_lncs[
    "dGTEx_brain_expression"
].sum()

gtex_lnc_count = human_lncs[
    "GTEx_brain_expression"
].sum()

both_count = (
    human_lncs["dGTEx_brain_expression"] &
    human_lncs["GTEx_brain_expression"]
).sum()

either_count = human_lncs[
    "GTEx_any_brain_expression"
].sum()


print("\nHuman lncRNA expression results:")

print(
    "lncRNAs expressed in dGTEx brain:",
    dgtex_lnc_count
)

print(
    "lncRNAs expressed in Adult GTEx brain:",
    gtex_lnc_count
)

print(
    "lncRNAs expressed in both:",
    both_count
)

print(
    "lncRNAs expressed in either dGTEx or Adult GTEx:",
    either_count
)

#creating a single catalogue of hum lncRNAs expressed in brain in GTEx and dGTEx

human_lncs.to_csv(
    "Human_lncRNAs_dGTEx_and_GTEx_brain_expression_status.tsv",
    sep="\t",
    index=False
)


#showing the lncRNAs only expressed in either dataset

expressed_lncs = human_lncs[
    human_lncs["GTEx_any_brain_expression"]
].copy()

expressed_lncs.to_csv(
    "Human_lncRNAs_expressed_in_dGTEx_or_GTEx_brain.tsv",
    sep="\t",
    index=False
)

print(
    "\nOutput: "
    "Human_lncRNAs_expressed_in_dGTEx_or_GTEx_brain.tsv"
)