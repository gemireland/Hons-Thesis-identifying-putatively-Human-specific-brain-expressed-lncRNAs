#!/usr/bin/env python3

import pandas as pd

#path to gene-level dGTEx expression file
expression_file = "dGTEx_Analysis_2025-07-21_v1_RNASeqV2.4.3_gene_tpm.gct"

#our 4 brain samples identified from the dGTEx sample metadata
brain_samples = [
    "DGTEX001-2BZ-01-SM-OLAQB",
    "DGTEX002-2BZ-01-SM-OLMRP",
    "DGTEX031-2BZ-01-SM-OLASE",
    "DGTEX047-2BZ-01-SM-OLATP"
]

#reading the expression file knowing that the first 2 rows have headers of metadata
expression = pd.read_csv(
    expression_file,
    sep="\t",
    skiprows=2
)

print("Total no. genes in dGTEx:", len(expression))
print("Total no. columns:", len(expression.columns))

#determining if the 4 brain samples as liste din metadata are present in the expression file
print("\nChecking brain samples:")

for sample in brain_samples:
    print(sample, ":", sample in expression.columns)

#stop run if one of the sample IDs cannot be found
missing_samples = [
    sample for sample in brain_samples
    if sample not in expression.columns
]

if missing_samples:
    raise ValueError(
        f"Brain samples not found: {missing_samples}"
    )

#retain the gene information and four brain samples
brain_expression = expression[
    ["Name", "Description"] + brain_samples
].copy()

print("\nFinal brain-only dataset:")
print("Genes:", len(brain_expression))
print("Brain samples:", len(brain_samples))

#create a brain gene-level expression file
brain_expression.to_csv(
    "dGTEx_onlybrain_genes_TPM.tsv",
    sep="\t",
    index=False
)

print("\nSaved: dGTEx_onlybrain_genes_TPM.tsv")