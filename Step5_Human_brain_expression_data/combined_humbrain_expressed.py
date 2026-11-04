#!/usr/bin/env python3

#this script is being written to take the identified lncRNAs as expressed in BrainSpan, dGTEx and GTEx to put out a catalogue of lncRNAs that are expressed during brain development
#this catalogue doesn't tell us which lncRNAs are expressed at certain times or brain reions yet but rather tells us the ones that are identified in development
#pandas to help understand data in the table format
import pandas as pd

#this the the genocde annotated human lncRNA blocks, we are using it to refer back to the ID 
human_lnc_blocks = pd.read_csv("Human_lncRNA_blocks.tsv", sep="\t")

#our input files are the expressed in brain human lncRNA catalogues
brainspan_expressed_lncs = pd.read_csv("/scratch/irlgem001/Step4-RNASeq/HumanExpression/BrainSpan/genes_matrix_csv/FinalHumanlncRNA_BrainSpan_expressed.tsv", sep="\t")

#for some reason I just kept both dGTEx and GTEx in the "dGTEx" folder under HumanExpression
dgtex_expressed_lncs = pd.read_csv("/scratch/irlgem001/Step4-RNASeq/HumanExpression/dGTEx/Human_lncRNAs_expressed_in_dGTEx_brain.tsv", sep="\t")

gtex_expressed_lncs = pd.read_csv("/scratch/irlgem001/Step4-RNASeq/HumanExpression/dGTEx/Human_lncRNAs_expressed_in_GTEx_brain.tsv", sep="\t")

#to make the ids comparable as some have a version number like .6 or .2 we need to strip the suffixes
human_lnc_blocks["gene_id"] = human_lnc_blocks["lncrna_gene_id"].str.split(".").str[0]
brainspan_expressed_lncs["gene_id"] = brainspan_expressed_lncs["lncrna_gene_id"].str.split(".").str[0]
dgtex_expressed_lncs["gene_id"] = dgtex_expressed_lncs["lncrna_gene_id"].str.split(".").str[0]
gtex_expressed_lncs["gene_id"] = gtex_expressed_lncs["lncrna_gene_id"].str.split(".").str[0]

#keep all the human lncRNAs that are expressed in atleast 1 dataset
expressed_ids = (set(brainspan_expressed_lncs["gene_id"]) | set(dgtex_expressed_lncs["gene_id"]) | set(gtex_expressed_lncs["gene_id"]))
brain_expressed_lncs = human_lnc_blocks[human_lnc_blocks["gene_id"].isin(expressed_ids)].copy()

print("BrainSpan", len(set(brainspan_expressed_lncs["gene_id"])))
print("dGTEx", len(set(dgtex_expressed_lncs["gene_id"])))
print("GTEx", len(set(gtex_expressed_lncs["gene_id"])))
print("Unique brain-expressed lncRNAs:", len(brain_expressed_lncs))

brain_expressed_lncs.to_csv(
    "Human_brain_expressed_lncRNA_catalogue.tsv",
    sep="\t",
    index=False
)

#to see the overlap between the sets we will make sets of gene IDs for each dataset
brainspan_ids = set(brainspan_expressed_lncs["gene_id"])
dgtex_ids = set(dgtex_expressed_lncs["gene_id"])
gtex_ids = set(gtex_expressed_lncs["gene_id"])


all_three = brainspan_ids & dgtex_ids & gtex_ids

brainspan_dgtex_only = (brainspan_ids & dgtex_ids) - gtex_ids
brainspan_gtex_only = (brainspan_ids & gtex_ids) - dgtex_ids
dgtex_gtex_only = (dgtex_ids & gtex_ids) - brainspan_ids

brainspan_only = brainspan_ids - dgtex_ids - gtex_ids
dgtex_only = dgtex_ids - brainspan_ids - gtex_ids
gtex_only = gtex_ids - brainspan_ids - dgtex_ids

print("\noverlap categories")
print("BrainSpan only", len(brainspan_only))
print("dGTEx only", len(dgtex_only))
print("GTEx only", len(gtex_only))
print("BrainSpan and dGTEx only", len(brainspan_dgtex_only))
print("BrainSpan and GTEx only", len(brainspan_gtex_only))
print("dGTEx and GTEx only", len(dgtex_gtex_only))
print("all three datasets", len(all_three))
