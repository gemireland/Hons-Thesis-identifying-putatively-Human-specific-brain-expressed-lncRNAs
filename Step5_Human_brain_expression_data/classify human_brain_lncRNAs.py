#!/usr/bin/env python3

#this scripts is being written to take the updated expanded chimp and gorilla conserved lncRNA files (from the synteny analysis on the novel lncRNAs from newly assemled transcriptomes)
#we are also looking at the original conserved mouse lncRNAs as these were not reanalysed as there was a decent number of conserved annotated human lncRNAs to mouse

import pandas as pd
#to start these are my input files that I will be comparing

human_brain_expressed = pd.read_csv("Human_brain_expressed_lncRNA_catalogue.tsv", sep="\t")
updated_chimp_conserved = pd.read_csv("UpdatedConserved_Chimp_lncRNAs.tsv", sep="\t")
updated_gorilla_conserved = pd.read_csv("UpdatedConserved_Gorilla_lncRNAs.tsv", sep="\t")
mouse_conserved = pd.read_csv("Conserved_Mouse_lncRNAs_Cleaned.tsv", sep="\t")

#once again some of my files have version numbers while other don't so we need to standardise it by stripping 
#the human one alreday has the cleaned gene id
human_brain_expressed["gene_id"] = human_brain_expressed["gene_id"].astype(str).str.strip()
updated_chimp_conserved["human_gene_id"] = updated_chimp_conserved["human_gene_id"].str.split(".").str[0]
updated_gorilla_conserved["human_gene_id"] = updated_gorilla_conserved["human_gene_id"].str.split(".").str[0]
mouse_conserved["human_gene_id"] = mouse_conserved["human_gene_id"].str.split(".").str[0]

#now we are creating sets where the human gene IDs are conserved
chimp_ids = set(updated_chimp_conserved["human_gene_id"].dropna())
gorilla_ids = set(updated_gorilla_conserved["human_gene_id"].dropna())
mouse_ids = set(mouse_conserved["human_gene_id"].dropna())

#seeing if our conserved IDs are in each species
human_brain_expressed["Chimp_conserved"] = human_brain_expressed["gene_id"].isin(chimp_ids)
human_brain_expressed["Gorilla_conserved"] = human_brain_expressed["gene_id"].isin(gorilla_ids)
human_brain_expressed["Mouse_conserved"] = human_brain_expressed["gene_id"].isin(mouse_ids)

#classifying conservation categories across the species
def classify(speciesconserved):
    if speciesconserved["Mouse_conserved"]:
        return "Mammalian-conserved"
    elif speciesconserved["Chimp_conserved"] or speciesconserved["Gorilla_conserved"]:
        return "Primate-conserved"
    else: 
        return "Putatively human brain specific"


human_brain_expressed["Conservation_category"] = human_brain_expressed.apply(classify, axis=1)

#creating a file where all conservation statuses are kept for human
human_brain_expressed.to_csv("Human_brain_lncRNAs_all_conservation.tsv", sep="\t", index=False)

#save each category independently
mammalian_conserved = human_brain_expressed[human_brain_expressed["Conservation_category"] == "Mammalian-conserved"]
primate_conserved = human_brain_expressed[human_brain_expressed["Conservation_category"] == "Primate-conserved"]
human_brain_specific = human_brain_expressed[human_brain_expressed["Conservation_category"] == "Putatively human brain specific"]

mammalian_conserved.to_csv("Human_brain_mammalian_connserved.tsv", sep="\t", index=False)
primate_conserved.to_csv("Human_brain_primate_conserved.tsv", sep="\t", index=False)
human_brain_expressed.to_csv("Human_brain_putatively_human_specific.tsv", sep="\t", index=False)

#printing the results to our terminal
#quick check to see the 18 859 human brain expressed lncs
print("total hum brain expressed lncs:expecting 18 859", len(human_brain_expressed))

print("Conserved hum brain in chimp:", human_brain_expressed["Chimp_conserved"].sum())
print("Conserved hum brain in gorilla:", human_brain_expressed["Gorilla_conserved"].sum())
print("Conserved hum brain in mouse:", human_brain_expressed["Mouse_conserved"].sum())

print("Mammalian-conserved:", len(mammalian_conserved))
print("Primate-conserved:", len(primate_conserved))
print("Putatively-human-brain-specific:", len(human_brain_specific))
