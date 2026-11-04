#!/usr/bin/env python3

#this is the script to compare the final retained conserved lncRNAs across human to chimp/gorilla from the original synteny and from the synteny analysis on StringTie assembled gtf
import csv

#these are the cleaned lnc files after the final round of orientation filtering (as described in supplementary material in thesis)

original_conserved_lncs = "Conserved_Gorilla_lncRNAs_Cleaned.tsv"
stringtie_conserved_lncs = "StringTieConserved_Gorilla_lncRNAs_Cleaned.tsv"

def get_ids(file):
    with open(file, "r") as infile:

        reader = csv.DictReader(infile, delimiter="\t")
        return {row["human_gene_id"] for row in reader}

#extracting the gene_ids from the original and stringtie files
original_ids = get_ids(original_conserved_lncs)
stringtie_ids = get_ids(stringtie_conserved_lncs)

#now to compare the sets of original vs newly identified 
both = original_ids & stringtie_ids

#take everything that is in original and remove what is also found in stringtie
original_only = original_ids - stringtie_ids

#take everything from stringtie and remove anything also in original
stringtie_only = stringtie_ids - original_ids

#the pipe function essentially calculates the original as: original + stringtie - both
total = original_ids | stringtie_ids

print("conserved from orig annotations", len(original_ids))
print("conserved to stringtie", len(stringtie_ids))
print("in both stringtie and orig", len(both))
print("original only:", len(original_only))
print("stringtie only:", len(stringtie_only))
print("final Gorilla lncs that are uniquely conserved", len(total))


