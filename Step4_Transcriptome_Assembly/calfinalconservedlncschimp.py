#!/usr/bin/env python3

#this is the script to compare the final retained conserved lncRNAs across human to chimp/gorilla from the original synteny and from the synteny analysis on StringTie assembled gtf
import csv

#these are the cleaned chimp lnc files after the final round of orientation filtering
original_conserved_lncs = "Conserved_Chimp_lncRNAs_Cleaned.tsv"
stringtie_conserved_lncs = "StringTieConserved_Chimp_lncRNAs_Cleaned.tsv"

#defining our function called get_ids
def get_ids(file):
    with open(file, "r") as infile:

        reader = csv.DictReader(infile, delimiter="\t")
        return {row["human_gene_id"] for row in reader}

#extracting the gene_ids from the original and stringtie files
original_ids = get_ids(original_conserved_lncs)
stringtie_ids = get_ids(stringtie_conserved_lncs)

#now to compare the sets
both = original_ids & stringtie_ids

#take everything that is in original and remove what is also in stringtie
original_only = original_ids - stringtie_ids

#take everything from stringtie and remove anything also in original
stringtie_only = stringtie_ids - original_ids

#the pipe function essentially calculates the original as: original + stringtie - both
total = original_ids | stringtie_ids

print("conserved in original", len(original_ids))
print("conserved in StringTie", len(stringtie_ids))
print("present in StringTie and orig", len(both))
print("original only", len(original_only))
print("stringtie only", len(stringtie_only))
print("final Chimp lncs unique", len(total))



