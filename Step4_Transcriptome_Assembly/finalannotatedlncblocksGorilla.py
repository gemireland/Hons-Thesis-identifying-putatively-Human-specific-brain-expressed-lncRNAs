#!/usr/bin/env python3

#this scipt is created to complie a list of the final annotated species lncRNA blocks
#it is combining the blocks created from the initial Ensembl from the species as well as from the merged stringtie lncRNA blocks
import csv

#input files
ensembl_file = "Gorilla_lncRNA_blocks.tsv"
stringtie_file = "Gorilla_final_novel_lncRNA_blocks.tsv"

#we want to read in the lncRNA coordinates from the ensembl file
ensembl = set()

with open(ensembl_file) as infile:

    reader = csv.DictReader(infile, delimiter="\t")

    for row in reader:
        chrom = row["chromosome"].replace("chr", "")
        start = int(row["start"])
        end = int(row["end"])

        ensembl.add((chrom, start, end))

#we want to read in the lncRNA coordinates from the stringtie file
stringtie = set()

with open(stringtie_file) as infile:

    reader = csv.DictReader(infile, delimiter="\t")
#remember stringtie headers were capitalised
    for row in reader:
        chrom = row["Chromosome"].replace("chr", "")
        start = int(row["Start"])
        end = int(row["End"])

        stringtie.add((chrom, start, end))

#to look at the differences between the files
same = ensembl & stringtie

print("Annotated ensembl lncRNA blocks:", len(ensembl))
print("Annotated stringtie lncRNA blocks:", len(stringtie))
print("LncRNAs with same coordinates in both files:", len(same))
