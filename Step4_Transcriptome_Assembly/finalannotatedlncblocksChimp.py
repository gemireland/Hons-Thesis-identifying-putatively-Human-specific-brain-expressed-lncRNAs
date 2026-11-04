#!/usr/bin/env python3

#this scipt is created to complie a list of the final annotated species lncRNA blocks (expanded lncRNAs)
#it is combining the blocks created from the initial Ensembl from the species as well as from the merged stringtie lncRNA blocks
import csv

#input files
ensembl_file = "Chimp_lncRNA_blocks.tsv"
stringtie_file = "Chimp_final_novel_lncRNA_blocks.tsv"

#we want to process the lncRNA coords from the ensembl file (extract the coords and put them in a set)
ensembl = set()

with open(ensembl_file) as infile:

    reader = csv.DictReader(infile, delimiter="\t")

    for row in reader:
        chrom = row["chromosome"].replace("chr", "")   #removing chr
        start = int(row["start"])
        end = int(row["end"])

        ensembl.add((chrom, start, end))

#we want to read in the lncRNA coordinates from the stringtie file (also store in a set)
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

print("lncRNA blocks annotated in Ensembl", len(ensembl))
print("Stringtie lncRNA blocks annotated", len(stringtie))
print("lncRNAs with exact coords in both files", len(same))
