#!/usr/bin/env python3 
import csv
#this script is to gorup all transcript isoforms into gene-level blcoks to use in syntenic analysis

#these are the files being used/created

#my input file is the final candidate novel lncRNAs that were identified from my merged StringTie gtf after cpc2 and feelnc were run
#they are all over 200nt, been classified as lncRNAs by FEELnc and as non-coding by cpc2

input_file = "/scratch/irlgem001/Step4-RNASeq/Chimp/Chimp_final_novel_lncRNAs.tsv"

#the output file will be one genomic block representing the interval of all transcript isoforms of a gene
output_file = "/scratch/irlgem001/Step4-RNASeq/Chimp/Chimp_final_novel_lncRNA_blocks.tsv"

#to read in the final lncRNA files that has multiple transcript isoforms per gene, we are grouping them according to gene_id

gene_blocks = {}

with open(input_file) as infile:
    reader = csv.DictReader(infile, delimiter="\t")

    for row in reader:

        chromosome = row["Chromosome"]
        start = int(row["Start"])
        end = int(row["End"])
        strand = row["Strand"]
        gene_id = row["Gene_ID"]

        #if this is the first transcript encountered for the gene create a new gene block 
        if gene_id not in gene_blocks:

            gene_blocks[gene_id] = {
                "chromosome": chromosome,
                "start": start,
                "end": end,
                "strand": strand,
                "transcript_count": 1
            }

        #if another transcript belongs to the same gene, allow adding to the gene block to include all isoform transcripts of the same gene
        else:
            gene_blocks[gene_id]["start"] = min(gene_blocks[gene_id]["start"], start)

            gene_blocks[gene_id]["end"] = min(gene_blocks[gene_id]["end"], end)

            gene_blocks[gene_id]["transcript_count"] += 1


#writing the lncRNA gene blocks where each row is one lncRNA gene not a transcript isoform

with open(output_file, "w") as outfile:

    outfile.write(
        "Chromosome\t"
        "Start\t"
        "End\t"
        "Strand\t"
        "Gene_ID\t"
        "Transcript_Count\n"
    )

    for gene_id, info in gene_blocks.items():

        outfile.write(
            f"{info['chromosome']}\t"
            f"{info['start']}\t"
            f"{info['end']}\t"
            f"{info['strand']}\t"
            f"{gene_id}\t"
            f"{info['transcript_count']}\n"
        )

#printing results to terminal

print("lncRNA gene block analysis")

print("candidate lncRNA transcripts:", sum(info["transcript_count"] for info in gene_blocks.values()))

print("lncRNA gene blocks:", len(gene_blocks))

print("Output file:", output_file)
                        