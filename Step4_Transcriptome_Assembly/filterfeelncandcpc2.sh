#!/usr/bin/env python3

#I am wrirting this script to filter the FEELnc and CPC2 results for lncRNAs that are predicted by both potential calculators
#I will use the GTF file from FEELnc and the txt file from CPC2 to create a new file with only the lncRNAs that are present in both files

#re is needed to parse the GTF file and extract the gene_id from the attributes field
import re
import csv

#these are my input files
initial_novel_intergenicgtf = "/scratch/irlgem001/Step4-RNASeq/Gorilla/Gorilla_NOVEL_INTERGENIC.gtf"

FEELNC_GTF = "/scratch/irlgem001/Step4-RNASeq/Gorilla/Gorilla_FEELnc_results/codpot/Gorilla_FEELnc.codpot.lncRNA.gtf"

CPC2_FILE = "/scratch/irlgem001/Step4-RNASeq/Chimp/cpc2gorilla_results.txt"

#these is my output file
newlnc_cpc2_feelnc_final = "/scratch/irlgem001/Step4-RNASeq/Gorilla/Gorilla_final_novel_lncRNAs.tsv"

## I will read the files of the original transcripst and calculate the length 
transcripts = {}

with open(initial_novel_intergenicgtf) as infile:

    for line in infile:
        if line.startswith("#"):
            continue

        fields = line.rstrip().split("\t")

        #only keep correct gtf lines that have 9 fields and that are exons (not introns!!!!)
        if len(fields) < 9 or fields[2] != "exon":
            continue

        chromosome = fields[0]
        start = int(fields[3])
        end = int(fields[4])
        strand = fields[6]
        attributes = fields[8]

            #we are extracting the gene_id and transcript_id from column 9 (attributes)
        gene_match = re.search(r'gene_id "([^"]+)"', attributes)
        transcript_match = re.search(r'transcript_id "([^"]+)"', attributes)

        #this means skip the row if either gene_id or transcript_id is not found in the attributes
        if not gene_match or not transcript_match:
            continue

        gene_id = gene_match.group(1)
        transcript_id = transcript_match.group(1)

        #if this is the first exon the computer has encounted for the transcript, create a new dictionary entry for this transcript and store its exons in a list. If the transcript has already been seen, append the new exon to the list of exons for that transcript.
        if transcript_id not in transcripts:
            transcripts[transcript_id] = {
                "chromosome": chromosome,
                "strand": strand,
                "gene_id": gene_id,
                "exons": []
            }

        #if the exon has already been seen for that transcript append it to the list
        #we want to store the start and end coordinates for the exons and all exons belonging to the same transcript are grouped together
        transcripts[transcript_id]["exons"].append((start, end))

#this is a filter to ensure all transcripts are longer than 200 nucleotides (this is the minimum length of a lncRNA by definition)
#we are first calculating the length of each transcript by summating the lengths of all its exons
#remember this is the SPLICED transcript length
long_transcripts = {}

for transcript_id, info in transcripts.items():

    #you have to +1 for gtf coordinates
    length = sum(end - start + 1 for start, end in info["exons"])

    #from here we only retain transcripts longer than 200nt we are also creating the genomic interval (min vs max) 
    if length >= 200:

        genomic_start = min(exon[0] for exon in info["exons"])

        genomic_end = max(exon[1] for exon in info["exons"])

        #we want to store the transcript and genomic info
        long_transcripts[transcript_id] = {
            "chromosome": info["chromosome"],
            "start": genomic_start,
            "end": genomic_end,
            "strand": info["strand"],
            "gene_id": info["gene_id"],
            "length": length
        }

#we will now find transcripts that are defined as lncRNAs by FEELnc
#we will use a set as each transcript may occur many times in the gtf for each of its exons but we only want the transcript id 
feelnc = set()

with open(FEELNC_GTF) as infile:

    for line in infile:

        #find the transcript_id in the attributes column
        match = re.search(  
            r'transcript_id "([^"]+)"', 
            line
        )

        if match:
            transcript_id = match.group(1)

            #we would like to add the id to the feelnc non-coding set
            feelnc.add(transcript_id)

#this section of the script looks to identify the transcripts that were predicted by cpc2 as non-coding
#we are storing the non-coding classified transcripts into another set

cpc2 = set()

with open(CPC2_FILE) as infile:

    for line in infile:

        if line.startswith("#"): #ignore headers or comment lines
            continue

        fields = line.rstrip().split("\t")

        #all the cpc2 rows should have 8 columns
        if len(fields) < 8:
            continue

        transcript_id = fields[0]
        #the cpc2 results has lowercase for coding status 
        classification = fields[7].strip().lower()

        #retain only the non-coding transcripts
        if classification == "noncoding":
            cpc2.add(transcript_id)

#this step will identify the candidate novel lncRNAs in species
#the criteria for selection are: found from novel intergenic transcript "u" set, spliced length of >= 200nt, feelnc classified as lncRNA
#and classified as noncoding by CPC2
#"&" finds the intersecion between the sets

final_transcripts = ( set(long_transcripts) & feelnc & cpc2)

#this section is to compile a table of the final lncRNA candidates from cpc2, feelnc and the 200 nt filtering step
#remember length is the sum of the exons so it doesn't inlcude introns!!!
#also remember that these candidates are all >200nt as classified as non-coding by BOTH tools indep.

with open(newlnc_cpc2_feelnc_final, "w") as outfile:

    #writing out the column headings
    outfile.write(
        "Chromosome\t"
        "Start\t"
        "End\t"
        "Strand\t"
        "Gene_ID\t"
        "Transcript_ID\t"
        "Length\n"
    )

    #to store information about the final candidate transcripts in the output file
    for transcript_id in sorted(final_transcripts):

        info = long_transcripts[transcript_id]

        outfile.write(
            f"{info['chromosome']}\t"
            f"{info['start']}\t"
            f"{info['end']}\t"
            f"{info['strand']}\t"
            f"{info['gene_id']}\t"
            f"{transcript_id}\t"
            f"{info['length']}\n")

#priniting the numbers of lncRNAs retained or stored at each stage
print("novel lncRNA analysis for Gorilla")
print("novel intergenic transcripts orig", len(transcripts))
print("nov. intergen transcrip >= 200 nt:", len(long_transcripts))
print("transcripts defined as lncRNAs by FEELnc", len(feelnc))
print("non-coding classified transcripts by cpc2", len(cpc2))
print("final transcripts >= 200nt and non-coding in both", len(final_transcripts))
print("final output " , newlnc_cpc2_feelnc_final)


