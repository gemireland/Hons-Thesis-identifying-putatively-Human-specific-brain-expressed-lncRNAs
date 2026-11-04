#!/usr/bin/env python3

import csv
import re

#this script is being used to extract the novel intergenic IDs from the created tsv file and comparing them to the StringTie merged gtf to get the coordinates of 
#the novel intergenic IDs so we can get their sequnces (making them fasta files) to then run cpc and feelnc to identify protein-coding ability

#the input files are
#change this for Chimpanzee !!!
novel_tsvfile = "Gorilla_NOVEL_INTERGENIC.tsv"
stringtie_gtf = "Gorilla_AllStringTiegtfs_merged_FPKM0.5.gtf"

#output file that will be written 
output_gtf = "Gorilla_NOVEL_INTERGENIC.gtf"

#to read in the transcript ids of the novel intergenic transcripts
novel_transcript_ids = set()

with open(novel_tsvfile, "r") as infile:
    reader = csv.DictReader(infile, delimiter ="\t")

    for row in reader:
        transcript_id = row["transcript_id"]

        #only keep trasncripts classified as novel intergenic or "u" by gffcompare as class code
        if row["class_code"] == "u":
            novel_transcript_ids.add(transcript_id)

#sanity check that we are picking up the correct no. novel intergenic transcripts
print("novel intergenic transcript IDs", len(novel_transcript_ids))

#to extract these novel intergenic transcript ids from stringtie gtfs (create a set)
found_transcript_ids = set()
written_lines = 0

with open(stringtie_gtf, "r") as infile, open(output_gtf, "w") as outfile:
    for line in infile:
        if line.startswith("#"):
            continue

        #remember that rstrip removes extra whitespace and ensures newlines don't start at the end of the line above (I noticed this didn't run initially with thsi step)
        fields = line.rstrip("\n").split("\t")

        #this is because most gtf files have 9 columns
        if len(fields) < 9:    
            continue   #skipped 

        attributes = fields[8]

        #extract the transcript_id from the attributes column of gtf
        #re means regular expression and is used to search for a pattern in a string
        match = re.search(r'transcript_id "([^"]+)"', attributes)

        #this means that if the transcript_id is found in the attributes column of the gtf, we will extract it and check if it is in the novel_transcript_ids set
        if match:
            transcript_id = match.group(1)

            #write trasncript lines that are in novel_transcript_ids to the novel output file
            if transcript_id in novel_transcript_ids:
                outfile.write(line)
                found_transcript_ids.add(transcript_id)
                written_lines += 1

#creating a summary of the results
print("transcript IDs novel", len(novel_transcript_ids))
print("novel transcript IDs found in GTF:", len(found_transcript_ids))
print("GTF lines written", written_lines)
print("saved result to", output_gtf)

#calculting the ones missed 
missing = novel_transcript_ids - found_transcript_ids

print("no. transcript ids not found", len(missing))



            

