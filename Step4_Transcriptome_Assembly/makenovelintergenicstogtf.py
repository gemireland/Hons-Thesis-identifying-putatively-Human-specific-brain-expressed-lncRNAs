#!/usr/bin/env python3

import csv
import re

#this script is being used to extract the novel intergenic IDs from the created tsv file and comparing them to the StringTie merged gtf to get the coordinates of 
#the novel intergenic IDs so we can get their sequnces (making them fasta files) to then run cpc and feelnc to identify protein-coding ability

#the input files are:
novel_tsvfile = "Gorilla_NOVEL_INTERGENIC.tsv"
stringtie_gtf = "Gorilla_AllStringTiegtfs_merged_FPKM0.5.gtf"

#output file
output_gtf = "Gorilla_NOVEL_INTERGENIC.gtf"

#to read in the transcript ids of the novel intergenic transcripts
novel_transcript_ids = set()

with open(novel_tsvfile, "r") as infile:
    reader = csv.DictReader(infile, delimiter ="\t")

    for row in reader:
        transcript_id = row["transcript_id"]

        #only keep trasncripts classified as noverl intergenic or "u"
        if row["class_code"] == "u":
            novel_transcript_ids.add(transcript_id)

#sanity check 
print("Novel intergenic transcript IDs:", len(novel_transcript_ids))

#to extract these novel intergenic transcript ids from stringtie gtfs
found_transcript_ids = set()
written_lines = 0

with open(stringtie_gtf, "r") as infile, open(output_gtf, "w") as outfile:
    for line in infile:
        if line.startswith("#"):
            continue

        #remember that rstrip removes extra whitespace and ensures newlines don't start at the end of the line above
        fields = line.rstrip("\n").split("\t")

        #this is because most gtf files have 9 columns
        if len(fields) < 9:    
            continue

        attributes = fields[8]

        #extract the transcript_id from the attributes column of gtf
        #re means regular expression and is used to search for a pattern in a string
        match = re.search(r'transcript_id "([^"]+)"', attributes)

        #this means that if the transcript_id is found in the attributes column of the gtf, we will extract it and check if it is in the novel_transcript_ids set
        if match:
            transcript_id = match.group(1)

            #write trasncript/exon lines that are in novel_transcript_ids to the novel output file
            if transcript_id in novel_transcript_ids:
                outfile.write(line)
                found_transcript_ids.add(transcript_id)
                written_lines += 1

#creating a summary of the results
print("Novel transcript IDs requested:", len(novel_transcript_ids))
print("Novel transcript IDs found in GTF:", len(found_transcript_ids))
print("GTF lines written:", written_lines)
print("Output saved to:", output_gtf)

missing = novel_transcript_ids - found_transcript_ids

print("Transcript IDs not found:", len(missing))

if missing:
    print("First 10 missing transcript IDs:")
    for transcript_id in sorted(missing)[:10]:
        print(transcript_id)

            

