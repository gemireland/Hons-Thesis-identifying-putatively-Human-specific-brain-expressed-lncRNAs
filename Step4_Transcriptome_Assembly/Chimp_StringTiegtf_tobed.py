#!/usr/bin/env python3

#this script serves as a tool to convert the gtf merged files from STringTie to .bed12 files
#we do this because the merged gtf files include trasncript strcutures that tell us where in introns and exons are
#however, bedtools getfasta just looks at BED coordinates and would extract everything from trancript start to end including introns (Which we don't want as we want the mature spliced RNA transcript)
#this means taht we will also need to use bedtools getfasta -split -s which will exctract only the exon structures and joins them together (removing the introns)

#This is to import modules we need to convert teh gtfs into bed files
import re

#defaultdict automatically creates a default value for a new key, whereas a normal dict requires the key to be created first.
from collections import defaultdict

#defining that the input gtf is the  merged StringTie GTF file (not the orginal GTF from Ensembl)
input_gtf = "Chimp_AllStringTiegtfs_merged_FPKM0.5.gtf"

#defining the output BED12 file
#We use BED12 as it stores the individual exon positions within each transcript
output_bed = "Chimp_StringTie_transcripts.bed12"

#we want to all exon coordinates belonging to each transcript
#this is because a transcript contains many exons and introns so we want the start adn end of every exon contained within a transcript
#for example: MSTRG.2.1 -> [(4306, 5396), (10268, 11782)]
transcripts = defaultdict(list)

#creating a dictionary to store the strand for each transcript
strands = {}

#creating a dictionary to store the chromosome for each transcript
chromosomes = {}

#we need to begin to read in the merged StringTie GTF
with open(input_gtf, "r") as infile:

    for line in infile:

        #firstly ignore the comment/header lines
        if line.startswith("#"):
            continue

        #we want to split the GTF row into its 9 tab-separated columns
        fields = line.rstrip().split("\t")

        #prompt the computer to ignore rows that do not contain all 9 GTF columns
        if len(fields) < 9:
            continue

        #we only look to use exon rows because these give the individual exon coordinates needed to exclude introns
        if fields[2] != "exon":
            continue
 
        #take teh value in column 0 and give teh variable name as chromosome
        chromosome = fields[0]

        #GTF coordinates are 1-based whereas BED coordinates start at 0
    #so the start of the exon is in column 4 in the gtf file but we subtract 1 as BED coordinates are gtf coordinates -1 (eg: 4089 GTF -> 4088 start in bed) (only do this for start not end)
    #we also need to make the text an integer for the computer to understand it as a number 
        start = int(fields[3]) - 1

        #The end coordinate can remain unchanged (not -1 as start already subtracted 1), it must still be changed from text into a number
        end = int(fields[4])

        #find the transcript strand
        strand = fields[6]

        #column 9 contains the gene and transcript information
        #however it is a list of things we will need to filter through and extract things form
        attributes = fields[8]

        #so we will extract the transcript ID from the attributes column
        #look for the text inside of the transcript id's " ", capture what is inside the barckets (), match anything that is not a ", match these caharcters + and then stop matching when you reach another "
        match = re.search(r'transcript_id "([^"]+)"', attributes)

        # Skip the row if there is no transcript ID
        if not match:
            continue
    
    #this will take teh transcript ID found by the previous re.search and save it in a variable called transcript_id 
        transcript_id = match.group(1)

        #add the exon coordinates to the appropriate transcript
        transcripts[transcript_id].append((start, end))

        #store the strand and chromosome for the transcript
        strands[transcript_id] = strand
        chromosomes[transcript_id] = chromosome


#now this is going to write all of the above transcript information to a BED12 file
with open(output_bed, "w") as outfile:

    for transcript_id, exons in transcripts.items():

        #we will start by sorting exons by genomic coordinate
        exons.sort()

        chrom = chromosomes[transcript_id]
        strand = strands[transcript_id]

        #this is to get the start and end of the complete transcript region (the beginning if the first exon of the transcript and the end of the ;ast exon in the transcript)
        transcript_start = exons[0][0]
        transcript_end = exons[-1][1]    
    #for the start you take the first exon pair [0] and take its start coord (4067, 4987) therefore teh stored value for the trsnacript beginning is 4067. -1 means the last exon pair and you take teh end coordinate


        #this is to store the size and relative starting position of each exon
        block_sizes = []
        block_starts = []

        for exon_start, exon_end in exons:

            #calculate exon length
            block_sizes.append(str(exon_end - exon_start))

            #calculate exon start relative to the start of the transcript
            block_starts.append(str(exon_start - transcript_start))

        #find the number of exons in the transcript
        block_count = len(exons)

        #writing the headers of the 12 BED columns
        outfile.write(
            f"{chrom}\t"                       # chromosome
            f"{transcript_start}\t"           # transcript start
            f"{transcript_end}\t"             # transcript end
            f"{transcript_id}\t"              # transcript ID
            f"0\t"                            # BED score
            f"{strand}\t"                     # transcript strand
            f"{transcript_start}\t"           # thickStart
            f"{transcript_end}\t"             # thickEnd (the way BED12 names the start and end of a feature, howver we want teh full transcript so we don't really use thick start and end)
            f"0\t"                            # display colour
            f"{block_count}\t"                # number of exons
            f"{','.join(block_sizes)}\t"      # exon lengths
            f"{','.join(block_starts)}\n"     # exon relative start positions
        )

# Print a message when the script has finished
print(f"BED12 written to: {output_bed}")
print(f"Transcripts written: {len(transcripts)}")
