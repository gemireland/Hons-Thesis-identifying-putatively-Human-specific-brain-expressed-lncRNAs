#!/usr/bin/env python3

#this script serves as a tool to convert the gtf merged files from STringTie to BED files
#we do this because the merged gtf files include trasncript strcutures that tell us where in introns and exons are
#however bedtools getfasta just looks at BED coordinates and would extract everything from trancript start to end including introns (Which we don't want as we want the mature spliced RNA transcript)
#This is to import modules we need to convert teh gtfs into bed files
import re

#defaultdict automatically creates a default value for a new key, whereas a normal dict requires the key to be created first
from collections import defaultdict

#defining that the input gtf is the  merged StringTie GTF file (not the orginal GTF from Ensembl)
input_gtf = "Gorilla_AllStringTiegtfs_merged_FPKM0.5.gtf"

#defining the output BED file
#We use BED12 as it stores the individual exon positions within each transcript
output_bed = "Gorilla_StringTie_transcripts.bed12"

#we want to all exon coordinates belonging to each transcript
#this is because a transcript contains many exons and introns so we want the start adn end of every exon contained within a transcript
#eg: MSTRG.2.1 has 2 sets of coords[(4306, 5396), (10268, 11782)]
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

        #we want to split the GTF row into its tab-separated columns (9 columns)
        fields = line.rstrip().split("\t")

        #tell python to ignore rows that do not contain all 9 GTF columns
        if len(fields) < 9:
            continue

        #we only look to use exon rows because these give the individual exon coordinates needed to exclude introns
        if fields[2] != "exon":
            continue
 
        #take teh value in column 0 and give the variable name as chromosome
        chromosome = fields[0]

        #GTF coordinates use 1-base whereas BED coordinates start at 0
        #so the start of the exon is in column 4 in the gtf file but we subtract 1 as BED coordinates are gtf coordinates -1 (only do this for start not end)
        #we also need to make the text an integer for the computer to understand it as a number 
        start = int(fields[3]) - 1

        #end coord stays the same (not -1 as start already subtracted 1), it must still be changed from text into a number using int
        end = int(fields[4])

        #find the transcript strand
        strand = fields[6]

        #column 9 contains the gene and transcript information
        #however it is a list of things we will need to filter through and extract things from
        attributes = fields[8]

        #so we will extract the transcript ID from the attributes column
        #look for the text inside of the transcript id's " ", get stuff inside the brackets (), match anything that is not a ", match these caharcters + and then stop matching when you reach another "
        match = re.search(r'transcript_id "([^"]+)"', attributes)

        #if there is no transcript ID continue which means skip
        if not match:
            continue
    
        #this will take the transcript ID found by the previous search (by re.) and save it in a variable called transcript_id 
        transcript_id = match.group(1)

        #add the exon coordinates to the corresponding transcript
        transcripts[transcript_id].append((start, end))

        #store the strand and chromosome for the transcript
        strands[transcript_id] = strand
        chromosomes[transcript_id] = chromosome


#now this is going to write all of the above transcript information to a BED file
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
            f"{chrom}\t"                      #chromosome
            f"{transcript_start}\t"           #the transcript start
            f"{transcript_end}\t"             #the transcript end
            f"{transcript_id}\t"              #the transcript ID
            f"0\t"                            #the BED score
            f"{strand}\t"                     #strand
            f"{transcript_start}\t"           #thickStart (see descrip below)
            f"{transcript_end}\t"             #thickEnd (the way BED12 names the start and end of a feature, howver we want teh full transcript so we don't really use thick start and end)
            f"0\t"                           
            f"{block_count}\t"               
            f"{','.join(block_sizes)}\t"    
            f"{','.join(block_starts)}\n"     
        )

print(f"BED written to {output_bed}")
print(f"no. transcripts written {len(transcripts)}")
