#!/usr/bin/env python3

#EXPLANATION OF THE SCRIPT 
#this script reads in the tab-delimited BLAST output (.tsv) that was in the format as stipulated in our BLASTN parameters:
# -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

#The output showed that a single human PCG often aligned to the several different mouse PCGs so we need to select which alignment was the bedt 

#import csv so that python can read and write tab-delimited files
import csv
#imort sys so filenames can be supplied into the terminal
import sys
#import math so that -log10(E-vlaue) can be calculated
import math
#defaultdict automatically creates nested dictionaries when a new query or subject gene is encountered
from collections import defaultdict
#path is used to check whether teh input file exists
from pathlib import Path


#these ouput column names are the headers of the output tsv file

OUTPUT_COLUMN_NAMES = [

    #Human query PCG info
    "qseqid",
    "query_chromosome",
    "query_gene_start",
    "query_gene_end",
    "query_strand",

    #info about our subject sequence (the species being comapred against the human query)
    "sseqid",
    "subject_chromosome",
    "subject_gene_start",
    "subject_gene_end",
    "subject_strand",

    #combined blast statistics
    "number_of_alignments",
    "total_alignment_length",
    "weighted_percent_identity",
    "total_identical_bases",
    "total_bitscore", 
    "weighted_log_evalue",
    "best_evalue",

    #final rank of the mouse candidate for this human PCG
    "rank_for_query"
]

#we use parse gene location to extract the chromosome, genomic start, genomic end and strand from a BLAST sequence identifier

def parse_gene_location(sequence_id):
    
    try: 
        #split teh gene name from the genomic coords (eg OR4F5::chr1:65418-71585(+) becomes:
        # OR4F5
        # chr1:65418-71585(+))
        location = sequence_id.split("::", 1)[1]

        #split chr from coords and strand
        chromosome, coordinates_and_strand= location.split(":",1)

        #split the genomic coords from strand
        coordinates, strand_text = coordinates_and_strand.rsplit("(", 1)
                                                             
        #sep the genomic start and end
        start, end = coordinates.split("-", 1)

        #remove the closing bracket from the strand
        strand = strand_text.rstrip(")")

        return chromosome, start, end, strand 
    

    except (IndexError, ValueError):
        #return blank fields if identifier does not match the expected format
        return "", "", "", ""


#to find top hits/matches for each human PCG
def find_top_hits(input_file, output_file):

    # we are creating a nested dictionary (this is a dictionary of dictionaries) to store our combined statistics
    #we use pair statistics  because we are storing statistics for each human-mouse gene pair
    pair_statistics = defaultdict(
        lambda: defaultdict(
            lambda: {
                #number of separate BLAST alignments for this human-mouse gene pair
                "number_of_alignments":0,

                #sum of all alignment lengths
                "total_alignment_length":0,

                #sum of percentage identity x alignment length
                #(this is later divided by the total alignment length to calculate the length-weighted pident)
                "weighted_identity_sum": 0.0, 

                #sum of all BLAST bitscores
                "total_bitscore": 0.0,

                #sum of:
                #-log10(e-value) x alignment length
                #this is later divided by total alignment length to give a combined length-weighted e-value score
                "weighted_log_evalue_sum": 0.0,

                #smallest individual e-value, kept to check incase selection is wrong
                "best_evalue": float("inf")
            }
        )
    )

    #count the total number of blast rows examined
    total_blast_lines = 0

    #count invalid rows that are skipped 
    invalid_blast_lines = 0

  #now we need to see the BLASTN output file and read it a line at a time
    with open(input_file, "r", newline="") as blast_file:

        #tell python that the file is tab-delimited 
        reader = csv.reader(blast_file, delimiter = "\t")

        #read one BLAST alignment at a time until the end of the file
        for line_number, row in enumerate(reader, start = 1):

            #ignore lines with ntohing in them (eg: blank lines at the end of the file)
            if not row:
                continue

            total_blast_lines +=1

            #the format outfmt 6 line (as specified in our BLAST .sh files) should contain 12 columns for this analysis
            if len(row) != 12:
                print(
                    f"skipped BLASTN line{line_number}:"
                    f"as we expected 12 columns but found {len(row)}"
                )

                invalid_blast_lines += 1
                continue

            #column 1 is the human query PCG
            query_id = row[0]

            #column 2 is the mouse subject PCG
            subject_id = row[1]

            #convert the numerical BLAST fields from text into python numbers
            try:
                #column 3 is the percentage identity (float is a way to store decimal numbers in python)
                pident = float(row[2])

                #column 4 is the alignment length
                alignment_length = int(row[3])

                #column 11 is the e-value
                evalue = float(row[10])

                #column 12 is the bitscore
                bitscore = float(row[11])

            except ValueError:
                print(
                    f"skipped BLASTN line{line_number}:"
                    f"a numeric field could not be converted"
                )

                invalid_blast_lines +=1
                continue

            #alignment length must be greater than 0 because it will later be used as a divisor to calculate length-weighted statistics
            if alignment_length <= 0:
                print(
                    f"skipped BLASTN line{line_number}:"
                    f"alignment length must be greater than zero"
                )

                invalid_blast_lines +=1
                continue

            #we need to access the stored statistics for this exact human and subject pair
            stats = pair_statistics[query_id][subject_id]
            
            #count this BLAST alignment fragment
            stats["number_of_alignments"] += 1

            #add this fragment's length to the total aligned length for the gene pair
            stats["total_alignment_length"] += alignment_length

            #TO COMBINE PERCENTAGE IDENTITY AND ALIGNMENT LENGTH (here is an example)
            #
            #1st alignment: 83% idnetity over 1000 bp
            #2nd alignmnet: 82% idnetity over 100bp
            #
            #weighted identity
            # [(83x1000)+(82x100)] / 1000 + 100 = 82.91%
            
            stats["weighted_identity_sum"] += (
                pident * alignment_length
            )

            stats["total_bitscore"] += bitscore

            #WE ARE COMBINING E VALUES AND ALIGNMENT LENGTHS
            #the raw e-values should not be averaged directly as it spans a large range of numbers so they must be converted to -log10(e-value)
            #a larger transformed value serves as evidence of a stronger alignmnet 

            #BLAST may report strong evalues as 0.0 so since log10(0) can't be calculated, use a very small positive evalue
            safe_evalue = max(evalue, 1e-300)

            #convert the raw e-value to -log10(evalue)
            log_evalue_score = -math.log10(safe_evalue)

            #weight the transformed e-value by alignment length
            stats["weighted_log_evalue_sum"] += (
                log_evalue_score * alignment_length
            )

            #keep the smallest individual e-value for reference
            stats["best_evalue"] = min(
                stats["best_evalue"],
                evalue        
            )

   #we want to rank the matches/the hits so we need to calculate some combined statistics 
    ranked_results = []

    #process one human query at a time
    for query_id, subject_hits in pair_statistics.items():

        #this list stores all comparison species candidates for the current human query PCG
        candidates = []

        #process one species candidate PCG at a time
        for subject_id, stats in subject_hits.items():

            #total aligned length across all fragments
            total_length = stats["total_alignment_length"]

            #calculate the length-weighted percentage identity
            weighted_percentage_identity = (
                stats["weighted_identity_sum"]/total_length
            )

            #working out the total number of identical bases across all alignments
            #an example would be: 83% identity over 1000bp => 83x1000/100 = 830 identical bases
            total_identical_bases = (
                stats["weighted_identity_sum"]/100.0
            )

            #calc the length-weighted log e-value
            weighted_log_evalue = (
                stats["weighted_log_evalue_sum"]/total_length
            )
 
            #human query gene coords
            (
                query_chromosome,
                query_gene_start,
                query_gene_end,
                query_strand
            ) = parse_gene_location(query_id)

            #extract the subject gene coords
            (
                subject_chromosome,
                subject_gene_start,
                subject_gene_end,
                subject_strand
            ) = parse_gene_location(subject_id)

            #store the combined result for this human vs comparison subject pair

            candidate = {

                #human PCG info
                "qseqid": query_id,
                "query_chromosome": query_chromosome,
                "query_gene_start": query_gene_start,
                "query_gene_end": query_gene_end,
                "query_strand": query_strand,

                #subject PCG info
                "sseqid": subject_id,
                "subject_chromosome": subject_chromosome,
                "subject_gene_start": subject_gene_start,
                "subject_gene_end": subject_gene_end,
                "subject_strand": subject_strand,

                #combined statistics
                "number_of_alignments": (
                    stats["number_of_alignments"]
                ),
                "total_alignment_length": total_length,
                "weighted_percent_identity": (
                    weighted_percentage_identity
                ),
                "total_identical_bases": (
                    total_identical_bases
                ),
                "total_bitscore": (
                    stats["total_bitscore"]
                ), 
                "weighted_log_evalue": (
                    weighted_log_evalue
                ),
                "best_evalue": (
                    stats["best_evalue"]
                )                
            }
            candidates.append(candidate)
                
                

        #This is the order of preference for ranking the candidates: 
        #python sorts by the first value and if equal then moves to the second for comparison and so on
        candidates.sort(
            key=lambda hit: (
                -hit["total_identical_bases"],
                -hit["total_bitscore"],
                -hit["weighted_percent_identity"],
                -hit["total_alignment_length"],
                -hit["weighted_log_evalue"]
            )
        )

        #keep only the best subject candidate for this human PCG
        for rank, candidate in enumerate(candidates, start=1):
            candidate["rank_for_query"] = rank
            ranked_results.append(candidate)

       #we wanted a separate file that has only the top ranked blast match

        top1_output = output_file.replace(".tsv", "_Top1.tsv")

        with open(top1_output, "w", newline="") as outfile:

            writer = csv.DictWriter(
                outfile,
                fieldnames=OUTPUT_COLUMN_NAMES,
                delimiter="\t"
            )

            writer.writeheader()

            for result in ranked_results:

                if result["rank_for_query"] == 1:

                    result_to_write = result.copy()

                    result_to_write["weighted_percent_identity"] = (
                        f'{result["weighted_percent_identity"]:.3f}'
                    )

                    result_to_write["total_identical_bases"] = (
                        f'{result["total_identical_bases"]:.2f}'
                    )

                    result_to_write["total_bitscore"] = (
                        f'{result["total_bitscore"]:.2f}'
                    )

                    result_to_write["weighted_log_evalue"] = (
                        f'{result["weighted_log_evalue"]:.3f}'
                    )

                    result_to_write["best_evalue"] = (
                        f'{result["best_evalue"]:.6g}'
                    )

                    writer.writerow(result_to_write)

    #order the output so that all results for the same human query are together, showing rank 1 first
    ranked_results.sort(
        key=lambda hit: (     
            hit["qseqid"],
            hit["rank_for_query"]
        )
    )
    
    #now we need to write the output file with the ranked results
    with open(output_file, "w", newline = "") as outfile:

        #create a tab-delimited dictionary writer
        writer = csv.DictWriter(
            outfile, 
            fieldnames= OUTPUT_COLUMN_NAMES,
            delimiter="\t"
        )

        #write the header row
        writer.writeheader()

        #write one combined human-subject gene pair per row
        for result in ranked_results:

            #copy the dictionary so that the formatting does not change the original values
            result_to_write = result.copy()

            #format weighted percentage identity to 3 decimals
            result_to_write["weighted_percent_identity"] =(
                f'{result["weighted_percent_identity"]:.3f}'
            )
            #format estimated identical bases to 2 decimals
            result_to_write["total_identical_bases"] = (
                f'{result["total_identical_bases"]:.2f}'
            )
            #format summed bitscore to 2 decimals
            result_to_write["total_bitscore"] = (
                f'{result["total_bitscore"]:.2f}'
            )
            #format transformed e-value score to 3 decimals
            result_to_write["weighted_log_evalue"] = (
                f'{result["weighted_log_evalue"]:.3f}'
            )
            #format the best raw evalue using scientific notation
            result_to_write["best_evalue"] = (
                f'{result["best_evalue"]:.6g}'
            )
            writer.writerow(result_to_write)


    number_of_queries = len(pair_statistics)
    number_of_gene_pairs = len(ranked_results)

    print("\nselection complete")
    print(f"lines skipped: {invalid_blast_lines}")
    print(f"human query PCGs identified {number_of_queries}")
    print(f"human-subject pairs ranked {number_of_gene_pairs}")
    print("to choose the top hit take the row where rank_for_query equals 1")


#need to write this if we want to run the script from teh command line to tell it that this si the main function 
#run it by python3 BLASTTop_hits.py input.tsv output.tsv

def main():

    if len(sys.argv) != 3:   #to check I have given this file name, the input and output file names 

        print("can't use all 3 files")

        sys.exit(1)

    #get the input and output file names from the command line arguments
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    #if the files don't exist then say error
    if not Path(input_file).is_file():

        print(
            f"Error: input file '{input_file}' does not exist"
        )

        sys.exit(1)

    # Run the BLAST grouping and ranking analysis.
    find_top_hits(input_file, output_file)

#need to write this if we want to run the script from teh command line to tell it that this si the main function 
#run it by python3 BLASTTop_hits.py input.tsv output.tsv
if __name__ == "__main__":
    main() 
