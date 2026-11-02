#!/usr/bin/env python3

#---------------- EXPLANATION OF THE SCRIPT ----------------
#this script reads in the tab-delimited BLAST output (.tsv) that was in the format :
# -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

#The output showed that a single human PCG often aligned to the several different mouse PCGs. 
#Also the same human-mouse PCG pair can have several separate BLAST alignment fragments

#so I have created this script to:
#1. Read every line of the BLATS output in the .tsv file
#2. Groups alignmnets by exact human query PCG and mouse subject PCG
# 3. Combines all BLAST alignment fragments that belong to the same human query PCG and the same mouse subject PCG
# example: Human OR4F5 -> Mouse OR4F5. may have several separate BLAST alignment fragments, but they will be combined into a single alignment
# 
# 4. Calculates combined similarity statistics 
# 5. Extracts the genomic coordinates and strand for both genes
# 6. Ranks all mouse candidates separately for each human query PCG based on the combined similarity statistics
#  (the top mouse candidate for each human PCG will have: rank_for_query =1)

#-------- summary done --------

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

# ------final output column layout---------
#these ouput column names are the headers of the output tsv file

OUTPUT_COLUMN_NAMES = [

    #Human query PCG info
    "qseqid",
    "query_chromosome",
    "query_gene_start",
    "query_gene_end",
    "query_strand",

    #Mouse subject info
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

#---using parse_gene_location -------

def parse_gene_location(sequence_id):
    
    """ Extract chrom, genomic start, genomic end and strand from a BLAST sequence identifier
    
    format:
        GeneName::chr1:1234-12345(+)
    
    """

    
    try: 
        #split teh gene name from the genomic location (eg OR4F5::chr1:65418-71585(+) becomes:
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


#--------- to find_top_hits for human PCG ---------

def find_top_hits(input_file, output_file):

    """ Read the full BLAST TVS file, combine all the alignments that belong to the same human 
    query PCG and mouse subject PCG, and rank the mouse candidates for every human query PCG """     

    #------------------------------------------
    # we are creating a NESTED DICTIONARY to store combined statistics
    #------------------------------------------
    #format: 
    #          pair_statistics[query_id][subject_id]
    # example:
    #         pair_statistics["OR4F5"]["Orf4f5"]

    #-----------------------------------

    pair_statistics = defaultdict(
        lambda: defaultdict(
            lambda: {
                #Number of separate BLAST fragments 
                "number_of_alignments":0,

                #sum of all alignment lengths
                "total_alignment_length":0,

                #sum of: pident x alignment length
                #(this is later divided by the total alignment length to calculate the length-weighted pident)
                "weighted_identity_sum": 0.0, 

                #sum of all BLAST bitscores
                "total_bitscore": 0.0,

                #sum of:
                #-log10(e-value) x alignment length
                #this is later divided by total alignment length to give a combined length-weighted e-value score
                "weighted_log_evalue_sum": 0.0,

                #smallest individual e-value, kept for reference
                "best_evalue": float("inf")
            }
        )
    )

    #count the total number of blast rows examined
    total_blast_lines = 0

    #count invalid rows that are skipped 
    invalid_blast_lines = 0

    #-------------------------
    #open and read the complete BLAST TSV file
    #--------------------------------------

    with open(input_file, "r", newline="") as blast_file:

        #tell python that the file is tab-delimited 
        reader = csv.reader(blast_file, delimiter = "\t")

        #read one BLAST alignment at a time until the end of the file
        for line_number, row in enumerate(reader, start = 1):

            #ignore empty lines
            if not row:
                continue

            total_blast_lines +=1

            #the format outfmt 6 line should contain 12 columns for this analysis
            if len(row) != 12:
                print(
                    f"Skipping invalid BLAST line{line_number}:"
                    f"expected 12 columns, found only {len(row)}"
                )

                invalid_blast_lines += 1
                continue

            #column 1 is the human query PCG
            query_id = row[0]

            #column 2 is the mouse subject PCG
            subject_id = row[1]

            #convert the numerical BLAST fields from text into python numbers
            try:
                #column 3:percentage identity 
                pident = float(row[2])

                #column 4: alignment length
                alignment_length = int(row[3])

                #column 11 is the e-value
                evalue = float(row[10])

                #column 12: bitscore
                bitscore = float(row[11])

            except ValueError:
                print(
                    f"Skipping invalid BLAST line{line_number}:"
                    "a numeric field could not be converted"
                )

                invalid_blast_lines +=1
                continue

            #alignment length must be greater than 0 because it will later be used as a divisor
            if alignment_length <= 0:
                print(
                    f"Skipping invalid BLAST line{line_number}:"
                    "alignment length must be greater than zero"
                )

                invalid_blast_lines +=1
                continue

            #we need to access the stored statistics for this exact human and subject pair
            stats = pair_statistics[query_id][subject_id]
            
            #count this BLAST alignment fragment
            stats["number_of_alignments"] += 1

            #add this fragment's length to the total aligned length for the gene pair
            stats["total_alignment_length"] += alignment_length

            #---------------------
            #COMBINING % IDENTITY
            #---------------------
            #
            #Alignment 1: 83% over 1000 bp
            #Alignment 2: 82% over 100bp
            #
            #weighted identity
            # ( 83 x 1000) + (82 x 100)
            # -------------------------
            #.    1000 + 100
            # --------------------------

            stats["weighted_identity_sum"] += (
                pident * alignment_length
            )

            stats["total_bitscore"] += bitscore

            #-------------------
            # combine E-values
            #-------------------
            #
            #raw v-values should not be averaged directly because they span a very large numerical range
            #instead e-values are convreted to -log10(e-value)
            #
            # 1e-10 is 10
            #1e -100 is 100
            #
            # a larger transformed value represents stronger alignmnet evidence

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

    #-------------------------------------------
    #CALCULATE COMBINED STATS TO RANK HITS
    #-------------------------------------------
    ranked_results = []

    #process one human query at a time
    for query_id, subject_hits in pair_statistics.items():

        #this list stores all mouse candidates for teh current human query PCG
        candidates = []

        #process one mouse candidate PCG at a time
        for subject_id, stats in subject_hits.items():

            #total aligned length across all fragments
            total_length = stats["total_alignment_length"]

            #calculate the length-eighted percentage identity
            weighted_percentage_identity = (
                stats["weighted_identity_sum"]/total_length
            )

            #estimate the total number of identical aligned bases
            #
            #example: 83% identity over 1000bp: 83x1000/100 = 830 identical bases
            total_identical_bases = (
                stats["weighted_identity_sum"]/100.0
            )

            #calc the length-weighted transformed evalue
            weighted_log_evalue = (
                stats["weighted_log_evalue_sum"]/total_length
            )
 
            #extract human query gene coords
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

            #store the combined result for this human-subject pair

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

                #combined BLAST statistics
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

        # ----------------------------------------------------
        # RANK THE MOUSE CANDIDATES FOR THIS HUMAN PCG
        # ----------------------------------------------------
        #
        # Ranking priority:
        #
        # 1. Highest weighted percentage identity
        # 2. Longest total alignment length
        # 3. Highest total bitscore
        # 4. Highest weighted -log10(E-value)
        #
        # Python sorts from smallest to largest by default.
        # Negative signs make larger values appear first.
        # ----------------------------------------------------

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

        # Keep only the best mouse candidate for this human PCG
        for rank, candidate in enumerate(candidates, start=1):
            candidate["rank_for_query"] = rank
            ranked_results.append(candidate)

        # ----------------------------------------
        # Write a second file containing only the
        # top-ranked mouse hit for each human PCG
        # ----------------------------------------

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

    # Sort the final output so that all results for the same
    # human query stay together and rank 1 is shown first.
    ranked_results.sort(
        key=lambda hit: (     
            hit["qseqid"],
            hit["rank_for_query"]
        )
    )

    #--------------------------------
    # WRITE THE OUTPUT TSV FILE
    #--------------------------------

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

    
    # --------------------------------------------------------
    # PRINT A SUMMARY TO THE TERMINAL
    # --------------------------------------------------------

    number_of_queries = len(pair_statistics)
    number_of_gene_pairs = len(ranked_results)

    print("\nFinished!")
    print(f"Input file: {input_file}")
    print(f"Output file: {output_file}")
    print(f"BLAST lines examined: {total_blast_lines}")
    print(f"Invalid lines skipped: {invalid_blast_lines}")
    print(f"Human query PCGs found: {number_of_queries}")
    print(f"Human-mouse gene pairs ranked: {number_of_gene_pairs}")
    print(
        "The top subject hit for each human PCG is the row "
        "where rank_for_query equals 1."
    )


# ============================================================
# MAIN FUNCTION
# ============================================================

def main():

    # The script should be run with:
    #
    # python3 script_name.py input.tsv output.tsv
    #
    # Therefore:
    #
    # sys.argv[0] = script filename
    # sys.argv[1] = input BLAST TSV file
    # sys.argv[2] = output ranked TSV file
    #
    # There must be exactly 3 command-line arguments in total.
    if len(sys.argv) != 3:

        print(
            "Usage: python3 Mtopblasthits_ranked.py "
            "input.tsv output.tsv"
        )

        sys.exit(1)

    # Read the input and output filenames.
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # Stop with a clear error if the input file does not exist.
    if not Path(input_file).is_file():

        print(
            f"Error: input file '{input_file}' does not exist"
        )

        sys.exit(1)

    # Run the BLAST grouping and ranking analysis.
    find_top_hits(input_file, output_file)


# ============================================================
# RUN THE SCRIPT
# ============================================================

if __name__ == "__main__":
    main()       
