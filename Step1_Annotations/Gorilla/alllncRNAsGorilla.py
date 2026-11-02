#this script is generated to be universal to extract lncRNA gene blocks from a GTF or gencode file. It will work for any species as long as the GTF file is formatted correctly and contains the necessary attributes to identify lncRNAs.
#this script is used on gorilla 

"""
Extract lncRNA gene blocks from a GTF file.

The script:

1. Reads transcript entries from a GTF file.
2. Identifies lncRNA transcripts using either:
       gene_type
       gene_biotype
       transcript_type
       transcript_biotype
3. Groups all transcript isoforms belonging to the same gene.
4. Creates one lncRNA block per gene:
       block start = smallest isoform start
       block end   = largest isoform end
5. Writes the lncRNA name, chromosome, start, end, strand
   and number of isoforms to a TSV file.
"""

import argparse        #this module is used to parse command line arguments (parse means to analyze and extract information from a string)
import csv             #this module is used to read and write tabular data in CSV format
import sys             #this module provides access to some variables used or maintained by the interpreter and to functions that interact strongly with the interpreter
from collections import defaultdict        #this module provides a dictionary-like object that allows you to create a dictionary with default values for non-existent keys
from pathlib import Path        #this module provides an object-oriented interface for working with file system paths

##### A lncRNA can be any of these biotypes/gene types
#####EDIT IF NEEDED

LNCRNA_BIOTYPES = {
    "lncRNA",
    "lincRNA",
    "antisense",
    "misc_RNA",
    "sense_intronic",
    "sense_overlapping",
    "3prime_overlapping_ncrna",
    "bidirectional_promoter_lncRNA",
    "non-coding",
    "miRNA",
}

#We need to convert the 9th column in the GTF file into a dictionary so we can easily access the attributes of each transcript. The 9th column is a semicolon-separated list of key-value pairs, where the key is the attribute name and the value is the attribute value. We will use a function to parse this column and return a dictionary of attributes.
'''Example input:
        gene_id "ENSMUSG00000102693";
        transcript_id "ENSMUST00000193812";
        gene_type "lncRNA";
        gene_name "4933401J01Rik";

    Example output:
        {
            "gene_id": "ENSMUSG00000102693",
            "transcript_id": "ENSMUST00000193812",
            "gene_type": "lncRNA",
            "gene_name": "4933401J01Rik"
        }
    '''
def parse_attributes(attribute_text):

    attributes = {}

    for item in attribute_text.strip().split(";"):
        item = item.strip()                            #this removes any leading or trailing whitespace from the item

        if not item:
            continue

        #the GTF typical format is key "value", so we need to split on the first space and remove the quotes from the value
        #each individual attribute is a key-value pair separated by a space, and the value is enclosed in double quotes. We will split the item into two parts: the key and the value, using the first whitespace as the delimiter. Then we will remove the quotes from the value and store it in the attributes dictionary.
        parts = item.split(None, 1)      #this splits the item into two parts: the key and the value, using the first whitespace as the delimiter

        if len(parts) != 2:
            continue

        key = parts[0]
        value = parts[1].strip('"')     #this removes the quotes from the value

        attributes[key] = value

    return attributes


#to ensure that each chromosome number has a chr at the start of the string, we will use a function to add the prefix "chr" to the chromosome number if it does not already have it. This is important because some GTF files may use different naming conventions for chromosomes, and we want to standardize them for downstream analysis.
def add_chr_prefix(chromosome):

    chromosome = chromosome.strip()  #this removes any leading or trailing whitespace from the chromosome string

    if chromosome.lower().startswith("chr"):
        return chromosome

    return "chr" + chromosome

def get_biotype(attributes):

    #get teh transcript or biotype from the gtf or gencode

    return(
        attributes.get("gene_type")
        or attributes.get("gene_biotype")
        or attributes.get("transcript_type")
        or attributes.get("transcript_biotype")
        or ""
    )

###this is only removing chr for sorting
def chromosome_sort_key(chromosome):

    chromosome_without_prefix = chromosome        #this removes the "chr" prefix from the chromosome string if it exists, so that we can sort the chromosomes numerically or alphabetically without the prefix interfering with the sorting order.

    if chromosome.lower().startswith("chr"):
        chromosome_without_prefix = chromosome[3:]

    if chromosome_without_prefix.isdigit():
        return 0, int(chromosome_without_prefix)

    chromosome_order = {
        "X": 23,
        "Y": 24,
        "M": 25,
        "MT": 25,
    }

    upper_name = chromosome_without_prefix.upper()        #chr:  x -> X

    if upper_name in chromosome_order:
        return 0, chromosome_order[upper_name]

    return 1, chromosome_without_prefix

def read_lncrna_isoforms(input_file):

    genes = {}

    total_lines = 0
    transcript_rows = 0
    lncrna_transcript_rows = 0
    missing_gene_id_rows = 0
    conflicting_records = 0

    with open(input_file, "r") as gtf:

        for line_number, line in enumerate(gtf, start=1):         #this iterates over each line in the GTF file, starting the line number count at 1

            total_lines +=1

            if not line.strip():       #this checks if the line is empty or contains only whitespace, and if so, it skips to the next iteration of the loop
                continue

            if line.startswith("#"):        #this checks if the line is a comment line (starting with "#"), and if so, it skips to the next iteration of the loop
                continue

            fields = line.rstrip("\n").split("\t")     #this splits the line into fields using the tab character as the delimiter, and removes the trailing newline character from the line 

            if len(fields) != 9:      #this checks if the line has exactly 9 fields, and if not, it skips to the next iteration of the loop
                continue 

            chromosome = add_chr_prefix(fields[0])
            feature_type = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attribute_text = fields[8]

            # Only transcript rows represent complete isoforms.
            if feature_type != "transcript":
                continue

            transcript_rows += 1

            if strand not in {"+", "-"}:
                continue

            try:
                start = int(fields[3])
                end = int(fields[4])

            except ValueError:

                print(
                    f"Warning: invalid coordinates on line {line_number}",
                    file=sys.stderr,
                )

                continue

            attributes = parse_attributes(attribute_text)

            biotype = get_biotype(attributes)

            if biotype not in LNCRNA_BIOTYPES:
                continue

            lncrna_transcript_rows += 1

            gene_id = attributes.get("gene_id")

            if not gene_id:
                missing_gene_id_rows +=1
                continue

            gene_name = attributes.get("gene_name", gene_id)     #this gets the gene name from the attributes dictionary, and if it does not exist, it uses the gene ID as the gene name
            transcript_id = attributes.get("transcript_id")          #the transcript id is a unique identifier for each transcript isoform, and it is used to group all isoforms belonging to the same gene

            if gene_id not in genes:
 
                genes[gene_id] = {            #this creates a new key in the genes dictionary for the gene ID, and initializes it with a dictionary containing the gene name, chromosome, strand, and an empty list of isoforms
                    "gene_name": gene_name,
                    "chromosome": chromosome,
                    "strand": strand,
                    "starts": [],
                    "ends": [],
                    "transcript_ids": set(),         #set() is a built-in Python data structure that stores unique elements, and it is used here to store the transcript IDs of all isoforms belonging to the same gene, ensuring that there are no duplicates
                    "biotypes": set(),
                }

            gene_record = genes[gene_id]

            #all isoforms of one gene should be on the same chromosome and strand
            if ( 
                gene_record["chromosome"] != chromosome
                or gene_record["strand"] != strand
            ):
                conflicting_records += 1

                print(
                    f"Warning: gene{gene_id} occurs on >1 chr or strand so skipping line {line_number}",
                    file=sys.stderr,
                )

                continue

            gene_record["starts"].append(start)
            gene_record["ends"].append(end)
            gene_record["transcript_ids"].add(transcript_id)
            gene_record["biotypes"].add(biotype)

            if transcript_id:
                gene_record["transcript_ids"].add(transcript_id)

    statistics = {
        "total_lines": total_lines,
        "transcript_rows": transcript_rows,
        "lncrna_transcript_rows": lncrna_transcript_rows,
        "missing_gene_id_rows": missing_gene_id_rows,
        "conflicting_records": conflicting_records,
    }

    return genes, statistics

def make_lncrna_blocks(genes):

    lncrna_blocks = []  #using [] to create an empty list to store the lncRNA blocks, which will be a list of dictionaries containing the gene name, chromosome, start, end, strand, number of isoforms, and biotypes for each lncRNA gene.
    for gene_id, gene_record in genes.items():

        if not gene_record["starts"]:           #this checks if the "starts" list in the gene record is empty, and if so, it skips to the next iteration of the loop
            continue

        if not gene_record["ends"]:             #this checks if the "ends" list in the gene record is empty, and if so, it skips to the next iteration of the loop
            continue

        block_start = min(gene_record["starts"])       #this calculates the start of the lncRNA block as the minimum of all isoform starts for the gene
        block_end = max(gene_record["ends"])           #this calculates the end of the

        if gene_record["transcript_ids"]:
            isoform_count = len(gene_record["transcript_ids"])        #this calculates the number of isoforms for the gene as the length of the "transcript_ids" set in the gene record
        else:
            isoform_count = len(gene_record["starts"])        #if the "transcript_ids" set is empty, it calculates the number of isoforms as the length of the "starts" list in the gene record

        lncrna_blocks.append(
            {
                "lncRNA_name": gene_record["gene_name"],
                "lncrna_gene_id": gene_id,
                "chromosome": gene_record["chromosome"],
                "start": block_start,
                "end": block_end,
                "strand": gene_record["strand"],
                "isoform_count": isoform_count,
                "biotypes": ",".join(sorted(gene_record["biotypes"])),        #this joins the biotypes in the "biotypes" set in the gene record into a comma-separated string, sorted alphabetically
            }
        )    

    lncrna_blocks.sort(
        key=lambda row: (
            chromosome_sort_key(row["chromosome"]),
            row["start"],
            row["end"],
            row["lncrna_gene_id"],
        )
    )

    return lncrna_blocks

def write_lncrna_blocks(lncrna_blocks, output_file):       #write teh lncRNA blocks to a TSV file, with the columns: lncRNA name, gene ID, chromosome, start, end, strand, number of isoforms, and biotypes.

    fieldnames = [
        "lncRNA_name",
        "lncrna_gene_id",
        "chromosome",
        "start",
        "end",
        "strand",
        "isoform_count",
        "biotypes",
    ]

    with open(
        output_file, 
        "w", 
        newline=""
    ) as output:


        writer = csv.DictWriter(
            output,
            fieldnames = fieldnames,
            delimiter="\t"
        )

        writer.writeheader()
        writer.writerows(lncrna_blocks)


def parse_arguments():      #this will read the input and output file names from the command line

    parser = argparse.ArgumentParser(
        description=(
            "Extract lncRNA isoforms from GTF file and condense them all into a single lncRNA gene block per lncRNA"
        )
    )

    parser.add_argument(
        "input_file",
        help = "Input GTF annotation file"
    )

    parser.add_argument(
        "-o",
        "--output",
        default = "lncRNA_gene_blocks.tsv",
        help = "Output TSV file"
    )

    return parser.parse_args()

def main():
    args = parse_arguments()

    gtf_path = Path(args.input_file)

    if not gtf_path.is_file():
        sys.exit(
            f"Error:GTF file not found: {gtf_path}"
        )

    print(f"Reading GTF file: {gtf_path}")

    genes, statistics = read_lncrna_isoforms(gtf_path)

    lncrna_blocks = make_lncrna_blocks(genes)

    write_lncrna_blocks(
        lncrna_blocks,
        args.output,
    )

    total_isoforms = sum(
        block["isoform_count"]
        for block in lncrna_blocks
    )

    print()
    print("finihsed extracting lncRNAs")
    print(
        f"total GTF lines read:"
        f"{statistics['total_lines']}"
    )

    print(
        f"Transcript rows examined: "
        f"{statistics['transcript_rows']}"
    )
    print(
        f"lncRNA transcript rows found: "
        f"{statistics['lncrna_transcript_rows']}"
    )
    print(
        f"lncRNA gene blocks written: "
        f"{len(lncrna_blocks)}"
    )
    print(
        f"Total lncRNA isoforms represented: "
        f"{total_isoforms}"
    )
    print(
        f"Rows missing gene_id: "
        f"{statistics['missing_gene_id_rows']}"
    )
    print(
        f"Conflicting records skipped: "
        f"{statistics['conflicting_records']}"
    )
    print(
        f"Results saved to: {args.output}"
    )


if __name__ == "__main__":
    main()