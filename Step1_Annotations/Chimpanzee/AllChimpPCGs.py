#this script is used as a secondary script to refine the first step, it builds on from the initial script finalchimpensstep1.py
#this script is written to extract all PCGs (not only the anchors that flank the lncRNA blocks) from Ensembl release 116 (PanTro3)
#here is the Ensembl file (Chimp's reference annotation)
input_file = 'Pan_troglodytes.Pan_tro_3.0.116.gtf'

#we want to create an output file that contains the PCGs in BED format so we can later use BEDTOOLS getfasta to extract sequences for use in BLASTN
output_file = "All_Chimp_PCGs.bed"

#used to extarct information about the gene from the attributes column within the GTF file
def parse_attributes(attribute_text):
    
    attributes = {}
    #there are different attributes within the attributes column so we want to separate the individual attributes 
    for attribute in attribute_text.split(";"):
        attribute = attribute.strip()

        #skip the lines where the attribute column is left empty as we won't know if it is coding or non-coding
        if not attribute:
            continue

        #this is done to separte the attribute name from its correpsonding value
        parts = attribute.split(" ", 1)

        if len(parts) ==2:
            key = parts[0]
            value = parts[1].strip().strip('"')
            attributes[key] = value

    return attributes

#=0 to initialise the count for the number of PCGs that will be identified, so it will add a value for each one identified
pcg_count = 0

with open(input_file, "r") as gtf, open(output_file, "w") as bed:

    #we want to process each line of the gtf file but ignore the header and the comment lines
    for line in gtf:

        if line.startswith("#"):
            continue

        #we then want to split each line into its 9 columns 
        fields = line.rstrip("\n").split("\t")

        #a useable GTF line must have 9 columns or it can not be processed 
        if len(fields) != 9:
            continue

        chrom = fields[0]
        feature_type = fields[2]
        start = int(fields[3])
        end = int(fields[4])
        strand = fields[6]
        attributes_text = fields[8]

        #use gene-level entries only to identify PCGs not transcripts or exons, so exclude the transcripts and exons 
        #keeping gene-level entires allows us to see the PCGs as a single genomic block that spans the entirety of the gene boundary rather 
        #we figured out this is an efficient method that encompasses what we did in our initial scripts (it is able to group the transcript isoforms already at gene-level)
        if feature_type != "gene":
            continue    

        attributes = parse_attributes(attributes_text)

        #ensembl uses gene_biotype, but this also supports files containing gene_type
        #this allows us to idnetify the gene biotype of the attribute column even if it is listed as gene type
        gene_biotype = attributes.get(
            "gene_biotype",
            attributes.get("gene_type", "")
        )

        #we are going to extract the genes that are protein-coding for this script
        if gene_biotype != "protein_coding":
            continue

        gene_id = attributes.get("gene_id", "unknown_gene_id")
        gene_name = attributes.get("gene_name", gene_id)

        #we know that GTF coords are 1-based starts but BED coords are 0-based so to change the GTF coords into BED format 
        #subtract 1 from start and leave the end unchanged
        bed_start = start - 1
        bed_end = int(end)

        #writing the bed files out to make our gene blocks 
        #essentially what we are doing is writing out each PCG as a single gene level block in BED format using its start and end coords
        bed.write(
            f"{chrom}\t{bed_start}\t{bed_end}\t"
            f"{gene_name}\t0\t{strand}\n"
        )

        pcg_count += 1

print(f"{output_file} created containing {pcg_count} PCGs")
print("total PCGs:", pcg_count)
