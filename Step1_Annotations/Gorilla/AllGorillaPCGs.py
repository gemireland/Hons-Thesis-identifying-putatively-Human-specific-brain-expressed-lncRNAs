#Input gorilla ensembl GTF file
input_file = 'Gorilla_gorilla.gorGor4.116.gtf'

#Output BED file containing every PCG
output_file = "All_Gorilla_PCGs.bed"

def parse_attributes(attribute_text):
    """ Convert the final GFT attributes column into a dictionary
    
    Example:
    gene_id "ENSPTRG00000012345";
    gene_biotype "protein_coding";
    """

    attributes = {}
    for attribute in attribute_text.split(";"):
        attribute = attribute.strip()
        
        if not attribute:
            continue

        parts = attribute.split(" ", 1)

        if len(parts) ==2:
            key = parts[0]
            value = parts[1].strip().strip('"')
            attributes[key] = value

    return attributes

pcg_count = 0

with open(input_file, "r") as gtf, open(output_file, "w") as bed:

    for line in gtf:

        #Ignore GTF header and comment lines
        if line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")

        #A valid GTF line must have 9 columns
        if len(fields) != 9:
            continue

        chrom = fields[0]
        feature_type = fields[2]
        start = int(fields[3])
        end = int(fields[4])
        strand = fields[6]
        attributes_text = fields[8]

        #use gene-level entries only to identify PCGs not transcripts or exons
        if feature_type != "gene":
            continue    

        attributes = parse_attributes(attributes_text)

        #ensembl commonly uses gene_biotype, but this also supports files containing gene_type
        gene_biotype = attributes.get(
            "gene_biotype",
            attributes.get("gene_type", "")
        )

        if gene_biotype != "protein_coding":
            continue

        gene_id = attributes.get("gene_id", "unknown_gene_id")
        gene_name = attributes.get("gene_name", gene_id)

        #GTF starts are 1-based; BD starts are 0-based, so subtract 1 from the GTF start coordinate
        bed_start = start - 1
        bed_end = int(end)

        #BED6 format: chrom, start, end, name, score, strand
        bed.write(
            f"{chrom}\t{bed_start}\t{bed_end}\t"
            f"{gene_name}\t0\t{strand}\n"
        )

        pcg_count += 1

print(f"Finished making {output_file} with {pcg_count} PCGs")
print("Total PCGs written:", pcg_count)