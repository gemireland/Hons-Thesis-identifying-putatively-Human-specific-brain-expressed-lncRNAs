#this script will be used as a sanity check to check that my list of lncRNAs that are conserved between 
#human and the other 3 species (mouse, chimp and gorilla) are correct
#I am mainly using this script to ensure that if the human lncRNA overlaps with its PCG that there is also an overlap with the mouse lncRNA and its PCG


import csv

#when I use this script on other species I must change this part to suit teh species
species = "chimp"
#species = "mouse"
#species = "gorilla"

conserved_file = f"Conserved_{species.capitalize()}_lncRNAs.tsv"
synteny_file = f"Human_{species.capitalize()}_lncRNA_Synteny_Results.tsv"
output_file = f"Conserved_{species.capitalize()}_lncRNAs_Cleaned.tsv"

def overlap_relationship(lnc_start, lnc_end, pcg_start, pcg_end):

    #If there is no overlap 
    if lnc_end < pcg_start or lnc_start > pcg_end:
        return "NO_OVERLAP"

    #a situation if the lncRNA is completely overlapped with the PCG
    if lnc_start >= pcg_start and lnc_end <= pcg_end:
        return "COMPLETE_OVERLAP"

    #situation of partial overlap (this is the last case: otherwise)
    return "PARTIAL_OVERLAP"


#read in the full synteny results file
synteny = {}

with open(synteny_file, "r", newline="") as infile:
    reader = csv.DictReader(infile, delimiter="\t")

    for row in reader:
        synteny[row["lncRNA_gene_id"]] = row

kept = 0
removed = 0
not_assessable = 0

with open(conserved_file, "r", newline = "") as infile, \
    open(output_file, "w", newline="") as outfile:

    reader = csv.DictReader(infile, delimiter="\t")

    fieldnames = reader.fieldnames + ["PCG_overlap_conserved"]

    writer = csv.DictWriter(
        outfile,
        fieldnames=fieldnames,
        delimiter="\t"
    )

    writer.writeheader()

    for row in reader:

        result = synteny.get(row["human_gene_id"])

        if result is None:
            continue

        species_pcg_chr = result[
            f"closest_PCG_{species}_chromosome"
        ]

        species_lnc_chr = result[
            f"selected_{species}_lncRNA_chromosome"
        ]

        #can't compare coordinates on different chromosomes
        if (
            species_pcg_chr == ""
            or species_lnc_chr == ""
            or species_pcg_chr != species_lnc_chr
        ):
            row["PCG_overlap_conserved"] = "NOT_ASSESSABLE"
            writer.writerow(row)

            not_assessable += 1
            kept +=1
            continue

        #human lncRNA-PCG overlap
        human_relationship = overlap_relationship (
            int(result["lncRNA_start"]),
            int(result["lncRNA_end"]),
            int(result["closest_PCG_start"]),
            int(result["closest_PCG_end"])
        )

        #Species lncRNA-PCG overlap
        species_relationship = overlap_relationship(
            int(result[f"selected_{species}_lncRNA_start"]),
            int(result[f"selected_{species}_lncRNA_end"]),
            int(result[f"closest_PCG_{species}_start"]),
            int(result[f"closest_PCG_{species}_end"])
        )

        #I only want to keep the conserved lncRNA pair when the relationship of the lncRNA to the gene of overlap or the non-overlap is conserved
        #if there is partial overlap in human then you can accept partial overlap or adjacent ones but you have to exclude complete overlap 
        #if there is complete overlap in humans you can only accept complete overlap in other species not partial or adjacent
        keep_pair = False

        #If human has no overlap species must also have no overlap
        if human_relationship == "NO_OVERLAP":
            if species_relationship == "NO_OVERLAP":
                keep_pair = True

        #human has partial overlap so we can accept only a partial overlap or no overlap (allowing adjacent to be kept)
        elif human_relationship == "PARTIAL_OVERLAP":
            if species_relationship in ("PARTIAL_OVERLAP", "NO_OVERLAP"):
                keep_pair = True

        #if teh human lncRNA completely overlaps the PCG then the species must show cpmplete oevrlap
        elif human_relationship == "COMPLETE_OVERLAP":
            if species_relationship == "COMPLETE_OVERLAP":
                keep_pair = True


        if keep_pair:

            row["PCG_overlap_conserved"] = "TRUE"
            writer.writerow(row)

            kept += 1

        else:
            removed +=1

print("Saved:", output_file)
print("Kept:", kept)
print("Removed", removed)
print("Not assesssable:", not_assessable)

