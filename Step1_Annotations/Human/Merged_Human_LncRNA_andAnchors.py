#!/usr/bin/env python3

import ast
import csv


# ------------------------------------------------------------
# INPUT AND OUTPUT FILES
# ------------------------------------------------------------

step1_file = "outputhumangencodestep1.txt"
anchor_file = "Human_lncRNA_PCG_anchors.txt"

output_file = "Human_lncRNA_Full_Anchor_Table.tsv"


# ------------------------------------------------------------
# READ THE HUMAN lncRNA COORDINATES FROM THE STEP 1 OUTPUT
# ------------------------------------------------------------

lncrna_coordinates = {}

inside_lncrna_blocks = False

with open(step1_file, "r") as infile:

    for raw_line in infile:

        line = raw_line.strip()

        if not line:
            continue

        # Start reading after this heading
        if line.startswith("lncRNA dictionary gene blocks"):
            inside_lncrna_blocks = True
            continue

        # Skip summary lines
        if line.startswith("Total lncRNA genes"):
            continue

        if line.startswith("Total lncRNA isoforms"):
            continue

        # Stop when the protein-coding gene section begins
        if line.startswith("PCG dictionary gene blocks"):
            inside_lncrna_blocks = False
            break

        if not inside_lncrna_blocks:
            continue

        fields = line.split()

        # Coordinate lines should look like:
        # LINC01409 chr1 778739 810066 + 46
        if len(fields) != 6:
            continue

        lncrna_name = fields[0]
        chromosome = fields[1]

        if not chromosome.startswith("chr"):
            continue

        try:
            start = int(fields[2])
            end = int(fields[3])
            strand = fields[4]
            isoform_count = int(fields[5])
        except ValueError:
            continue

        record = {
            "lncrna_name": lncrna_name,
            "lncrna_chr": chromosome,
            "lncrna_start": start,
            "lncrna_end": end,
            "lncrna_strand": strand,
            "lncrna_isoform_count": isoform_count
        }

        # Some lncRNA names may occur at more than one locus.
        # Store every record under that name.
        if lncrna_name not in lncrna_coordinates:
            lncrna_coordinates[lncrna_name] = []

        lncrna_coordinates[lncrna_name].append(record)


# ------------------------------------------------------------
# FUNCTION TO READ AN UPSTREAM OR DOWNSTREAM PCG TUPLE
# ------------------------------------------------------------

def parse_anchor(anchor_text):
    """
    Converts:

    ('OR4F5', 54758, 'chr1', 65419, 71585, '+')

    into a dictionary.

    If the anchor is None, all values are returned as blank.
    """

    if anchor_text == "None":
        return {
            "pcg_name": "",
            "distance": "",
            "chromosome": "",
            "start": "",
            "end": "",
            "strand": ""
        }

    try:
        anchor = ast.literal_eval(anchor_text)
    except (ValueError, SyntaxError):
        return {
            "pcg_name": "",
            "distance": "",
            "chromosome": "",
            "start": "",
            "end": "",
            "strand": ""
        }

    if not isinstance(anchor, tuple) or len(anchor) != 6:
        return {
            "pcg_name": "",
            "distance": "",
            "chromosome": "",
            "start": "",
            "end": "",
            "strand": ""
        }

    return {
        "pcg_name": anchor[0],
        "distance": anchor[1],
        "chromosome": anchor[2],
        "start": anchor[3],
        "end": anchor[4],
        "strand": anchor[5]
    }


# ------------------------------------------------------------
# READ THE HUMAN lncRNA–PCG ANCHOR FILE
# ------------------------------------------------------------

anchor_records = []

current_lncrna = None
current_upstream = None

with open(anchor_file, "r") as infile:

    for raw_line in infile:

        line = raw_line.strip()

        if not line:
            continue

        if line.startswith("lncRNA:"):
            current_lncrna = line.split(":", 1)[1].strip()

        elif line.startswith("upstream PCG:"):
            anchor_text = line.split(":", 1)[1].strip()
            current_upstream = parse_anchor(anchor_text)

        elif line.startswith("downstream PCG:"):
            anchor_text = line.split(":", 1)[1].strip()
            current_downstream = parse_anchor(anchor_text)

            if current_lncrna is not None:
                anchor_records.append({
                    "lncrna_name": current_lncrna,
                    "upstream": current_upstream,
                    "downstream": current_downstream
                })

            current_lncrna = None
            current_upstream = None


# ------------------------------------------------------------
# MERGE THE TWO DATASETS
# ------------------------------------------------------------

merged_records = []

missing_coordinate_count = 0
multiple_coordinate_count = 0

for anchor_record in anchor_records:

    lncrna_name = anchor_record["lncrna_name"]

    coordinate_matches = lncrna_coordinates.get(lncrna_name, [])

    if len(coordinate_matches) == 0:

        missing_coordinate_count += 1

        merged_records.append({
            "lncrna_name": lncrna_name,
            "lncrna_chr": "",
            "lncrna_start": "",
            "lncrna_end": "",
            "lncrna_strand": "",
            "lncrna_isoform_count": "",
            "upstream": anchor_record["upstream"],
            "downstream": anchor_record["downstream"],
            "coordinate_match_status": "NO_COORDINATE_MATCH"
        })

    elif len(coordinate_matches) == 1:

        coordinate_record = coordinate_matches[0]

        merged_records.append({
            **coordinate_record,
            "upstream": anchor_record["upstream"],
            "downstream": anchor_record["downstream"],
            "coordinate_match_status": "UNIQUE_MATCH"
        })

    else:

        multiple_coordinate_count += 1

        # Write one row for every genomic locus with that name.
        for coordinate_record in coordinate_matches:

            merged_records.append({
                **coordinate_record,
                "upstream": anchor_record["upstream"],
                "downstream": anchor_record["downstream"],
                "coordinate_match_status": "MULTIPLE_NAME_MATCH"
            })


# ------------------------------------------------------------
# WRITE THE FINAL TAB-DELIMITED TABLE
# ------------------------------------------------------------

column_names = [
    "human_lncrna_name",
    "human_lncrna_chr",
    "human_lncrna_start",
    "human_lncrna_end",
    "human_lncrna_strand",
    "human_lncrna_isoform_count",

    "human_upstream_pcg",
    "human_upstream_distance",
    "human_upstream_chr",
    "human_upstream_start",
    "human_upstream_end",
    "human_upstream_strand",

    "human_downstream_pcg",
    "human_downstream_distance",
    "human_downstream_chr",
    "human_downstream_start",
    "human_downstream_end",
    "human_downstream_strand",

    "coordinate_match_status"
]


with open(output_file, "w", newline="") as outfile:

    writer = csv.DictWriter(
        outfile,
        fieldnames=column_names,
        delimiter="\t"
    )

    writer.writeheader()

    for record in merged_records:

        upstream = record["upstream"]
        downstream = record["downstream"]

        writer.writerow({
            "human_lncrna_name": record["lncrna_name"],
            "human_lncrna_chr": record["lncrna_chr"],
            "human_lncrna_start": record["lncrna_start"],
            "human_lncrna_end": record["lncrna_end"],
            "human_lncrna_strand": record["lncrna_strand"],
            "human_lncrna_isoform_count":
                record["lncrna_isoform_count"],

            "human_upstream_pcg": upstream["pcg_name"],
            "human_upstream_distance": upstream["distance"],
            "human_upstream_chr": upstream["chromosome"],
            "human_upstream_start": upstream["start"],
            "human_upstream_end": upstream["end"],
            "human_upstream_strand": upstream["strand"],

            "human_downstream_pcg": downstream["pcg_name"],
            "human_downstream_distance": downstream["distance"],
            "human_downstream_chr": downstream["chromosome"],
            "human_downstream_start": downstream["start"],
            "human_downstream_end": downstream["end"],
            "human_downstream_strand": downstream["strand"],

            "coordinate_match_status":
                record["coordinate_match_status"]
        })


# ------------------------------------------------------------
# PRINT A SUMMARY
# ------------------------------------------------------------

print("Finished merging the human lncRNA files.")
print(f"Coordinate records read: {sum(len(v) for v in lncrna_coordinates.values())}")
print(f"Anchor records read: {len(anchor_records)}")
print(f"Rows written: {len(merged_records)}")
print(f"Anchor records without coordinate matches: {missing_coordinate_count}")
print(f"Names with multiple coordinate matches: {multiple_coordinate_count}")
print(f"Output file: {output_file}")