#I want to extract the lncRNA coordinates and the PCG coordinates for chimp from the first output file I created

import ast
import pandas as pd


# File containing the printed lncRNA blocks and their PCG anchors
input_file = "outputgorillastep1"

# New full gorilla anchor table
output_file = "Gorilla_lncRNA_Full_Anchor_Table.tsv"


# Store each completed table row
rows = []


# Open the Step 1 output file
with open(input_file, "r") as infile:

    for line in infile:

        # Only use lines containing an lncRNA and its two anchors
        if ":::" not in line:
            continue

        # Divide the line into:
        # lncRNA block ::: upstream PCG ::: downstream PCG
        parts = line.strip().split(":::")

        lnc = ast.literal_eval(parts[0].strip())

        upstream_text = parts[1].strip()
        downstream_text = parts[2].strip()

        # Convert the upstream PCG text into a Python tuple
        if upstream_text == "None":
            upstream = ("NA", "NA", "NA", "NA", "NA", "NA")
        else:
            upstream = ast.literal_eval(upstream_text)

        # Convert the downstream PCG text into a Python tuple
        if downstream_text == "None":
            downstream = ("NA", "NA", "NA", "NA", "NA", "NA")
        else:
            downstream = ast.literal_eval(downstream_text)

        # Combine the lncRNA and its two PCG anchors into one row
        rows.append([
            lnc[0],          # Gorilla lncRNA gene ID
            lnc[1],          # lncRNA chromosome
            lnc[2],          # lncRNA start
            lnc[3],          # lncRNA end
            lnc[4],          # lncRNA strand
            lnc[5],          # Number of lncRNA isoforms

            upstream[0],     # Upstream PCG gene ID
            upstream[1],     # Distance to upstream PCG
            upstream[2],     # Upstream PCG chromosome
            upstream[3],     # Upstream PCG start
            upstream[4],     # Upstream PCG end
            upstream[5],     # Upstream PCG strand

            downstream[0],   # Downstream PCG gene ID
            downstream[1],   # Distance to downstream PCG
            downstream[2],   # Downstream PCG chromosome
            downstream[3],   # Downstream PCG start
            downstream[4],   # Downstream PCG end
            downstream[5]    # Downstream PCG strand
        ])


# Create the full table
df = pd.DataFrame(
    rows,
    columns=[
        "gorilla_lncrna_name",
        "gorilla_lncrna_chr",
        "gorilla_lncrna_start",
        "gorilla_lncrna_end",
        "gorilla_lncrna_strand",
        "gorilla_lncrna_isoform_count",

        "gorilla_upstream_pcg",
        "gorilla_upstream_distance",
        "gorilla_upstream_chr",
        "gorilla_upstream_start",
        "gorilla_upstream_end",
        "gorilla_upstream_strand",

        "gorilla_downstream_pcg",
        "gorilla_downstream_distance",
        "gorilla_downstream_chr",
        "gorilla_downstream_start",
        "gorilla_downstream_end",
        "gorilla_downstream_strand"
    ]
)


# Gorilla has chromosomes 2A and 2B, like chimpanzee
standard_chromosomes = {
    "1", "2A", "2B", "3", "4", "5", "6", "7", "8", "9",
    "10", "11", "12", "13", "14", "15", "16", "17", "18",
    "19", "20", "21", "22", "X", "Y", "MT"
}


# Add "chr" only to standard chromosomes.
# Scaffold names such as CABD030... remain unchanged.
for column in [
    "gorilla_lncrna_chr",
    "gorilla_upstream_chr",
    "gorilla_downstream_chr"
]:

    df[column] = df[column].apply(
        lambda value: (
            f"chr{value}"
            if (
                pd.notna(value)
                and str(value) in standard_chromosomes
            )
            else value
        )
    )


# This table is not performing coordinate matching,
# so do not incorrectly label every row UNIQUE_MATCH.
df["coordinate_match_status"] = "NOT_APPLICABLE"


# Save as a tab-separated table
df.to_csv(output_file, sep="\t", index=False)


print(f"Finished: {output_file}")
print(f"Total gorilla lncRNAs written: {len(df)}")
print(
    "Total gorilla lncRNA isoforms:",
    pd.to_numeric(
        df["gorilla_lncrna_isoform_count"],
        errors="coerce"
    ).sum()
)