#I want to extract the lncRNA coordinates and the PCG coordinates for chimp from the first output file I created

#Importing pandas libarary to read, edit and save the table 
import pandas as pd

#this input file has existing chimp PCG anchors 
input_file = "Chimp_anchors_output.txt"

#this output file will be a tab delimited file that has the lncRNAs and the PCGs 
output_file = "Chimp_lncRNA_full_Anchor_Table.tsv"

#i want to read the input into a pandas (pd) dataframe (df) as a dataframe tabel has rows and columns
df = pd.read_csv(input_file, sep="\t")

#creating columns names that are chimp specific 
df.columns =[
    "chimp_lncrna_name",             # Chimp lncRNA gene ID
    "chimp_lncrna_chr",              # Chromosome containing the lncRNA
    "chimp_lncrna_start",            # Start coordinate of the lncRNA
    "chimp_lncrna_end",              # End coordinate of the lncRNA
    "chimp_lncrna_strand",           # Strand of the lncRNA: + or -
    "chimp_lncrna_isoform_count",    # Number of transcript isoforms for the lncRNA
    "chimp_upstream_pcg",            # Gene ID of the upstream PCG anchor
    "chimp_upstream_distance",       # Distance between the lncRNA and upstream PCG
    "chimp_upstream_chr",            # Chromosome containing the upstream PCG
    "chimp_upstream_start",          # Start coordinate of the upstream PCG
    "chimp_upstream_end",            # End coordinate of the upstream PCG
    "chimp_upstream_strand",         # Strand of the upstream PCG
    "chimp_downstream_pcg",          # Gene ID of the downstream PCG anchor
    "chimp_downstream_distance",     # Distance between the lncRNA and downstream PCG
    "chimp_downstream_chr",          # Chromosome containing the downstream PCG
    "chimp_downstream_start",        # Start coordinate of the downstream PCG
    "chimp_downstream_end",          # End coordinate of the downstream PCG
    "chimp_downstream_strand"        # Strand of the downstream PCG
]

#we need to go through each chr column and ensure that there is "chr" before the number
standard_chromosomes = {
     "1", "2A", "2B", "3", "4", "5", "6", "7", "8", "9",
    "10", "11", "12", "13", "14", "15", "16", "17", "18",
    "19", "20", "21", "22", "X", "Y", "MT"
}

for column in [
    "chimp_lncrna_chr",
    "chimp_upstream_chr",
    "chimp_downstream_chr"
]:
    df[column] = df[column].apply(
        lambda value: (
            f"chr{value}"
            if pd.notna(value) and str(value) in standard_chromosomes
            else value
        )
    )

#add a new column called coord_match_status 
# Count how many rows have the same lncRNA chromosome, start, end and strand

coordinate_counts = df.groupby([
    "chimp_lncrna_chr",
    "chimp_lncrna_start",
    "chimp_lncrna_end",
    "chimp_lncrna_strand"
])["chimp_lncrna_name"].transform("count")

# Label coordinates that occur once as UNIQUE_MATCH.
# Label coordinates that occur more than once as MULTIPLE_MATCH.

df["coordinate_match_status"] = coordinate_counts.apply(
    lambda count: "UNIQUE_MATCH" if count == 1 else "MULTIPLE_MATCH"
)

#save the completed dataframe as a tsv file  (index=Flase prevents pandas adding an extra numbered column)
df.to_csv(output_file, sep="\t", index = False)

print(f"Finihsed: {output_file}")
print(f"Total chimp lncRNAs written: {len(df)}")