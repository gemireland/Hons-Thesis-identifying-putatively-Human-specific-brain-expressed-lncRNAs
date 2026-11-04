#this script has been written to find the PCG blocks coords previously extracted and converts them into BED format
input_file = "outputhumangencodestep1.txt"
output_file = "Human_PCG_Blocks.bed"

f = open(input_file, "r")
w = open(output_file, "w")

in_pcg_section = False

for line in f:

    line = line.strip()

    if line.startswith("total pcgs isoforms"):
        in_pcg_section = True
        continue

    if line.startswith("lncRNA:"):
        in_pcg_section = False
        break

    if in_pcg_section:
        fields = line.split()

        if len(fields) >= 5:
            gene = fields[0]
            chrom = fields[1]
            start = fields[2]
            end = fields[3]
            strand = fields[4]
            

            if chrom.startswith("chr"):
                bed_start = int(start) - 1
                bed_end = int(end)

                w.write(f"{chrom}\t{bed_start}\t{bed_end}\t{gene}\t0\t{strand}\n")

f.close()
w.close()

print("Human_PCG_Blocks.bed made")

