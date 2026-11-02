input_file = "outputfinalchimpensstep1"
output_file = "Chimp_PCG_Blocks.bed"

f = open(input_file, "r")
w = open(output_file, "w")

in_pcg_section = False

for line in f:

    line = line.strip()

    if line[0]!='[':
        continue


    clean_line = line.replace("[", "").replace("]", "").replace("'", "").replace(",", "")
    fields = clean_line.split()

    if len(fields) >= 5:
        gene = fields[0]
        chrom = fields[1]
        start = fields[2]
        end = fields[3]
        strand = fields[4]

        bed_start = int(start) - 1
        bed_end = int(end)

        w.write(f"{chrom}\t{bed_start}\t{bed_end}\t{gene}\t0\t{strand}\n")


w.close()

print("Finished making Chimp_PCG_Blocks.bed")
