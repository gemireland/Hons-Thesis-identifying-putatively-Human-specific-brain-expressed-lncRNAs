#this script has been written to find the PCG blocks coords previously extracted and converts them into BED format
input_file = "outputgorillastep1"
output_file = "Gorilla_PCG_Blocks.bed"

f = open(input_file, "r")
w = open(output_file, "w")

#in_pcg_section = False
count = 0
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

        count += 1
        #this is because we originally identified this no. of PCGs
        if count == 21587:
            break


w.close()

print("Gorilla_PCG_Blocks.bed made")
