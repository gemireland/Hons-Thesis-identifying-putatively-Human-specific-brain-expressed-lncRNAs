#chimpensstep1.py
#this is very basic code and the script was refined in later scripts

#to read the chimp gtf file into a pandas DataFrame
f = open('Pan_troglodytes.Pan_tro_3.0.116.gtf', 'r')
chimplines = f.readlines() #to read the gtf file into a list of lines, which we can then process to extract the information we need for our analysis
f.close()

print(len(chimplines)) #to check how many lines we have in the gtf file
'''print(chimplines[0:2])''' #to check the structure of the gtf file by looking at the first 2 lines, to see which columns we need to extract for our analysis

print("lncRNA dictionary gene blocks:")

all_biotypes = {}

chimplncRNAs = {'+':{}, '-':{}}
for line in chimplines:

    if line.startswith('#'):
        continue #to skip the header lines that start with #, which are not valid gtf lines


    columns = line.strip().split('\t')

    if len(columns) !=9:
        continue #to skip any lines that don't have 9 columns, which are not valid gtf lines

    chrom = columns[0]
    database = columns[1]
    feature = columns[2]
    start = int(columns[3])
    end = int(columns[4])
    score = columns[5]
    strand = columns[6]
    frame = columns[7]
    attributes = columns[8]



    #keep transcript rows because these represent isoforms
    if feature != "transcript":
        continue

    #added code to count the different types of biotypes
    bio = attributes.split('transcript_biotype')[1].split(';')[0].strip().strip('"')
    try:
        all_biotypes[bio]+=1
    except:
        all_biotypes[bio]=1

    #keep only lncRNA transcripts (also included pseudogenes and mics_RNAs)
    '''lncRNA_biotypes = [
    'gene_biotype "lncRNA"',
    'transcript_biotype "lncRNA"',
    'gene_biotype "lincRNA"',
    'transcript_biotype "lincRNA"',
    'gene_biotype "antisense"',
    'transcript_biotype "antisense"',
    'gene_biotype "processed_transcript"',
    'transcript_biotype "processed_transcript"',
    'pseudogene',
    'processed_pseudogene'
    'misc_RNA']
        
    if not any(biotype in attributes for biotype in lncRNA_biotypes):
        continue'''
    
    lncRNA_biotypes = ['lincRNA','lncRNA','pseudogene','processed_pseudogene','misc_RNA']
    if bio not in lncRNA_biotypes:
        continue

    #the following code needs to be indented to be contained within the for loop (it was not indented)
    gene_id = attributes.split('gene_id "')[1].split('"')[0]
    transcript_id = attributes.split('transcript_id "')[1].split('"')[0]


    #to store the lncRNA transcripts in a dictionary with the gene name as the key and a list of coordinates and transcript ids as the value, separated by strand
    if gene_id not in chimplncRNAs[strand]:
        chimplncRNAs[strand][gene_id] = []

    #to add the coordinates and transcript id of each lncRNA transcript to the list of values for the corresponding gene name in the dictionary, separated by strand
    chimplncRNAs[strand][gene_id].append(
        [chrom, start, end, transcript_id]
        )


total_isoforms = 0

for strand in chimplncRNAs:
    for gene_id in chimplncRNAs[strand]:
        total_isoforms += len(chimplncRNAs[strand][gene_id])
    
print("Total chimp lncRNA genes:", len(chimplncRNAs["+"]) + len(chimplncRNAs["-"]))
print("Total chimp lncRNA isoforms:", total_isoforms)

print(all_biotypes)

#to count isoform numbers per gene 
#for gene_id in chimplncRNAs['+']:
    #print(gene_id, len(chimplncRNAs['+'][gene_id]))


#to see a gene with a lot of isoforms (this was done to check as the script wasn't working)
#for gene_id in chimplncRNAs['+']:
  #print(gene_id)
  #print(chimplncRNAs['+'][gene_id])
   #break

print("Creating lncRNA gene blocks")
#now we want to create gene blocks for the lncRNAs, which are the coordinates that encompass all the isoforms of a gene, so we can compare them to the PCG gene blocks. We will create a dictionary of lncRNA gene blocks keyed by gene name, with the value being a list of coordinates for each isoform of that gene. We will also filter out any lncRNA genes that have isoforms on different chromosomes, as this would make it difficult to compare them to PCG gene blocks.
lncRNAs_gene_blocks = []

for strand in chimplncRNAs:
    for gene_id in chimplncRNAs[strand]:
        isoforms = chimplncRNAs[strand][gene_id]

        chroms = [x[0] for x in isoforms] #go through every isoform and collect the chromosome value
        starts = [x[1] for x in isoforms] #go through every isoform and collect the start value
        ends = [x[2] for x in isoforms]  #go through every isoform and collect the end value

        chrom = chroms[0]
        gene_block_start = min(starts)
        gene_block_end = max(ends)

        lncRNAs_gene_blocks.append([gene_id, chrom, gene_block_start, gene_block_end, strand, len(isoforms)])
       #print(gene_id, chrom, gene_block_start, gene_block_end, len(isoforms)) #to check the gene blocks we have created for the lncRNAs


#now we need to find the nearest PCG to each lncRNA gene block
#make a new dictionary for PCGs, with the same structure as the lncRNA dictionary 
print("PCG dictionary gene blocks:")
pcgsDict = {'+': {}, '-': {}}

for line in chimplines:
    columns = line.strip().split('\t')

    if len(columns) != 9:
        continue

    chrom = columns[0]
    database = columns[1]
    feature = columns[2]
    start = int(columns[3])
    end = int(columns[4])
    score = columns[5]
    strand = columns[6]
    frame = columns[7]
    attributes = columns[8]

    #keep transcript rows because these represent isoforms
    if feature != "transcript":
        continue

    #keep only protein coding transcripts
    if 'gene_biotype "protein_coding"' not in attributes:
        continue

    #extract gene name and transcript ID
    #extract the unique Ensembl gene ID
    gene_id = attributes.split('gene_id "')[1].split('"')[0]

    #extract the gene symbol where available
    if 'gene_name "' in attributes:
        gene_name = attributes.split('gene_name "')[1].split('"')[0]
    else:
        gene_name = gene_id

    transcript_id = attributes.split('transcript_id "')[1].split('"')[0]

    #IMPORTANT:
    # Group transcripts using gene_id because gene_id is unique.
    # Store gene_name separately for the final output and BLAST matching.
    if gene_id not in pcgsDict[strand]:
        pcgsDict[strand][gene_id] = {
            "gene_name": gene_name,
            "isoforms": []
        }

    pcgsDict[strand][gene_id]["isoforms"].append(
        [chrom, start, end, transcript_id]
    )


total_isoforms = 0

for strand in pcgsDict:
    for gene_id in pcgsDict[strand]:
        total_isoforms += len(
            pcgsDict[strand][gene_id]["isoforms"]
        )
print("total pcgs", len(pcgsDict["+"]) + len(pcgsDict["-"]))
print("total pcgs isoforms", total_isoforms)

# count isoform numbers per gene on + strand: COMMENTING OUT AS NOT ESSENTIAL TO CODE AND WANT TO TRACK OUTPUT MORE EASILY
'''for gene_id in pcgsDict['+']:
    print(gene_id, len(pcgsDict['+'][gene_id]))
 
for gene_id in pcgsDict['-']:
    print(gene_id, len(pcgsDict['-'][gene_id]))

#to see a PCG with a lot of isoforms
for gene_id in pcgsDict['+']:
    print(gene_id)
    print(pcgsDict['+'][gene_id])
    break #to see a PCG with a lot of isoforms, we can break after the first one we find'''

#now we want to create gene blocks for the lncRNAs, which are the coordinates that encompass all the isoforms of a gene, so we can compare them to the PCG gene blocks. We will create a dictionary of lncRNA gene blocks keyed by gene name, with the value being a list of coordinates for each isoform of that gene. We will also filter out any lncRNA genes that have isoforms on different chromosomes, as this would make it difficult to compare them to PCG gene blocks.
pcgs_gene_blocks = []

for strand in pcgsDict:
    for gene_id in pcgsDict[strand]:

        gene_name = pcgsDict[strand][gene_id]["gene_name"]
        pcgsisoforms = pcgsDict[strand][gene_id]["isoforms"]

        chroms = [x[0] for x in pcgsisoforms]
        starts = [x[1] for x in pcgsisoforms]
        ends = [x[2] for x in pcgsisoforms]

        # Normally all transcripts belonging to one gene_id
        # should be on the same chromosome.
        chrom = chroms[0]

        gene_block_start = min(starts)
        gene_block_end = max(ends)

        pcgs_gene_blocks.append([
            gene_id,
            gene_name,
            chrom,
            gene_block_start,
            gene_block_end,
            strand,
            len(pcgsisoforms)
        ])

#okay now we want to see the nearest PCG blocks that flank the lncRNA blocks (compare lncRNAs_gene_blocks to pcg_gene_blocks) and the distance to those PCG blocks)

#CREATE A VIARABLE TO STORE THE UPSTREAM AND DOWNSTREAM

anchors = []

for lnc in lncRNAs_gene_blocks:

    lnc_name = lnc[0]
    lnc_chrom = lnc[1]
    lnc_start = lnc[2]
    lnc_end = lnc[3]
    lnc_strand = lnc[4]
    lnc_nisoforms = lnc[5]

    upstream_pcgs = []
    downstream_pcgs = []


#this will print the number of PCGs only for the last lnc_chrom value in the above for loop - this is not essential to the code though, commenting it out for now, AS IT IS NOT INDENTED
#print("lnc chrom:", lnc_chrom)
#print("number of PCGs on this chrom:", len([p for p in pcgs_gene_blocks if p[1] == lnc_chrom]))


#OK THIS PCG PART NEEDS TO BE WITHIN THE ABOVE FOR LOOP, AS YOU WANT TO CHECK THIS FOR EVERY LNCRNA   
    for pcg in pcgs_gene_blocks:

        pcg_gene_id = pcg[0]
        pcg_name = pcg[1]
        pcg_chrom = pcg[2]
        pcg_start = pcg[3]
        pcg_end = pcg[4]
        pcg_strand = pcg[5]
        pcg_nisoforms = pcg[6]

       #only compare the pcg on the same chromosome as the lncRNA blocks
        if pcg_chrom != lnc_chrom:
            continue

        #PCG is upstream if it STARTS before the lncRNA starts
        if pcg_start <= lnc_start:
            distance = lnc_start - pcg_end #its fine if this is negative, will still give the closest PCG
            upstream_pcgs.append((pcg_name, distance, pcg_chrom, pcg_start, pcg_end, pcg_strand))
            
        #PCG is downstream if it starts after the lncRNA ends
        if pcg_start >= lnc_start:
            distance = pcg_start - lnc_end #its fine if this is negative, will still give the closest PCG
            downstream_pcgs.append((pcg_name, distance, pcg_chrom, pcg_start, pcg_end, pcg_strand))
            
    #THIS NEEDS TO OUTDENTED TO OUTISDE THE FOR LOOP. ONLY SORT ONCE YOU HAVE ADDED ALL PCGS TO UPSTREAM AND DOWNSTREAM LISTS
    #sort the upstream and downstream PCGs by distance to the lncRNA block
    upstream_pcgs.sort(key=lambda x: x[1]) #sort by distance, which is the second element #for each row x, use the second element (distance) as the key for sorting
    downstream_pcgs.sort(key=lambda x: x[1]) #sort by distance, which is the second element   

    nearest_upstreamgene = upstream_pcgs[0] if len(upstream_pcgs) > 0  else None
    nearest_downstreamgene = downstream_pcgs[0] if len(downstream_pcgs) > 0 else None

    #NOW WE CAN STORE THIS IN THE VARIABLE ANCHORS
    anchors.append((lnc,nearest_upstreamgene,nearest_downstreamgene))


#TO PRINT EACH LNC AND ITS NEIGHBOURS, LOOP THROUGH ANCHORS

for a in anchors:
    lnc = a[0]
    upstream = a[1]
    downstream = a[2]
    print (lnc,":::", upstream,":::", downstream)

output_file = "Chimp_lncRNA_Full_Anchor_Table.tsv"

with open(output_file, "w") as w:

    w.write(
        "chimp_lncrna_name\t"
        "chimp_lncrna_chr\t"
        "chimp_lncrna_start\t"
        "chimp_lncrna_end\t"
        "chimp_lncrna_strand\t"
        "chimp_lncrna_isoform_count\t"

        "chimp_upstream_pcg\t"
        "chimp_upstream_distance\t"
        "chimp_upstream_chr\t"
        "chimp_upstream_start\t"
        "chimp_upstream_end\t"
        "chimp_upstream_strand\t"

        "chimp_downstream_pcg\t"
        "chimp_downstream_distance\t"
        "chimp_downstream_chr\t"
        "chimp_downstream_start\t"
        "chimp_downstream_end\t"
        "chimp_downstream_strand\n"
    )

    for a in anchors:

        lnc = a[0]
        upstream = a[1]
        downstream = a[2]

        if upstream is None:
            upstream = (
                "", "", "",
                "", "", ""
            )

        if downstream is None:
            downstream = (
                "", "", "",
                "", "", ""
            )

        w.write(
            f"{lnc[0]}\t"
            f"{lnc[1]}\t"
            f"{lnc[2]}\t"
            f"{lnc[3]}\t"
            f"{lnc[4]}\t"
            f"{lnc[5]}\t"

            f"{upstream[0]}\t"
            f"{upstream[1]}\t"
            f"{upstream[2]}\t"
            f"{upstream[3]}\t"
            f"{upstream[4]}\t"
            f"{upstream[5]}\t"

            f"{downstream[0]}\t"
            f"{downstream[1]}\t"
            f"{downstream[2]}\t"
            f"{downstream[3]}\t"
            f"{downstream[4]}\t"
            f"{downstream[5]}\n"
        )

print("made", output_file)
print("total chimp lncRNA blocks written", len(anchors))
