#mousegencodestep1.py
#this script was written in very basic code, I refined it in later scripts 
#to read the mouse gtf file into a pandas DataFrame
f = open('gencode.vM38.annotation.gtf', 'r')
mouselines = f.readlines() #to read the gtf file into a list of lines, which we can then process to extract the information we need for our analysis
f.close()

print(len(mouselines)) #to check how many lines we have in the gtf file
'''print(mouselines[0:2])''' #to check the structure of the gtf file by looking at the first 2 lines, to see which columns we need to extract for our analysis

print("lncRNA dictionary gene blocks:")

mouselncRNAs = {'+':{}, '-':{}}
for line in mouselines:

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

    #keep only lncRNA transcripts
    if 'gene_type "lncRNA"' not in attributes:
        continue 

    #extract gene name from attributes
    gene_name = attributes.split('gene_name "')[1].split('"')[0] #to extract the gene name from the attributes column, which is in the format gene_name "GENE
    transcript_id = attributes.split("transcript_id")[1].split('"')[1] #to extract the transcript id from the attributes column, which is in the format transcript_id "TRANSCRIPT_ID"; we split the string at 'transcript_id "' and take the second part, then split that at '"' and take the first part to get just the transcript id

#to store the lncRNA transcripts in a dictionary with the gene name as the key and a list of coordinates and transcript ids as the value, separated by strand
    if gene_name not in mouselncRNAs[strand]:
        mouselncRNAs[strand][gene_name] = []

    #to add the coordinates and transcript id of each lncRNA transcript to the list of values for the corresponding gene name in the dictionary, separated by strand
    mouselncRNAs[strand][gene_name].append(
        [chrom, start, end, transcript_id]
        )


total_isoforms = 0

for strand in mouselncRNAs:
    for gene_name in mouselncRNAs[strand]:
        total_isoforms += len(mouselncRNAs[strand][gene_name])
    
print("Total mouse lncRNA genes:", len(mouselncRNAs["+"]) + len(mouselncRNAs["-"]))
print("Total mouse lncRNA isoforms:", total_isoforms)

#to count isoform numbers per gene 
for gene_name in mouselncRNAs['+']:
    print(gene_name, len(mouselncRNAs['+'][gene_name]))


#to see a gene with a lot of isoforms
for gene_name in mouselncRNAs['+']:
    print(gene_name)
    print(mouselncRNAs['+'][gene_name])
    break

print("Creating lncRNA gene blocks:")
#now we want to create gene blocks for the lncRNAs, which are the coordinates that encompass all the isoforms of a gene, so we can compare them to the PCG gene blocks. We will create a dictionary of lncRNA gene blocks keyed by gene name, with the value being a list of coordinates for each isoform of that gene. We will also filter out any lncRNA genes that have isoforms on different chromosomes, as this would make it difficult to compare them to PCG gene blocks.
lncRNAs_gene_blocks = []

for strand in mouselncRNAs:
    for gene_name in mouselncRNAs[strand]:
        isoforms = mouselncRNAs[strand][gene_name]

        chroms = [x[0] for x in isoforms] #go through every isoform and collect the chromosome value
        starts = [x[1] for x in isoforms] #go through every isoform and collect the start value
        ends = [x[2] for x in isoforms]  #go through every isoform and collect the end value

        chrom = chroms[0]
        gene_block_start = min(starts)
        gene_block_end = max(ends)

        lncRNAs_gene_blocks.append([gene_name, chrom, gene_block_start, gene_block_end, strand, len(isoforms)])
        print(gene_name, chrom, gene_block_start, gene_block_end, strand, len(isoforms)) #to check the gene blocks we have created for the lncRNAs


#now we need to find the nearest PCG to each lncRNA gene block
#make a new dictionary for PCGs, with the same structure as the lncRNA dictionary 
print("PCG dictionary gene blocks:")
pcgsDict = {'+': {}, '-': {}}

for line in mouselines:
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

    # keep transcript rows because these represent isoforms
    if feature != "transcript":
        continue

    # keep only protein coding transcripts
    if 'gene_type "protein_coding"' not in attributes:
        continue

    # extract gene name and transcript ID
    gene_name = attributes.split('gene_name "')[1].split('"')[0]
    transcript_id = attributes.split('transcript_id "')[1].split('"')[0]

    # store protein coding transcripts
    if gene_name not in pcgsDict[strand]:
        pcgsDict[strand][gene_name] = []

    pcgsDict[strand][gene_name].append([chrom, start, end, transcript_id])


total_isoforms = 0

for strand in pcgsDict:
    for gene_name in pcgsDict[strand]:
        total_isoforms += len(pcgsDict[strand][gene_name])

print("Total pcgs:", len(pcgsDict["+"]) + len(pcgsDict["-"]))
print("Total pcgs isoforms:", total_isoforms)

# count isoform numbers per gene on + strand
for gene_name in pcgsDict['+']:
    print(gene_name, len(pcgsDict['+'][gene_name]))

for gene_name in pcgsDict['-']:
    print(gene_name, len(pcgsDict['-'][gene_name]))

#to see a PCG with a lot of isoforms
for gene_name in pcgsDict['+']:
    print(gene_name)
    print(pcgsDict['+'][gene_name])
    break

#now we want to create gene blocks for the lncRNAs, which are the coordinates that encompass all the isoforms of a gene, so we can compare them to the PCG gene blocks. We will create a dictionary of lncRNA gene blocks keyed by gene name, with the value being a list of coordinates for each isoform of that gene. We will also filter out any lncRNA genes that have isoforms on different chromosomes, as this would make it difficult to compare them to PCG gene blocks.
pcgs_gene_blocks = []

for strand in pcgsDict:
    for gene_name in pcgsDict[strand]:
        pcgsisoforms = pcgsDict[strand][gene_name]

        chroms = [x[0] for x in pcgsisoforms] #go through every isoform and collect the chromosome value
        starts = [x[1] for x in pcgsisoforms] #go through every isoform and collect the start value
        ends = [x[2] for x in pcgsisoforms]  #go through every isoform and collect the end value

        chrom = chroms[0]
        gene_block_start = min(starts)
        gene_block_end = max(ends)

        pcgs_gene_blocks.append([gene_name, chrom, gene_block_start, gene_block_end, strand, len(pcgsisoforms)])
        print(gene_name, chrom, gene_block_start, gene_block_end, strand, len(pcgsisoforms)) 

#okay now we want to see the nearest PCG blocks that flank the lncRNA blocks (compare lncRNAs_gene_blocks to pcg_gene_blocks) and the distance to those PCG blocks)
for lnc in lncRNAs_gene_blocks:

    lnc_name = lnc[0]
    lnc_chrom = lnc[1]
    lnc_start = lnc[2]
    lnc_end = lnc[3]
    lnc_strand = lnc[4]
    lnc_nisoforms = lnc[5]

    upstream_pcgs = []
    downstream_pcgs = []

    for pcg in pcgs_gene_blocks:

        pcg_name = pcg[0]
        pcg_chrom = pcg[1]
        pcg_start = pcg[2]
        pcg_end = pcg[3]
        pcg_strand = pcg[4]
        pcg_nisoforms = pcg[5]

       #only compare the pcg on the same chromosome as the lncRNA blocks
        if pcg_chrom != lnc_chrom:
            continue

        #PCG is upstream if it ends before the lncRNA starts
        if pcg_end <= lnc_end:
            distance = lnc_end - pcg_end
            upstream_pcgs.append((pcg_name, distance, pcg_chrom, pcg_start, pcg_end, pcg_strand))
            
        #PCG is downstream if it starts after the lncRNA ends
        if pcg_start >= lnc_start:
            distance = pcg_start - lnc_start
            downstream_pcgs.append((pcg_name, distance, pcg_chrom, pcg_start, pcg_end, pcg_strand))

        #sort the upstream and downstream PCGs by distance to the lncRNA block
        upstream_pcgs.sort(key=lambda x: x[1]) #sort by distance, which is the second element #for each row x, use the second element (distance) as the key for sorting
        downstream_pcgs.sort(key=lambda x: x[1]) #sort by distance, which is the second element   

        nearest_upstreamgene = upstream_pcgs[0] if len(upstream_pcgs) > 0  else None
        nearest_downstreamgene = downstream_pcgs[0] if len(downstream_pcgs) > 0 else None

        upstream_pcgs.sort(key=lambda x: x[1])
        downstream_pcgs.sort(key=lambda x: x[1])

    nearest_upstream = upstream_pcgs[0] if len(upstream_pcgs) > 0 else None
    nearest_downstream = downstream_pcgs[0] if len(downstream_pcgs) > 0 else None

    print("lncRNA:", lnc_name)
    print("upstream PCG:", nearest_upstream)
    print("downstream PCG:", nearest_downstream)
    print()

#save mouse PCG blocks

w = open("Mouse_PCG_Blocks.bed", "w")

for pcg in pcgs_gene_blocks:

    gene_name = pcg[0]
    chrom = pcg[1]
    start = pcg[2]
    end = pcg[3]
    strand = pcg[4]

    w.write(
        f"{chrom}\t{start}\t{end}\t{gene_name}\t0\t{strand}\n"
    )

w.close()

print("Mouse PCG blocks saved")
print("Total PCG blocks:", len(pcgs_gene_blocks))

