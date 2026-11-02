#gorillastep1.py

#to read the gorilla gtf file into a pandas DataFrame
f = open('Gorilla_gorilla.gorGor4.116.gtf', 'r')
gorillalines = f.readlines() #to read the gtf file into a list of lines, which we can then process to extract the information we need for our analysis
f.close()

print(len(gorillalines)) #to check how many lines we have in the gtf file
'''print(gorillalines[0:2])''' #to check the structure of the gtf file by looking at the first 2 lines, to see which columns we need to extract for our analysis

print("lncRNA dictionary gene blocks:")

all_biotypes = {}

gorillalncRNAs = {'+':{}, '-':{}}
for line in gorillalines:

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
    if gene_id not in gorillalncRNAs[strand]:
        gorillalncRNAs[strand][gene_id] = []

    #to add the coordinates and transcript id of each lncRNA transcript to the list of values for the corresponding gene name in the dictionary, separated by strand
    gorillalncRNAs[strand][gene_id].append(
        [chrom, start, end, transcript_id]
        )


total_isoforms = 0

for strand in gorillalncRNAs:
    for gene_id in gorillalncRNAs[strand]:
        total_isoforms += len(gorillalncRNAs[strand][gene_id])
    
print("Total gorilla lncRNA genes:", len(gorillalncRNAs["+"]) + len(gorillalncRNAs["-"]))
print("Total gorilla lncRNA isoforms:", total_isoforms)

print(all_biotypes)

#to count isoform numbers per gene 
#for gene_id in gorillalncRNAs['+']:
    #print(gene_id, len(gorillalncRNAs['+'][gene_id]))


#to see a gene with a lot of isoforms
#for gene_id in gorillalncRNAs['+']:
  #  print(gene_id)
  #  print(gorillalncRNAs['+'][gene_id])
   # break

print("Creating lncRNA gene blocks:")
#now we want to create gene blocks for the lncRNAs, which are the coordinates that encompass all the isoforms of a gene, so we can compare them to the PCG gene blocks. We will create a dictionary of lncRNA gene blocks keyed by gene name, with the value being a list of coordinates for each isoform of that gene. We will also filter out any lncRNA genes that have isoforms on different chromosomes, as this would make it difficult to compare them to PCG gene blocks.
lncRNAs_gene_blocks = []

for strand in gorillalncRNAs:
    for gene_id in gorillalncRNAs[strand]:
        isoforms = gorillalncRNAs[strand][gene_id]

        chroms = [x[0] for x in isoforms] #go through every isoform and collect the chromosome value
        starts = [x[1] for x in isoforms] #go through every isoform and collect the start value
        ends = [x[2] for x in isoforms]  #go through every isoform and collect the end value

        chrom = chroms[0]
        gene_block_start = min(starts)
        gene_block_end = max(ends)

        lncRNAs_gene_blocks.append([gene_id, chrom, gene_block_start, gene_block_end, strand, len(isoforms)])
       # print(gene_id, chrom, gene_block_start, gene_block_end, len(isoforms)) #to check the gene blocks we have created for the lncRNAs


#now we need to find the nearest PCG to each lncRNA gene block
#make a new dictionary for PCGs, with the same structure as the lncRNA dictionary 
print("PCG dictionary gene blocks:")
pcgsDict = {'+': {}, '-': {}}

for line in gorillalines:
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
    if 'gene_biotype "protein_coding"' not in attributes:
        continue

    # extract gene name and transcript ID
    gene_id = attributes.split('gene_id "')[1].split('"')[0]
    transcript_id = attributes.split('transcript_id "')[1].split('"')[0]


    # store protein coding transcripts
    if gene_id not in pcgsDict[strand]:
        pcgsDict[strand][gene_id] = []

    pcgsDict[strand][gene_id].append([chrom, start, end, transcript_id])


total_isoforms = 0

for strand in pcgsDict:
    for gene_id in pcgsDict[strand]:
        total_isoforms += len(pcgsDict[strand][gene_id])

print("Total pcgs:", len(pcgsDict["+"]) + len(pcgsDict["-"]))
print("Total pcgs isoforms:", total_isoforms)

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
        pcgsisoforms = pcgsDict[strand][gene_id]

        chroms = [x[0] for x in pcgsisoforms] #go through every isoform and collect the chromosome value
        starts = [x[1] for x in pcgsisoforms] #go through every isoform and collect the start value
        ends = [x[2] for x in pcgsisoforms]  #go through every isoform and collect the end value

        chrom = chroms[0]
        gene_block_start = min(starts)
        gene_block_end = max(ends)

        pcgs_gene_blocks.append([gene_id, chrom, gene_block_start, gene_block_end, strand, len(pcgsisoforms)])
        print([gene_id, chrom, gene_block_start, gene_block_end, len(pcgsisoforms)])


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

        pcg_name = pcg[0]
        pcg_chrom = pcg[1]
        pcg_start = pcg[2]
        pcg_end = pcg[3]
        pcg_strand = pcg[4]
        pcg_nisoforms = pcg[5]

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

        #

#here you have repetitive code so I am commenting it out - note the indentations are also wrong here, so you can delete this once you understand the error
'''upstream_pcgs.sort(key=lambda x: x[1])
downstream_pcgs.sort(key=lambda x: x[1])

nearest_upstream = upstream_pcgs[0] if len(upstream_pcgs) > 0 else None
nearest_downstream = downstream_pcgs[0] if len(downstream_pcgs) > 0 else None''' 

'''print("lncRNA:", lnc_name)
print("upstream PCG:", nearest_upstream)
print("downstream PCG:", nearest_downstream)
print()'''

#TO PRINT EACH LNC AND ITS NEIGHBOURS YOU CAN NOW LOOP THROUGH ANCHORS

for a in anchors:
    lnc = a[0]
    upstream = a[1]
    downstream = a[2]
    print (lnc,":::", upstream,":::", downstream)

w = open("gorilla_PCG_anchors.txt", "w")

w.write("upstream_pcg\tupstream_distance\tupstream_chr\tupstream_start\tupstream_end\tupstream_strand\t")
w.write("downstream_pcg\tdownstream_distance\tdownstream_chr\tdownstream_start\tdownstream_end\tdownstream_strand\n")

for a in anchors:

    upstream = a[1]
    downstream = a[2]

    if upstream is None:
        upstream = ("NA","NA","NA","NA","NA","NA")

    if downstream is None:
        downstream = ("NA","NA","NA","NA","NA","NA")

    w.write(
        f"{upstream[0]}\t{upstream[1]}\t{upstream[2]}\t{upstream[3]}\t{upstream[4]}\t{upstream[5]}\t"
        f"{downstream[0]}\t{downstream[1]}\t{downstream[2]}\t{downstream[3]}\t{downstream[4]}\t{downstream[5]}\n"
    )

w.close()