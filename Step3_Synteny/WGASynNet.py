#Script to process WGA SenNet.Maf
synfile = "hg38.mm39.synNet.maf"
species = 'mm39'
f = open(synfile,'r')
lines = f.readlines()
f.close()

alignedRegions = {}
count=0
countS=0
for n,line in enumerate(lines):
    if 'a score=' in line:

        hregion = ''
        hchr = ''
        hstrand = ''
        hstart = 0
        rlength = 0
        hend = 0
        human_line = lines[n+1]



        
        if 'hg38.' not in human_line:
            count+=1
        else:
            hregion = human_line.split()[1:6]
            hchr = hregion[0].split('.')[1]
            hstrand = hregion[3]
            hstart = int(hregion[1])
            rlength = int(hregion[2]) 
            hend = hstart+rlength
            clength = int(hregion[4])
            if hstrand =='-':
                hend = clength-hstart
                hstart = hend - rlength
            

        sregion = ''
        schr = ''
        sstrand = ''
        sstart = 0
        slength = 0
        send = 0   
        species_line = lines[n+2]
        
        if species+'.' not in species_line:
            countS+=1
        else:
            sregion = species_line.split()[1:6]
            schr = sregion[0].split('.')[1]
            sstrand = sregion[3]
            sstart = int(sregion[1])
            rlength = int(sregion[2]) 
            send = sstart+rlength
            clength = int(sregion[4])
            if sstrand =='-':                           #this si where i have inverted the coordinate sif on teh negative strand
                send = clength-sstart
                sstart = send - rlength


        try:
            alignedRegions[hchr][(hstart,hend)].append([hstrand,schr,sstart,send,sstrand])
        except:
            try:
                alignedRegions[hchr][(hstart,hend)]= [[hstrand,schr,sstart,send,sstrand]]
            except:
                alignedRegions[hchr] = {(hstart,hend):[[hstrand,schr,sstart,send,sstrand]]}
            

#lets check for multiple-to-multiple alignments

countmulti = 0
for chrm in alignedRegions:
    for region in alignedRegions[chrm]:
        if len(alignedRegions[chrm][region])>1:
            countmulti+=1
if countmulti>0:
    print(countmulti)
else:
    print("WGAs do not have multi-multi mapped regions")


#At this stage we have aligned conserved regions from WGA, we will now add additional aligned regions from blast results



#this is a script to process the top blast hits into the same format as the WGA synteny scripts
import csv
species = "mouse"
blastfile = "Human_top_hits_Mouse_Top1.tsv"

########################################################Change filenames above as needed

f = open(blastfile, "r", newline="")
reader = csv.DictReader(f, delimiter="\t")
blastRegions = {}
count = 0
countS =0 #s = species

for row in reader:

    human_pcg_name = row["qseqid"].split("::")[0].strip()
    mouse_pcg_name = row["sseqid"].split("::")[0].strip()


    #HUMAN REGION that looks the same format as the human_lines in the WGA synteny

    hregion =""
    hchr =""
    hstrand =""
    hstart =0
    rlength=0
    hend=0

    if (
        row["query_chromosome"].strip()==""
        or row["query_gene_start"].strip()==""
        or row["query_gene_end"].strip()==""
        or row["query_strand"].strip()==""  
    ):
        count += 1

    else:
        hregion =[
            row["query_chromosome"],
            row["query_gene_start"],
            row["query_gene_end"],
            row["query_strand"]
        ]

        hchr = hregion[0].strip()

        #if the chr doesnt start with "chr" then add it to the front of the string
        if not hchr.startswith("chr"):
            hchr = "chr" + hchr

        hstart = int(hregion[1])
        hend = int(hregion[2])
        hstrand = hregion[3].strip()

        rlength = hend - hstart 

    #SPECIES REGION that looks the same format as the species_lines in the WGA synteny

    sregion =""
    schr =""
    sstrand =""
    sstart =0
    rlength=0
    send=0

    if (
        row["subject_chromosome"].strip()==""
        or row["subject_gene_start"].strip()==""
        or row["subject_gene_end"].strip()==""
        or row["subject_strand"].strip()==""
    ):
        countS += 1

    else:
        sregion =[
            row["subject_chromosome"],
            row["subject_gene_start"],
            row["subject_gene_end"],
            row["subject_strand"]
        ]

        schr = sregion[0].strip()

        #if the chr doesnt start with "chr" then add it to the front of the string
        if not schr.startswith("chr"):
            schr = "chr" + schr

        sstart = int(sregion[1])
        send = int(sregion[2])
        sstrand = sregion[3].strip()

        rlength = send - sstart

    #skip rows where coords are missing for either species
    if hregion == "" or sregion == "":
        #print(hregion,sregion)
        continue    

    #EXACT SAME STRUCTURE AS THE WGA 
    try:
        blastRegions[hchr][(hstart,hend)].append(
            [hstrand,schr,sstart,send,sstrand, human_pcg_name, mouse_pcg_name]
        )
    except:
        try: 
            blastRegions[hchr][(hstart,hend)]=[
                [hstrand,schr,sstart,send,sstrand, human_pcg_name, mouse_pcg_name]
            ]

        except:
            blastRegions[hchr]= {
                (hstart,hend): [
                    [hstrand,schr,sstart,send,sstrand, human_pcg_name, mouse_pcg_name]
                ]
            }

f.close()

countmulti = 0

for chrm in blastRegions:
    for region in blastRegions[chrm]:
        if len(blastRegions[chrm][region]) > 1:
            print(region,blastRegions[chrm][region]) 
            countmulti += 1

print(countmulti)

output_file = "Human_Mouse_BLAST_Aligned_Regions.tsv"              #is this writing the correct thing?? or do I need to move it elsewhere?

with open(output_file, "w", newline="") as outfile:

    writer = csv.writer(outfile, delimiter="\t")

    writer.writerow([
        "human_chromosome",
        "human_start",
        "human_end",
        "human_strand",
        "species_chromosome",
        "species_start",
        "species_end",
        "species_strand"
    ])

    for hchr in blastRegions:

        for hregion in blastRegions[hchr]:

            hstart = hregion[0]
            hend = hregion[1]

            for match in blastRegions[hchr][hregion]:

                hstrand = match[0]
                schr = match[1]
                sstart = match[2]
                send = match[3]
                sstrand = match[4]

                writer.writerow([
                    hchr,
                    hstart,
                    hend,
                    hstrand,
                    schr,
                    sstart,
                    send,
                    sstrand
                ])

print("Output saved to:", output_file)

#####READING IN THE MOUSE AND HUMAN lncRNA BLCOKS

human_lncRNA_file = "Human_lncRNA_blocks.tsv"
mouse_lncRNA_file = "Mouse_lncRNA_blocks.tsv"

def read_lncRNA_blocks(filename):
    ###I want to read in each of the files and sort the lncRNAs by chromosome
    '''
    this will allow me hopefully to sort them and get an output liek:

    lncRNAs["chr1"] =[
        {
        "name": "  ",
        "gene_id": " ",
        "start": "",
        "end": "",
        "strand": "+/-"
        }
    ]    '''

    lncRNAs = {}

    with open(filename, "r", newline ="") as infile:      #newline= "" means that it won't change teh format of teh ends of the lines from what they were in teh file

        reader = csv.DictReader(infile, delimiter="\t")

        for row in reader:

            chromosome = row["chromosome"].strip()  #remove the extra space
            start = int(row["start"])
            end = int(row["end"])
            strand = row["strand"].strip()
            name = row["lncRNA_name"].strip()
            gene_id = row["lncrna_gene_id"].strip()

            #now we need to ensure that we have teh prefix "chr" before every chromosome value
            if not chromosome.startswith("chr"):
                chromosome = "chr" + chromosome

            lncRNA_record = {
                "name": name,
                "gene_id": gene_id,
                "start":start,
                "end": end,
                "strand": strand
            }

            if chromosome not in lncRNAs:
                lncRNAs[chromosome] = []

            lncRNAs[chromosome].append(lncRNA_record)

    #to sort the lncRNAs on each chr by position of their start
    for chromosome in lncRNAs:
        lncRNAs[chromosome].sort(
            key = lambda lncRNA: lncRNA["start"]
        )
      # print lncRNAs[chromosome]

    return lncRNAs

humanLncRNAs = read_lncRNA_blocks(human_lncRNA_file)
mouseLncRNAs = read_lncRNA_blocks(mouse_lncRNA_file)

#I AM CREATING A COMBINED DICTIONARY FOR THE BLAST AND THE WGA ALIGNED REGIONS
totalAlignedRegions = {}

#### alignedRegions.keys will return all of the chr values as these are keys in the dictionary
all_chromosomes = set(alignedRegions.keys()) | set(blastRegions.keys())

for chromosome in all_chromosomes:

    blast_chromosome_regions = blastRegions.get(chromosome, {})           #I used .get instead of blastRegions["chr1"] as .get will return an empty dictionary instead of throwing an error if the chromosome number is misisng
    wga_chromosome_regions = alignedRegions.get(chromosome, {})

    #I am creating an empty list for this chromosome
    totalAlignedRegions[chromosome] = []

    #add the WGA regions
    for human_region in wga_chromosome_regions:

        hstart = human_region[0]
        hend = human_region[1]

        for match in wga_chromosome_regions[human_region]:

            totalAlignedRegions[chromosome].append(
                {
                    "human_start": hstart,
                    "human_end": hend,
                    "human_strand": match[0],
                    "species_chromosome": match[1],
                    "species_start": match[2],
                    "species_end": match[3],
                    "species_strand": match[4],
                    "source": "WGA"
                }
            )

    #add the blast regions
    for human_region in blast_chromosome_regions:                 #looking for human_regions as the WGA is indexed by human_coordinates

        hstart = human_region[0]
        hend = human_region[1]

        for match in blast_chromosome_regions[human_region]:

            totalAlignedRegions[chromosome].append(
               {
                    "human_start": hstart,
                    "human_end": hend,
                    "human_strand": match[0],

                    "species_chromosome": match[1],
                    "species_start": match[2],
                    "species_end": match[3],
                    "species_strand": match[4],

                    "human_pcg_name": match[5],
                    "species_pcg_name": match[6],

                    "source": "BLAST"
                } 
            )

    #now i am going to sort all of the matched regions on the chromosome by the start position of human
    totalAlignedRegions[chromosome].sort(
        key=lambda region: region["human_start"]                            #means take one dictionary called region and return its human_start value
    )

#####SANITY CHECK 
print("\nCombined aligned regions:")
selectedChrm = ['chr1','chr2','chr3','chr4','chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
for chromosome in sorted(totalAlignedRegions):
    if chromosome in selectedChrm:
        print(
            chromosome,
            len(totalAlignedRegions[chromosome])
    )

########NOW I AM STARTING TO LOOK THROUGH THE DICTIONARIES TO MAKE A TEXT FILE TO FIND EACH LNCRNA WITH ITS NEAREST PCG AND ITS NEAREST UPSTREAM AND DOWNSTREAM ALIGNED REGIONS


import bisect        #bisect means that python will go into the values at a value near what it thinks and will move backwards towards it like when you open a book in the middle and page back to page 1


def distance_between_regions(
    first_start,
    first_end,
    second_start,
    second_end
):
    """
    Return the distance between two genomic regions.

    If the regions overlap, the distance is 0.
    """

    # teh second region is upstream of the first region
    if second_end < first_start:
        return first_start - second_end

    # the second region is downstream of the first region
    if second_start > first_end:
        return second_start - first_end

    # if the regionsegions overlap
    return 0


def relative_position(
    lnc_start,
    lnc_end,
    region_start,
    region_end
):
    """
    ######this will describe the position of a region relative to a lncRNA.
    """

    if region_end < lnc_start:
        return "UPSTREAM"

    if region_start > lnc_end:
        return "DOWNSTREAM"

    return "OVERLAPS"

def transcription_relative_to_anchor(
        lnc_start,
        lnc_end,
        lnc_strand,
        region_start,
        region_end
):

    #to determine if the lncRNA transcribes away from or towards the conserved anchor (be it: PCG, NUA or NDA)

    #if the conserved region lies on teh left of the lncRNA
    if region_end <= lnc_start:

        if lnc_strand == "-":
            return "TOWARD"
        else:
            return "AWAY_FROM"


        #if the conserved region lies on the right of the lncRNA
    elif region_start >= lnc_end:

        if lnc_strand == "+":
            return "TOWARD"
        else:
            return "AWAY_FROM"

    #conserved region overlaps the lncRNA
    return "OVERLAPS"

# This list will contain one result dictionary for every lncRNA
lncRNA_neighbour_results = []

conservedlncrnas = {}
# Process every chromosome containing human lncRNAs
for chromosome in humanLncRNAs:

    chromosome_lncRNAs = humanLncRNAs[chromosome]

    # Get all combined WGA and BLAST regions on this chromosome
    chromosome_regions = totalAlignedRegions.get(
        chromosome,
        []
    )

    
    #### WE ARE GETTING THE ALIGNED REGIONS


    #These are used to find the closest downstream region
    regions_sorted_by_start = sorted(
        chromosome_regions,
        key=lambda region: region["human_start"]
    )

    region_start_values = [
        region["human_start"]
        for region in regions_sorted_by_start
    ]

    #These are used to find the closest upstream region
    regions_sorted_by_end = sorted(
        chromosome_regions,
        key=lambda region: region["human_end"]
    )

    region_end_values = [
        region["human_end"]
        for region in regions_sorted_by_end
    ]

    # Keep only BLAST records because these represent PCGs
    chromosome_pcg_regions = [
        region
        for region in chromosome_regions
        if region["source"] == "BLAST"
    ]

    
    # PROCESS EVERY lncRNA ON THIS CHROMOSOME
    

    for lncRNA in chromosome_lncRNAs:

        lnc_start = lncRNA["start"]
        lnc_end = lncRNA["end"]

        
        # FIND THE CLOSEST BLAST PCG
        

        closest_pcg = None
        closest_pcg_distance = None
        closest_pcg_position = ""

        for pcg_region in chromosome_pcg_regions:

            pcg_distance = distance_between_regions(
                lnc_start,
                lnc_end,
                pcg_region["human_start"],
                pcg_region["human_end"]
            )

            if (
                closest_pcg_distance is None
                or pcg_distance < closest_pcg_distance
            ):
                closest_pcg = pcg_region
                closest_pcg_distance = pcg_distance

                closest_pcg_position = relative_position(
                    lnc_start,
                    lnc_end,
                    pcg_region["human_start"],
                    pcg_region["human_end"]
                )

                # Distance 0 is the smallest possible distance
                if pcg_distance == 0:
                    break

        
        # FIND THE CLOSEST UPSTREAM ALIGNED REGION
        

        closest_upstream_region = None
        closest_upstream_distance = None

        # Find the final region whose end is <= lncRNA start
        upstream_index = bisect.bisect_right(
            region_end_values,
            lnc_start
        ) - 1

        if upstream_index >= 0:

            closest_upstream_region = regions_sorted_by_end[
                upstream_index
            ]

            closest_upstream_distance = (
                lnc_start
                - closest_upstream_region["human_end"]
            )

        #  FIND THE CLOSEST DOWNSTREAM ALIGNED REGION
        

        closest_downstream_region = None
        closest_downstream_distance = None

        # Find the first region whose start is >= lncRNA end
        downstream_index = bisect.bisect_left(
            region_start_values,
            lnc_end
        )

        if downstream_index < len(regions_sorted_by_start):

            closest_downstream_region = regions_sorted_by_start[
                downstream_index
            ]

            closest_downstream_distance = (
                closest_downstream_region["human_start"]
                - lnc_end
            )

       
        # 4. CREATE THE OUTPUT RECORD
    

        result = {
            "human_chromosome": chromosome,
            "lncRNA_name": lncRNA["name"],
            "lncRNA_gene_id": lncRNA["gene_id"],
            "lncRNA_start": lnc_start,
            "lncRNA_end": lnc_end,
            "lncRNA_strand": lncRNA["strand"],

            # Closest PCG fields
            "closest_PCG_human_name": "",
            "closest_PCG_start": "",
            "closest_PCG_end": "",
            "closest_PCG_strand": "",
            "closest_PCG_position": "",
            "closest_PCG_distance": "",
            "closest_PCG_mouse_name": "",
            "closest_PCG_mouse_chromosome": "",
            "closest_PCG_mouse_start": "",
            "closest_PCG_mouse_end": "",
            "closest_PCG_mouse_strand": "",

            # Upstream aligned-region fields
            "upstream_region_start": "",
            "upstream_region_end": "",
            "upstream_region_human_strand": "",
            "upstream_region_source": "",
            "upstream_region_distance": "",
            "upstream_region_mouse_chromosome": "",
            "upstream_region_mouse_start": "",
            "upstream_region_mouse_end": "",
            "upstream_region_mouse_strand": "",

            # Downstream aligned-region fields
            "downstream_region_start": "",
            "downstream_region_end": "",
            "downstream_region_human_strand": "",
            "downstream_region_source": "",
            "downstream_region_distance": "",
            "downstream_region_mouse_chromosome": "",
            "downstream_region_mouse_start": "",
            "downstream_region_mouse_end": "",
            "downstream_region_mouse_strand": "",

            "selected_mouse_lncRNA_name": "",
            "selected_mouse_lncRNA_gene_id": "",
            "selected_mouse_lncRNA_chromosome": "",
            "selected_mouse_lncRNA_start": "",
            "selected_mouse_lncRNA_end": "",
            "selected_mouse_lncRNA_strand": "",
            "selected_mouse_lncRNA_anchor": "",
            "selected_mouse_lncRNA_distance": "",
            "selected_mouse_lncRNA_orientation_conserved": ""

        }

        
        # ADD THE CLOSEST PCG INFORMATION

        if closest_pcg is not None:

            result["closest_PCG_human_name"] = closest_pcg.get(
                "human_pcg_name",
                ""
            )

            result["closest_PCG_start"] = closest_pcg[
                "human_start"
            ]

            result["closest_PCG_end"] = closest_pcg[
                "human_end"
            ]

            result["closest_PCG_strand"] = closest_pcg[
                "human_strand"
            ]

            result["closest_PCG_position"] = (
                closest_pcg_position
            )

            result["closest_PCG_distance"] = (
                closest_pcg_distance
            )

            result["closest_PCG_mouse_name"] = closest_pcg.get(
                "species_pcg_name",
                ""
            )

            result["closest_PCG_mouse_chromosome"] = (
                closest_pcg["species_chromosome"]
            )

            result["closest_PCG_mouse_start"] = closest_pcg[
                "species_start"
            ]

            result["closest_PCG_mouse_end"] = closest_pcg[
                "species_end"
            ]

            result["closest_PCG_mouse_strand"] = closest_pcg[
                "species_strand"
            ]

        
        # ADD THE UPSTREAM REGION INFORMATION

        if closest_upstream_region is not None:

            result["upstream_region_start"] = (
                closest_upstream_region["human_start"]
            )

            result["upstream_region_end"] = (
                closest_upstream_region["human_end"]
            )

            result["upstream_region_human_strand"] = (
                closest_upstream_region["human_strand"]
            )

            result["upstream_region_source"] = (
                closest_upstream_region["source"]
            )

            result["upstream_region_distance"] = (
                closest_upstream_distance
            )

            result["upstream_region_mouse_chromosome"] = (
                closest_upstream_region["species_chromosome"]
            )

            result["upstream_region_mouse_start"] = (
                closest_upstream_region["species_start"]
            )

            result["upstream_region_mouse_end"] = (
                closest_upstream_region["species_end"]
            )

            result["upstream_region_mouse_strand"] = (
                closest_upstream_region["species_strand"]
            )

       
        # ADD THE DOWNSTREAM REGION INFORMATION
       

        if closest_downstream_region is not None:

            result["downstream_region_start"] = (
                closest_downstream_region["human_start"]
            )

            result["downstream_region_end"] = (
                closest_downstream_region["human_end"]
            )

            result["downstream_region_human_strand"] = (
                closest_downstream_region["human_strand"]
            )

            result["downstream_region_source"] = (
                closest_downstream_region["source"]
            )

            result["downstream_region_distance"] = (
                closest_downstream_distance
            )

            result["downstream_region_mouse_chromosome"] = (
                closest_downstream_region["species_chromosome"]
            )

            result["downstream_region_mouse_start"] = (
                closest_downstream_region["species_start"]
            )

            result["downstream_region_mouse_end"] = (
                closest_downstream_region["species_end"]
            )

            result["downstream_region_mouse_strand"] = (
                closest_downstream_region["species_strand"]
            )

       #I am deleting this as it was too early... '''lncRNA_neighbour_results.append(result)''' 
        #nua = result[] such that nua is the cordinate
        #nda = result["hdjfhjd"] nda in cordinate
        # need to get chromo, start, end, strand
        # get all lncRNAs in mouse on the corresponf chrom
        # get the lncRNA located between nua and nda region
        # get lncRNa that is closest to the nearest align region, and check orientation is conserved
        # if statements to check orientation 

        ####WE ARE NOW LOOKING FOR A MOUSE LNCRNA BETWEEN THE MAPPED NEAREST UPSTREAM ALIGNMENT (NUA) AND THE NEAREST DOWNSTREAM ALIGNMENT (NDA)

        conserved_mouse_lncRNA = []

        #here we are getting the coordinates of the nearest upstream alignment (NUA)
        nua_chr = result["upstream_region_mouse_chromosome"]
        nua_start = result["upstream_region_mouse_start"]
        nua_end = result["upstream_region_mouse_end"]

        #mouse coordinates of nearest downstream alignment (NDA)
        nda_chr = result["downstream_region_mouse_chromosome"]
        nda_start = result["downstream_region_mouse_start"]
        nda_end = result["downstream_region_mouse_end"]

        #the mouse interval must have the alignments on the same chromosome on mouse
        if (
            nua_chr != ""
            and nda_chr != ""                                      #add in somewhere around lines 750 to 840 that we need to account for the lncRNA being inverted (the lncRNA but either converge or diverge from the nda or the nua in the same way relative to the human nda or nua and its lncRNA)
            and nua_chr == nda_chr
        ):

            #now we need to see which mapped alignent is closest/comes first in the mouse (ensuring it still works if the human region is inverted in mouse)
            if nua_start < nda_start:

                mouse_interval_start = nua_end
                mouse_interval_end = nda_start

            else:
                mouse_interval_start = nda_end
                mouse_interval_end = nua_start

            #we have the lncRNAs in blocks of isoforms so we must find all of the mouse lncRNAs annotated on this chromosome
            mouse_lncRNAs_on_chr = mouseLncRNAs.get(nua_chr, [])

            #looking for a lncRNA mapped within this interval
            for mouse_lncRNA in mouse_lncRNAs_on_chr:

                if (
                    mouse_lncRNA["start"] <= mouse_interval_end
                    and mouse_lncRNA["end"] >= mouse_interval_start
                ):
                    conserved_mouse_lncRNA.append(mouse_lncRNA)

            '''#I want to save teh possible mouse lncRNA candidates
            if len(conserved_mouse_lncRNA) !=0:

                result["mouse_lncRNA_candidates"] = ";".join(
                    (
                        f'{candidate["name"]}|'
                        f'{candidate["gene_id"]}|'
                        f'{nua_chr}:'
                        f'{candidate["start"]}-'
                        f'{candidate["end"]}|'
                        f'{candidate["strand"]}|'
                    )
                    for candidate in conserved_mouse_lncRNA
                )''' #this is not neccessary for the final script but interesting to see how it works 
                
        #I need to save the result 
        if len(conserved_mouse_lncRNA) !=0:
            #select the closest lncRNA 
            #------            hlnc1.---------
            #------- mlnc1.     mlnc2----------    #I would say the the mlnc2 is the closest and conserved one to the hlnc1 as it is in the same position relative to the nda in the human 

            #now I am selecting the closest mouse lncRNA
            #first i must initialise variables used to select the best mouse lncRNA
            selected_mouse_lncRNA = None
            selected_mouse_distance = None
            selected_orientation_conserved = False
        
            # FIND WHICH CONSERVED FEATURE IS CLOSEST TO THE HUMAN lncRNA

            anchor_distances = {}

            #this condition means thar you only let teh PCG compete for being teh conserevd anchor is the blast mapped PCG for mouse is on the same chromosome as the WGA interval
            #if teh PCG is on a different chromosome to teh GWA the script ignores it for anchro selectiona and uses teh NUA or the NDA
            if (
                closest_pcg_distance is not None
                and result["closest_PCG_mouse_chromosome"] != ""
                and result["closest_PCG_mouse_chromosome"] == nua_chr
                and result["closest_PCG_mouse_chromosome"] == nda_chr
            ):
                anchor_distances["PCG"] = closest_pcg_distance

            if closest_upstream_distance is not None:
                anchor_distances["NUA"] = closest_upstream_distance

            if closest_downstream_distance is not None:
                anchor_distances["NDA"] = closest_downstream_distance


            # Get the closest human conserved feature
            closest_anchor = min(
                anchor_distances,
                key=anchor_distances.get
            )
     
            # FIND THE MOUSE lncRNA CLOSEST TO THAT SAME CONSERVED FEATURE
            
            for mouse_lncRNA in conserved_mouse_lncRNA:

                # --------------------------------------------
                # CASE 1: closest conserved feature is the PCG
                # --------------------------------------------
                if closest_anchor == "PCG":

                    anchor_start = result["closest_PCG_mouse_start"]
                    anchor_end = result["closest_PCG_mouse_end"]

                    mouse_distance = distance_between_regions(
                        mouse_lncRNA["start"],
                        mouse_lncRNA["end"],
                        anchor_start,
                        anchor_end
                    )

                    # Human lncRNA relative to human PCG
                    human_same_orientation = (
                        lncRNA["strand"]
                        == result["closest_PCG_strand"]
                    )

                    # Mouse lncRNA relative to mouse PCG
                    mouse_same_orientation = (
                        mouse_lncRNA["strand"]
                        == result["closest_PCG_mouse_strand"]
                    )

                    orientation_conserved = (
                    human_same_orientation
                    == mouse_same_orientation
                    )

                    # If human lncRNA overlaps the PCG,
                    # mouse lncRNA must also overlap its mapped PCG.
                    if closest_pcg_distance == 0 and mouse_distance != 0:
                        orientation_conserved = False

               
                # CASE 2: closest conserved feature is the NUA
                elif closest_anchor == "NUA":

                    anchor_start = nua_start
                    anchor_end = nua_end

                    mouse_distance = distance_between_regions(
                        mouse_lncRNA["start"],
                        mouse_lncRNA["end"],
                        anchor_start,
                        anchor_end
                    )
                                                                                    
                    # Is the human lncRNA transcribing toward or
                        # away from its human NUA?
                    human_relationship = transcription_relative_to_anchor(
                        lncRNA["start"],
                        lncRNA["end"],
                        lncRNA["strand"],
                        result["upstream_region_start"],                            result["upstream_region_end"]
                    )

                    # Is this mouse lncRNA transcribing toward or
                    # away from the corresponding mouse NUA?
                    mouse_relationship = transcription_relative_to_anchor(
                        mouse_lncRNA["start"],
                        mouse_lncRNA["end"],
                        mouse_lncRNA["strand"],
                        nua_start,                            
                        nua_end
                    )

                    orientation_conserved = (
                        human_relationship
                        == mouse_relationship
                    )


                # --------------------------------------------
                # CASE 3: closest conserved feature is the NDA
                # --------------------------------------------
                else:

                    anchor_start = nda_start
                    anchor_end = nda_end

                    mouse_distance = distance_between_regions(
                        mouse_lncRNA["start"],
                        mouse_lncRNA["end"],
                        anchor_start,
                        anchor_end
                    )

                # Human lncRNA relative to human NDA
                    human_relationship = transcription_relative_to_anchor(
                        lncRNA["start"],
                        lncRNA["end"],
                        lncRNA["strand"],
                        result["downstream_region_start"],
                        result["downstream_region_end"]
                    )

                    # Mouse lncRNA relative to mapped mouse NDA
                    mouse_relationship = transcription_relative_to_anchor(
                        mouse_lncRNA["start"],
                        mouse_lncRNA["end"],
                        mouse_lncRNA["strand"],
                        nda_start,
                        nda_end
                    )

                    orientation_conserved = (
                        human_relationship
                        == mouse_relationship
                    )

                # SELECT THE BEST MOUSE lncRNA
            

                # first candidate
                if selected_mouse_distance is None:
                    selected_mouse_distance = mouse_distance
                    selected_mouse_lncRNA = mouse_lncRNA
                    selected_orientation_conserved = orientation_conserved

                # orientation-conserved candidate beats a
                # non-orientation-conserved candidate
                elif orientation_conserved and not selected_orientation_conserved:
                    selected_mouse_distance = mouse_distance
                    selected_mouse_lncRNA = mouse_lncRNA
                    selected_orientation_conserved = True

                # if both have the same orientation status,
                # choose whichever is closer
                elif (
                    orientation_conserved == selected_orientation_conserved
                    and mouse_distance < selected_mouse_distance
                ):
                    selected_mouse_distance = mouse_distance
                    selected_mouse_lncRNA = mouse_lncRNA
                    selected_orientation_conserved = orientation_conserved


            if selected_mouse_lncRNA is not None:
            

                # always save the best selected mouse lncRNA candidate in the main result table
                result["selected_mouse_lncRNA_name"] = (
                    selected_mouse_lncRNA["name"]
                )

                result["selected_mouse_lncRNA_gene_id"] = (
                    selected_mouse_lncRNA["gene_id"]
                )

                result["selected_mouse_lncRNA_chromosome"] = nua_chr

                result["selected_mouse_lncRNA_start"] = (
                    selected_mouse_lncRNA["start"]
                )

                result["selected_mouse_lncRNA_end"] = (
                    selected_mouse_lncRNA["end"]
                )

                result["selected_mouse_lncRNA_strand"] = (
                    selected_mouse_lncRNA["strand"]
                )

                result["selected_mouse_lncRNA_anchor"] = (
                    closest_anchor
                )

                result["selected_mouse_lncRNA_distance"] = (
                    selected_mouse_distance
                )

                result["selected_mouse_lncRNA_orientation_conserved"] = (
                    selected_orientation_conserved
                )

                #only call it conserved if orientation is conserved
                if selected_orientation_conserved:

                    conservedlncrnas[lncRNA["gene_id"]] = {
                        "human_name": lncRNA["name"],
                        "mouse_name": selected_mouse_lncRNA["name"],
                        "mouse_gene_id": selected_mouse_lncRNA["gene_id"],
                        "closest_anchor": closest_anchor,
                        "distance": selected_mouse_distance,
                    }
                
        #save teh completed result for this human lncRNA
        lncRNA_neighbour_results.append(result)

#print the conservedlncrnas in a format that is useful
###to complete the script we now need to write a single row for every human lncRNA
output_file = "Conserved_Mouse_lncRNAs.tsv"

with open(output_file, "w", newline="") as outfile:

    writer = csv.writer(outfile, delimiter="\t")

    writer.writerow([
        "human_gene_id",
        "human_name",
        "mouse_gene_id",
        "mouse_name",
        "closest_anchor",
        "distance"
    ])

    for gene_id in sorted(conservedlncrnas):

        record = conservedlncrnas[gene_id]

        writer.writerow([
            gene_id,
            record["human_name"],
            record["mouse_gene_id"],
            record["mouse_name"],
            record["closest_anchor"],
            record["distance"]
        ])

print("Saved conserved lncRNAs to:", output_file)

output_file = "Human_Mouse_lncRNA_Synteny_Results.tsv"

fieldnames = [
    "human_chromosome",
    "lncRNA_name",
    "lncRNA_gene_id",
    "lncRNA_start",
    "lncRNA_end",
    "lncRNA_strand",

    "closest_PCG_human_name",
    "closest_PCG_start",
    "closest_PCG_end",
    "closest_PCG_strand",
    "closest_PCG_position",
    "closest_PCG_distance",

    "closest_PCG_mouse_name",
    "closest_PCG_mouse_chromosome",
    "closest_PCG_mouse_start",
    "closest_PCG_mouse_end",
    "closest_PCG_mouse_strand",

    "upstream_region_start",
    "upstream_region_end",
    "upstream_region_human_strand",
    "upstream_region_source",
    "upstream_region_distance",
    "upstream_region_mouse_chromosome",
    "upstream_region_mouse_start",
    "upstream_region_mouse_end",
    "upstream_region_mouse_strand",

    "downstream_region_start",
    "downstream_region_end",
    "downstream_region_human_strand",
    "downstream_region_source",
    "downstream_region_distance",
    "downstream_region_mouse_chromosome",
    "downstream_region_mouse_start",
    "downstream_region_mouse_end",
    "downstream_region_mouse_strand",

    "selected_mouse_lncRNA_name",
    "selected_mouse_lncRNA_gene_id",
    "selected_mouse_lncRNA_chromosome",
    "selected_mouse_lncRNA_start",
    "selected_mouse_lncRNA_end",
    "selected_mouse_lncRNA_strand",
    "selected_mouse_lncRNA_anchor",
    "selected_mouse_lncRNA_distance",
    "selected_mouse_lncRNA_orientation_conserved"
    ]

with open(output_file, "w", newline="") as outfile:

    writer = csv.DictWriter(
        outfile,
        fieldnames=fieldnames,
        delimiter="\t"
    )

    writer.writeheader()

    for result in lncRNA_neighbour_results:
        writer.writerow(result)

print("\nNeighbour results saved to:", output_file)
print(
    "Number of human lncRNAs written:",
    len(lncRNA_neighbour_results)
)