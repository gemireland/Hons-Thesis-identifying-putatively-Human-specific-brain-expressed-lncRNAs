#!/usr/bin/env python3

#this script is looking at all annotated human lncRNAs to see if they are expressed in BrainSpan's dataset
#we define expression as having an RPKM (as this is a gene-level set) > 1

#pandas is a tool that lets us read in a table and allows for easy editing and filtering of the data
import pandas as pd

#defining our files for the expression analysis 

human_lncRNA_file = "/scratch/irlgem001/Step4-RNASeq/HumanExpression/BrainSpan/genes_matrix_csv/Human_lncRNA_blocks.tsv"
brainspan_expression_file = "/scratch/irlgem001/Step4-RNASeq/HumanExpression/BrainSpan/genes_matrix_csv/expression_matrix.csv"
#the rows metadata file is a table that has the gene_id and gene_name for each row in the expression matrix (downloaded from same link as the actual expression dataset)
rows_metadata_file = "/scratch/irlgem001/Step4-RNASeq/HumanExpression/BrainSpan/genes_matrix_csv/rows_metadata.csv"

Humlnc_vs_BrainSpan_output = "HumanlncRNA_BrainSpan_expression_output.tsv"
final_humlncvsbrainspan_expression = "FinalHumanlncRNA_BrainSpan_expressed.tsv"

#we need to define our threshold RPKM as 1 (eliminates low level expression that is actually transcriptional noise)
threshold_RPKM = 1

#to read in human lncRNA blocks 
human_lncRNAs = pd.read_csv(human_lncRNA_file, sep="\t")

#the gencode gene_id in the human_lncRNA_blocks has a suffix of .2 that the brainspan gene_id doesn't have so we need to remove it
#this is to allow for comparison across datasets
human_lncRNAs["ensembl_gene_id"] = (human_lncRNAs["lncrna_gene_id"].astype(str).str.split(".").str[0])

#to see how many annotated lncRNA blocks 
print("human lnc blocks", len(human_lncRNAs))

#loading the BrainSpan gene metadata dataset
rows = pd.read_csv(rows_metadata_file)
print("no. Brainspan genes", len(rows))


#loading brainspan expression file dataset
expression = pd.read_csv(brainspan_expression_file)
#the shape tells us the number of rows and columns are found in the expression matrix
print("expression matrix shape", expression.shape)


#for each gene found in BrainSpan, calculate its expression 
#the first row in the expression dataset is the row number of brainspan so the rest of the columns have the RPKM numbers for each gene in the rows metadata file
#iloc means to use the index of the rows and columns to select stuff
expression_values = expression.iloc[:, 1:]

#defining the maximum expression value for each gene in the expression matrix
rows["max_RPKM"] = expression_values.max(axis=1)

#idnetifying the number of samples wehre the gene is expressed abover the 0.1 threshold
rows["num_samples_expressedabove1"] = (expression_values > threshold_RPKM).sum(axis=1)

#we define expression as having >0.1 in any sample (doesn't have to be every sample/tissue as lncRNAs are tissue-specific)
rows["BrainSpan_expressed"] = rows["max_RPKM"] > threshold_RPKM

#now looking for human lncRNAs that match expression in BrainSpan dataset
brainspan = rows[["ensembl_gene_id", "max_RPKM", "num_samples_expressedabove1", "BrainSpan_expressed"]]

#on left means we retain all rows from left dataframe (the human lncRNA blocks) and only those rows from the right datafram(the brainspan expression) with gene_id matching that in lncRNA blocks
results = human_lncRNAs.merge(brainspan, on="ensembl_gene_id", how="left")

#idnetify lncRNAs that are not expressed in BrainSpan
results["BrainSpan_present"] = results["max_RPKM"].notna()

#write the results to a csv file
#we want to keep the index as false so that we don't have an extra column in the output file
#function is to csv but we sep by tab as it is a tsv file
results.to_csv(Humlnc_vs_BrainSpan_output, sep="\t", index=False)

#saving only the lncRNAs expressed in BrainSpan to their own file (will be incorporated when comparing the 3 datasets)
expressed_lncRNAs = results[results["BrainSpan_expressed"] ==True].copy()

expressed_lncRNAs.to_csv(final_humlncvsbrainspan_expression, sep="\t", index=False)

#printing results to terminal for collection of results
print("lncs in BrainSpan", results["BrainSpan_present"].sum())

print("lncs not in BrainSpan", (~results["BrainSpan_present"]).sum())

print("lncs expressed in BrainSpan >1 RPKM", len(expressed_lncRNAs))

print("lncs present but not >1 RPKM", ((results["BrainSpan_present"]) & (results["BrainSpan_expressed"] == False)).sum())

print("\ncreating file to save output")
print(Humlnc_vs_BrainSpan_output)
print(final_humlncvsbrainspan_expression)
      
