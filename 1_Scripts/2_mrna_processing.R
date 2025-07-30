# This script is meant as the mRNA pre-processing guideline. Normalization, feature selection, and other
# processing techniques will be applied to ensure the mRNA expresison values across every sample.

#--------------------Load libraries-----------------------
library(TCGAbiolinks)            # Version: 2.36.0
library(SummarizedExperiment)    # Version: 1.48.1
library(tidyverse)               # Version: 2.0.0
library(biomaRt)                 # Version: 2.64.0

#--------------------Load object--------------------------
# To understabd where this object came from, check the 1_get_data.R script
samples_data <- read.table("3_Data/samples_data.tsv", header = T, sep='\t')

#--------------------Prepare data-------------------------
# Get mRNA expression values
exp_query <- GDCquery(project = "TCGA-LIHC",                         # Liver hepatocellular carcinoma project
                      data.category = "Transcriptome Profiling",     # Refers to RNA
                      data.type = "Gene Expression Quantification",  # Gene expression 
                      workflow.type = "STAR - Counts",               # How reads are set as counts per gene
                      barcode = Barcode)                             # Barcode
exp_data <- GDCprepare(exp_query,                                    # Query object
                       directory = "3_Data/GDCdata/",                # Directory where files are stored
                       summarizedExperiment = F)                     # Do not create a summarized experiment

# Keep only protein coding RNAs
exp_data <- exp_data[which(exp_data$gene_type == "protein_coding"),] # 19,962 genes identified

# Remove rows with no transcript reads 
exp_data <- exp_data[rowSums(exp_data[,-c("gene_id",                 # 129,507 genes remain
                                          "gene_name", 
                                          "gene_type")]) != 0,]

# Set ensembl gene id as row names
rownames(exp_data) <- exp_data$gene_id

# Remove ensembl id column as well as the gene name and gene type
exp_data <- dplyr::select(exp_data, -c("gene_id", "gene_name", "gene_type"))

# Keep only unstranded reads
keep_cols <- colnames(exp_data)[sapply(colnames(exp_data), function(col){
  strsplit(col, "_")[[1]][1] == "unstranded"})]
exp_data <- exp_data %>% dplyr::select(keep_cols)

# Remove the "unstranded_" part of the columns' names
colnames(exp_data) <- unlist(strsplit(colnames(exp_data), "_"))[
  unlist(strsplit(colnames(exp_data), "_")) != "unstranded"]

# Remove the transcript version and only keep ensembl gene id
rownames(exp_data) <- sapply(strsplit(rownames(exp_data),".",fixed=T),
                             function(x) x[1])

#--------------------Annotation data----------------------
# Get annotation data
mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
ann_data <- getBM(attributes = c("ensembl_gene_id", 
                                 "percentage_gene_gc_content", 
                                 "gene_biotype",
                                 "start_position",
                                 "end_position",
                                 "hgnc_id",
                                 "hgnc_symbol"),
                  filters = "ensembl_gene_id", 
                  values=rownames(exp_data), 
                  mart=mart)
ann_data$length <- abs(ann_data$end_position - ann_data$start_position)

# Remove non protein coding genes from annotation
ann_data <- ann_data[which(ann_data$gene_biotype == "protein_coding"),]

# Remove transcripts with no annotation data
exp_data <- exp_data[which(rownames(exp_data) %in% ann_data$ensembl_gene_id),]

#--------------------Check for biases---------------------



