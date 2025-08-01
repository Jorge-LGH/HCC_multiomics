# This script is meant as the mRNA pre-processing guideline. Normalization, feature selection, and other
# processing techniques will be applied to ensure the mRNA expression values across every sample.

#--------------------Load libraries-----------------------
library(TCGAbiolinks)            # Version: 2.36.0
library(SummarizedExperiment)    # Version: 1.48.1
library(tidyverse)               # Version: 2.0.0
library(biomaRt)                 # Version: 2.64.0
library(NOISeq)                  # Version: 2.52.0
library(edgeR)                   # Version: 4.6.2

#--------------------Load object--------------------------
# To understand where this object came from, check the 1_get_data.R script
samples_data <- read.table("3_Data/samples_data.tsv", header = T, sep='\t')

#--------------------Prepare data-------------------------
# Get mRNA expression values
exp_query <- GDCquery(project = "TCGA-LIHC",                         # Liver hepatocellular carcinoma project
                      data.category = "Transcriptome Profiling",     # Refers to RNA
                      data.type = "Gene Expression Quantification",  # Gene expression 
                      workflow.type = "STAR - Counts",               # How reads are set as counts per gene
                      barcode = samples_data$Barcode)                # Barcode
exp_data <- GDCprepare(exp_query,                                    # Query object
                       directory = "3_Data/GDCdata/",                # Directory where files are stored
                       summarizedExperiment = F)                     # Do not create a summarized experiment

# Keep only protein coding RNAs
exp_data <- exp_data[which(exp_data$gene_type == "protein_coding"),] # 19,962 genes identified

# Remove rows with no transcript reads 
exp_data <- exp_data[rowSums(exp_data[,-c("gene_id",                 # 19,507 genes remain
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

# Remove rows with no transcript reads again
exp_data <- as.data.frame(exp_data,row.names = rownames(exp_data))
exp_data <- exp_data[rowSums(exp_data) != 0, , drop = FALSE]

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
ann_data$length <- abs(ann_data$end_position - ann_data$start_position)        # Add length

# Remove non protein coding genes from annotation
ann_data <- ann_data[which(ann_data$gene_biotype == "protein_coding"),]        # 19,397 genes remain

# Remove transcripts with no annotation data
exp_data <- exp_data[which(rownames(exp_data) %in% ann_data$ensembl_gene_id),] # Still 19,937

#--------------------Check for biases---------------------
# Create noiseq object for it to be compatible with selected workflow
noiseqData <- readData(exp_data, 
                       factors = samples_data,
                       gc = ann_data[, c("ensembl_gene_id", "percentage_gene_gc_content")],
                       biotype = ann_data[, c("ensembl_gene_id", "gene_biotype")],
                       length = ann_data[, c("ensembl_gene_id", "length")])

counts_data <- dat(noiseqData, type = "countsbio", factor = "sample_type")

# Expression values for cancer samples and controls
explo.plot(counts_data, plottype = "boxplot")

# Expression values in CPM for cancer samples and controls
explo.plot(counts_data, plottype = "barplot")

# Visualize low CPM values
cpm_hist <- ggplot(exp_data, aes(x = rowMeans(cpm(exp_data, log = T)))) + 
  geom_histogram(colour = "blue", fill = "lightblue") + xlab("CPM") + 
  ylab("Genes") + 
  theme_classic() + 
  geom_vline(aes(xintercept = 0), linetype = "dashed", colour = "red")
cpm_hist
sum(rowMeans(cpm(exp_data, log = T))>0)/nrow(exp_data)*100                 # ~64% genes have CPM>0

# Check for transcript composition bias
cd_data <- dat(noiseqData, type = "cd", norm = F)                          # Reference sample is: TCGA-ZS-A9CG-01A-11R-A37K-07
table(cd_data@dat$DiagnosticTest[, "Diagnostic Test"])                     # 374 Failed and 32 Passed (01-08-2025)
explo.plot(cd_data, samples = sample(1:ncol(exp_data),10))                 # The result clearly shows composition bias

# Check for GC bias
# The results show a little effect on expression values based on the GC content, however in neither case do the values
# have a fit of 50% or over. The fit for the cancer samples is 49.13% with a p value of 2.3e-08. The fit for the controls
# is of 38.55% and a p value of 2.1e-05
gc_content <- dat(noiseqData, type = "GCbias", k = 0, factor = "sample_type")
explo.plot(gc_content)  

# Check for length bias
# The results show a little effect on expression values based on the GC content. The fit for the cancer samples is 51.36% 
# with a p value of 1.5e-09. The fit for the controls is of 38.88% and a p value of 7.4e-06
len_bias <- dat(noiseqData, k = 0, type = "lengthbias", factor = "sample_type")
explo.plot(len_bias)

# Check for batch effect
# The samples do aggregate mostly by cancerous and control samples
myPCA = dat(noiseqData, type = "PCA", norm = F, logtransf = F)
explo.plot(myPCA, factor = "sample_type")