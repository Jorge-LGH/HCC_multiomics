# This script is meant for downloading all the necessary data fro the project. Keep in mind that, while this script
# is designed to work with the current (03-07-2025) data available in the TCGA platform, it will most likely be
# outdated in some aspects in the future. Please keep this in mind and modify it if necessary.

#--------------------Load libraries-----------------------
library(TCGAbiolinks)            # Version: 2.36.0
library(SummarizedExperiment)    # Version: 1.48.1

#--------------------Query preparation--------------------
# This section is in charge of making the queries to the TCGA database
## mRNA
exp_query <- GDCquery(project = "TCGA-LIHC",                         # Liver hepatocellular carcinoma project
                      data.category = "Transcriptome Profiling",     # Refers to RNA
                      data.type = "Gene Expression Quantification",  # Gene expression 
                      workflow.type = "STAR - Counts")               # How reads are set as counts per gene

## miRNA
mir_query <- GDCquery(project = "TCGA-LIHC",                         # Liver hepatocellular carcinoma project
                      data.category = "Transcriptome Profiling",     # Refers to RNA
                      data.type = "miRNA Expression Quantification") # miRNA gen expression 

## Methylation data
met_query <- GDCquery(project = "TCGA-LIHC",                         # Liver hepatocellular carcinoma project
                      data.category = "DNA Methylation",             # DNA methylation data
                      platform="Illumina Human Methylation 450")     # CpG detection platform

#--------------------Getting query results----------------
# These functions get the results for the query in a manner that can be later parsed and used to download data
exp_res <- getResults(exp_query)
mir_res <- getResults(mir_query)
met_res <- getResults(met_query)

#--------------------Data pre-selection-------------------
# Identify the cases that have available data for mRNA, CpG islands, and miRNA
cases <- exp_res$cases

# Visit the 

#--------------------Data download------------------------






