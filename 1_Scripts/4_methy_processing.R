# This script is meant as the methylation values pre-processing guideline

#--------------------Load libraries-----------------------
library(TCGAbiolinks)            # Version: 2.36.0
library(SummarizedExperiment)    # Version: 1.48.1
library(tidyverse)               # Version: 2.0.0

#--------------------Load object--------------------------
# To understand where this object came from, check the 1_get_data.R script
samples_data <- read.table("3_Data/samples_data.tsv", header = T, sep='\t')

#--------------------Prepare data-------------------------
# Get miRNA expression values
met_query <-  GDCquery(project = "TCGA-LIHC",                        # Liver hepatocellular carcinoma project
                       data.category = "DNA Methylation",            # DNA methylation data
                       platform = "Illumina Human Methylation 450",  # CpG detection platform
                       data.type = "Methylation Beta Value",         # Data type
                       barcode = samples_data$Barcode)               # Barcode

# The resulting object already has the methylation beta values by sample and probe as row names
met_data <- GDCprepare(met_query,                                    # Query object
                       directory = "3_Data/GDCdata/",                # Directory where files are stored
                       summarizedExperiment = T)                     # Create summarized experiment (includes little metadata)

# Filter data with just 0's asidde from NAs
met_data <- met_data[rowSums(met_data@assays@data@listData[[1]],
                             na.rm=T) != 0, ,drop = FALSE]

# Filter data that is missing in 20% or more of the samples
met_data <- met_data[which(!rowSums(
  is.na(met_data@assays@data@listData[[1]])) >
    (ncol(met_data@assays@data@listData[[1]])*0.2)),]

#--------------------Differential methylation-------------
TCGAanalyze_DMC(met_data,
                groupCol = "definition",
                p.cut = 0.5,
                title = "Differential methylation",
                save.directory = "3_Data/")
