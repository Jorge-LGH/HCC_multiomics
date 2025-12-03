# This script is meant as the methylation values pre-processing guideline

#--------------------Load libraries-----------------------
library(TCGAbiolinks)            # Version: 2.36.0
library(SummarizedExperiment)    # Version: 1.48.1
library(tidyverse)               # Version: 2.0.0
library(methyLImp2)              # Version: 1.2.0

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

# Starting object qualities
dim(met_data)                                                        # 485,577 probes across 407 samples
sum(is.na(assay(met_data)))                                          # 31,589,946 NA's (need to impute values after some filtering)
sum(colSums(assay(met_data), na.rm = T) == 0)                        # No empty columns (Should worry if it wasn't the case)
sum(rowSums(assay(met_data), na.rm = T) == 0)                        # There are 64,213 rows with only 0's if we drop the NA's

# Filter data with just 0's aside from NAs
met_data <- met_data[rowSums(met_data@assays@data@listData[[1]],     # 421,364 probes remain
                             na.rm=T) != 0, ,drop = FALSE]

# Filter probes with missing values in 20% or more of the samples
met_data <- met_data[which(!rowSums(                                 # 396,413 probes remain
  is.na(met_data@assays@data@listData[[1]])) >
    (ncol(met_data@assays@data@listData[[1]])*0.2)),]

# Check for probes with ambiguous chromosome mapping
seqnames(rowRanges(met_data)) %>% table()                            # No ambiguous mapping apparently (25-11-2025)

#--------------------QC-----------------------------------
# Check beta values distribution


#--------------------Methylation data imputation----------
# Impute missing values with regression based method. See: https://doi.org/10.1186/s12859-020-03592-5
met_data_imputed <- methyLImp2(met_data,
                               type = "450K",
                               groups = met_data@colData@listData[["definition"]],
                               BPPARAM = BiocParallel::SnowParam(workers = 30))

#--------------------Differential methylation-------------
TCGAanalyze_DMC(met_data,
                groupCol = "definition",
                p.cut = 0.5,
                title = "Differential methylation",
                save.directory = "3_Data/")
