# This script is meant for downloading all the necessary data fro the project. Keep in mind that, while this script
# is designed to work with the current (03-07-2025) data available in the TCGA platform, it will most likely be
# outdated in some aspects in the future. Please keep this in mind and modify it if necessary.

#--------------------Load libraries-----------------------
library(TCGAbiolinks)            # Version: 2.36.0
library(SummarizedExperiment)    # Version: 1.48.1
library(tidyverse)               # Version: 2.0.0

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
# The current TCGA cases' ids are composed with various characters (Eg."TCGA-G3-A3CJ-01A-11R-A213-07"). 
# The way to make sure the same samples are selected despite their data type (RNA or CpG) is cropping 
# the id and selecting the common part. Right now, the first 19 characters are enough to provide this information.
# Additionally, the work will only focus on primary tumors and any available controls, so filtering is required.
## Check for sample types
table(exp_res$sample_type) # (371 primary, 3 recurrent, and 50 normal) (08-07-2025)

## Keep only primary and normal
exp_res <- exp_res[which(exp_res$sample_type != "Recurrent Tumor"),]

## Extract initial samples' names
cases <- exp_res$cases

## Crop the samples' ids
cases <- substr(cases, 1, 19)

## Check how many are of the cases are shared among the three data types
Barcode <- cases[cases %in% substr(met_res$cases, 1, 19)] %>% 
  .[. %in% substr(mir_res$cases, 1 ,19)]

## 407 cases are shared among all three data types
length(Barcode)                  

#--------------------Data download------------------------
# This section is in charge of downloading the queries to the TCGA database. The barcode tells the TCGA
# platform which specific samples I want to request.
## mRNA
exp_query <- GDCquery(project = "TCGA-LIHC",                         # Liver hepatocellular carcinoma project
                      data.category = "Transcriptome Profiling",     # Refers to RNA
                      data.type = "Gene Expression Quantification",  # Gene expression 
                      workflow.type = "STAR - Counts",               # How reads are set as counts per gene
                      barcode = Barcode)                             # Barcode
GDCdownload(exp_query, "3_Data/GDCdata")                             # Downloading files

## miRNA
mir_query <- GDCquery(project = "TCGA-LIHC",                         # Liver hepatocellular carcinoma project
                      data.category = "Transcriptome Profiling",     # Refers to RNA
                      data.type = "miRNA Expression Quantification", # miRNA gen expression
                      barcode = Barcode)                             # Barcode
GDCdownload(mir_query, "3_Data/GDCdata")                             # Downloading files

## Methylation data
met_query <- GDCquery(project = "TCGA-LIHC",                         # Liver hepatocellular carcinoma project
                      data.category = "DNA Methylation",             # DNA methylation data
                      platform = "Illumina Human Methylation 450",   # CpG detection platform
                      data.type = "Methylation Beta Value",          # Data type
                      barcode = Barcode)                             # Barcode
GDCdownload(met_query, "3_Data/GDCdata", files.per.chunk = 50)       # Downloading files. I do not know why
                                                                     # I had to use the files.per.chunk argument

#--------------------Clinical data------------------------
# Available clinical data is also important as it can provide more insight into each sample. Not every
# characteristic will be used, but be sure to check which columns you wish to keep for further analysis.
## Get clinical data
cli_data <- GDCquery_clinic(project = "TCGA-LIHC",                   # Liver hepatocellular carcinoma project 
                            type = "clinical")                       # Acquire clinical data

## Select clinical features
cli_data <- cli_data %>% select(c("bcr_patient_barcode",
                                  "synchronous_malignancy",
                                  "prior_malignancy",
                                  "prior_treatment",
                                  "tumor_grade",
                                  "race",
                                  "gender",
                                  "vital_status"))

## Additional clinical data
subtypes <- TCGAquery_subtype("lihc")

## Select specific data
subtypes <- subtypes %>% select(c("Barcode",
                                  "ObesityClass1",
                                  "Alcoholic liver disease",
                                  "Hepatitis C",
                                  "Hepatitis B",
                                  "NAFLD",
                                  "Ploidy"))

## Formatting data
subtypes[which(subtypes$ObesityClass1 == "---"),]$ObesityClass1 <- "NA"
subtypes[which(subtypes$`Alcoholic liver disease` == "---"),]$`Alcoholic liver disease` <- "NA"
subtypes[which(subtypes$`Hepatitis C` == "---"),]$`Hepatitis C`<- "NA"
subtypes[which(subtypes$`Hepatitis B` == "---"),]$`Hepatitis B` <- "NA"
subtypes[which(subtypes$NAFLD == "---"),]$NAFLD <- "NA"
cli_data[which(cli_data$synchronous_malignancy == "Not Reported"),]$synchronous_malignancy <- "NA"
cli_data[which(cli_data$prior_malignancy == "not reported"),]$prior_malignancy <- "NA"
cli_data[which(cli_data$race == "not reported"),]$race <- "Unknown"

#--------------------Save objects-------------------------
## Check for missing data or duplicates
sum(substr(Barcode, 1, 12) %in%             # All samples are present in the clinical data
      cli_data$bcr_patient_barcode)

sum(substr(Barcode, 1, 16) %in%             # Only 189 samples have the subtype data
      subtypes$Barcode)

## Create data frame and combine data
sample_data <- data.frame(cbind(Barcode,    # Barcode, patient id, sample type
                                substr(Barcode, 1, 12),
                                exp_res[substr(exp_res$cases, 1, 19) %in% Barcode,]$sample_type))


