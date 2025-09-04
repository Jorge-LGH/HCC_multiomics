# Inferencia de agentes reguladores en el carcinoma hepatocelular mediante la integración de datos multiómicos

This repository will contain all the work pertaining the Hepatocellular Carcinoma (HCC) multiomic integration project.

## Pipeline Overview

**1. Data Acquisition:** Download mRNA, miRNA, and methylation values for samples of HCC available in the TCGA databases.

- [1_get_data](./1_Scripts/1_get_data.R): Download data of selected cases.

**2. Data pre-processing:** Filter data, reduce present biases, normalize data, perform differential expression analyses, and store pre-processed data for all three data types.

- [2_mrna_processing.R](./1_Scripts/2_mrna_processing.R): Pre-process mRNA data.
- [3_mirna_processing.R](./1_Scripts/3_mirna_processing.R): Pre-process miRNA data.
- [4_methy_processing.R](./1_Scripts/4_methy_processing.R): Pre-process CpG data.