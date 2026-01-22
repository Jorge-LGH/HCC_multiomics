# This script is intended as a concatenation for the selected omic data types. In this case, they are
# miRNA, mRNA, and CpG sites. This concatenation of data types is based on the SGCCA method.

#-----------------------Load libraries-----------------------

#-----------------------Load data----------------------------
# Load only the blocks which have been normalized, scaled, and centered
exp_data <- read.table("3_Data/eig_exp_comparable.tsv",sep=',',row.names=T)     # mRNA expression data
mir_data <- read.table("3_Data/eig_mir_comparable.tsv",sep=',',row.names=T)     # miRNA expression data
                                                                                # CpG methylation data

#-----------------------Concatenate--------------------------
# The data concatenation is performed as to have the samples as columns and ALL features as rows
