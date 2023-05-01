---
title: "GxEprs"
author: "Dovini Jayasinghe and Hong Lee"
date: "01/05/2023" #to be edited
output: pdf_document
---


# GxEprs
The 'GxEprs' is an R package for **genomic prediction** that uses a sophisticated method that has been enhanced for its prediction accuracy. It performs Genome Wide Assosciation Studies (GWAS) and Genome Wide Envinonment Interaction Studies (GWEIS) using plink2 which is easily accessible from within R.


# Package installation
The current github version of **GxEprs** can be installed via:
```
library(devtools)
install_github("DoviniJ/GxEprs") 
```
The current CRAN version of **GxEprs** can be installed via:
```
install.packages("GxEprs")
```

# Quick start
We will illustrate the usage of **GxEprs** using a few example datasets downloaded and white labelled from ukbiobank. Follow the step wise process:

Step 1: Download and install plink2.exe in your working directory. Link: https://www.cog-genomics.org/plink/2.0/

Step 2: Download all the provided data files (in **DoviniJ/GxEprs/inst** directory) into the working directory:
1) DummyData.fam
2) DummyData.bim
3) DummyData.bed
4) Bcov_discovery.txt
5) Bphe_discovery.txt
6) Bcov_target.txt
7) Bphe_target.txt
8) Qcov_discovery.txt
9) Qphe_discovery.txt
10) Qcov_target.txt
11) Qphe_target.txt

Step 3: Run the following code (functions should be run in the given order) to obtain the risk scores of individuals in the target dataset:
```
DummyData <- "DummyData"
Bphe_discovery <- "Bphe_discovery.txt"
Bcov_discovery <- "Bcov_discovery.txt"
Bphe_target <- "Bphe_target.txt"
Bcov_target <- "Bcov_target.txt"
n_confounders = 14
thread = 20
GWAS_binary(DummyData, Bphe_discovery, Bcov_discovery, n_confounders, thread)
GWEIS_binary(DummyData, Bphe_discovery, Bcov_discovery, n_confounders, thread)
summary_binary(n_confounders)
PRS_binary(DummyData)
results_binary(n_confounders)
```

# Data preparation

## File formats
1) DummyData <- "DummyData" #binary file name (DummyData.fam, DummyData.bim, DummyData.bed)
2) Bphe_discovery <- "Bphe_discovery.txt" #File format: FID, IID, PHE (1 = controls, 2 = cases)
3) Bcov_discovery <- "Bcov_discovery.txt" #File format: FID, IID, scaled_COV, scaled_COV_sq, all the confounders as separate columns
4) Bphe_target <- "Bphe_target.txt" #File format: FID, IID, PHE (0 = controls, 1 = cases)
5) Bcov_target <- "Bcov_target.txt" #File format: FID, IID, COV, COV_sq, all the confounders as separate columns

## Tasks of each function
1) GWAS_binary: This performs GWAS and outputs the file B_trd.sum which contains GWAS summary statistics of all additive SNP effects
2) GWEIS_binary: This performs GWEIS and outputs the files B_add.sum and G_gxe.sum which contain GWEIS summary statistics of all additive and interaction SNP effects
3) PRS_binary: This computes polygenic risk scores for each individual in the target dataset and outputs the files B_trd.sscore, B_trd.sscore and B_trd.sscore  
4) results_binary: This outputs the file Individual_risk_values_without_confounders.txt or Individual_risk_values.txt (depending on the number of confounders given as argument)containing all the calculated individual risk scores using the relevant genomic prediction model


# Contact 
Please contact Hong Lee (hong.lee@unisa.edu.au) or Dovini Jayasinghe (dovinij@gmail.com) for queries.
