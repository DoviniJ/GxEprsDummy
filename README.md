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
install_github("DoviniJ/GxEprsDummy") 
```
# Load the library
```
library(GxEprsDummy)
```

# Quick start
We will illustrate the usage of **GxEprs** using a few example datasets downloaded and white labelled from ukbiobank. Follow the step wise process:

##### Step 1: Download and install plink2.exe in your machine. The package supports both Linux and Windows version.
Link: https://www.cog-genomics.org/plink/2.0/

##### Set 2: Obtain the path to the executable plink application <plink_path>

##### Step 3.1: Run the following code (functions should be run in the given order) to obtain the risk scores of individuals in the target dataset, when the outcome variable is binary (using regular model for binary outcome):
```
plink_path <- "<plink_path>/plink2" #give the path where plink executable file is located
inst_path <- system.file(package = "GxEprsDummy")
DummyData <- paste0(inst_path, "/DummyData")
Bphe_discovery <- paste0(inst_path, "/Bphe_discovery.txt")
Bcov_discovery <- paste0(inst_path, "/Bcov_discovery.txt")
Bphe_target <- paste0(inst_path, "/Bphe_target.txt")
Bcov_target <- paste0(inst_path, "/Bcov_target.txt")
n_confounders = 14
thread = 20
setwd("<path to working directory>") #set the working directory where you need to save the output files
GWAS_binary(plink_path, DummyData, Bphe_discovery, Bcov_discovery, n_confounders, thread)
GWEIS_binary(plink_path, DummyData, Bphe_discovery, Bcov_discovery, n_confounders, thread)
PRS_binary(plink_path, DummyData)
summary_regular_binary(Bphe_target, Bcov_target, n_confounders)
summary_permuted_binary(Bphe_target, Bcov_target, n_confounders)
results_regular_binary(Bphe_target, Bcov_target, n_confounders)
```
##### Step 3.2: Run the following code (functions should be run in the given order) to obtain the risk scores of individuals in the target dataset, when the outcome variable is binary (using permuted model for binary outcome):
```
results_permuted_binary(Bphe_target, Bcov_target, n_confounders)
```

Note: It is recommended to fit both regular and permuted models and obtain the summary of both fitted models (using ```summary_regular_binary(n_confounders)``` and ```summary_permuted_binary(n_confounders)```. If the 'PRS_gxe x E' term is significant of insignificant in both the models, any model could be used to obtain results (i.e. ```results_regular_binary(n_confounders)``` or ```results_permuted_binary(n_confounders)```). If the 'PRS_gxe x E' term is significant in one model, and insignificant in other model, it is advised to use the permuted model to obtain results (i.e. ```results_permuted_binary(n_confounders)```).





<!--- 


##### Step 3.3: Run the following code (functions should be run in the given order) to obtain the risk scores of individuals in the target dataset, when the outcome variable is quantitative (using regular model for quantitative outcome):
```
DummyData <- "DummyData"
Qphe_discovery <- "Qphe_discovery.txt"
Qcov_discovery <- "Qcov_discovery.txt"
Qphe_target <- "Qphe_target.txt"
Qcov_target <- "Qcov_target.txt"
n_confounders = 14
thread = 20
GWAS_quantitative(DummyData, Bphe_discovery, Bcov_discovery, n_confounders, thread)
GWEIS_quantitative(DummyData, Bphe_discovery, Bcov_discovery, n_confounders, thread)
PRS_quantitative(DummyData)
summary_regular_quantitative(Qphe_target, Qcov_target, n_confounders)
summary_permuted_quantitative(Qphe_target, Qcov_target, n_confounders)
results_regular_quantitative(Qphe_target, Qcov_target, n_confounders)
```
##### Step 3.4: Run the following code (functions should be run in the given order) to obtain the risk scores of individuals in the target dataset, when the outcome variable is quantitative (using permuted model for quantitative outcome):
```
results_permuted_quantitative(Qphe_target, Qcov_target, n_confounders)
```

Note: It is recommended to fit both regular and permuted models and obtain the summary of both fitted models (using ```summary_regular_quantitative(n_confounders)``` and ```summary_permuted_quantitative(n_confounders)```. If the 'PRS_gxe x E' term is significant of insignificant in both the models, any model could be used to obtain results (i.e. ```results_regular_quantitative(n_confounders)``` or ```results_permuted_quantitative(n_confounders)```). If the 'PRS_gxe x E' term is significant in one model, and insignificant in other model, it is advised to use the permuted model to obtain results (i.e. ```results_permuted_quantitative(n_confounders)```). 


-->





# Data preparation

## File formats
### Input files
1) DummyData.fam - This is one of the binary files which contains family ID (FID), individual ID (IID), father's ID, mother's ID, sex and phenotype value as columns. The example dataset has 10,000 individuals. Note that the file has no column headings.   
![image](https://user-images.githubusercontent.com/131835334/236634478-99a10ec5-2e05-4259-981e-d67562b1a06a.png)
  
2) DummyData.bim - This is one of the binary files which contains chromosome code, SNP ID, position of centimorgans, base-pair coordinate, minor allele and reference allele as columns. The example dataset has 10,000 SNPs. Note that the file has no column headings.
![image](https://user-images.githubusercontent.com/131835334/236634694-5dbe6a29-5ae0-44c9-b076-b80fcabb7144.png)

3) DummyData.bed - This is also a binary file which cannot be read by humans.
4) Bphe_discovery.txt - This is a .txt file which contains FID, IID and binary phenotype (1=controls, 2=cases) of the discovery sample as columns. The discovery dataset has 7916 individuals. Note that the file has no column headings.    
![image](https://user-images.githubusercontent.com/131835334/236635016-88560176-a22a-4863-b200-4ddca8ca6980.png)

5) Bcov_discovery.txt - This is a .txt file which contains FID, IID, standardized covariate, square of the standardized covariate and 14 confounders of the discovery sample as columns. The discovery dataset has 7916 individuals. Note that the file has no column headings.    
![image](https://user-images.githubusercontent.com/131835334/236635276-7e1c6d92-3a84-4f9e-b68b-a171d9684da3.png)

6) Bphe_target.txt - This is a .txt file which contains FID, IID and binary phenotype (0=controls, 1=cases) of the target sample as columns. The target dataset has 1939 individuals who are independent from the discovery dataset. Note that the file has no column headings.   
![image](https://user-images.githubusercontent.com/131835334/236635388-53c3ff05-ae8b-498c-8354-0e1419aaf56f.png)

7) Bcov_target.txt - This is a .txt file which contains FID, IID, standardized covariate, square of the standardized covariate and 14 confounders of the target sample as columns. The target dataset has 1939 individuals who are independent from the discovery dataset. Note that the file has no column headings.   
![image](https://user-images.githubusercontent.com/131835334/236635508-d08cec8f-38bb-4008-8bc5-18ad87c8eece.png)


### Output files
1) B_trd.sum - This contains GWAS summary statistics of all additive SNP effects, when the outcome is binary. 
2) B_add.sum - This contains GWEIS summary statistics of all additive SNP effects, when the outcome is binary.
3) B_gxe.sum - This contains GWEIS summary statistics of all interaction SNP effects, when the outcome is binary.
4) B_trd.sscore - This contains
5) B_add.sscore - This contains
6) B_gxe.sscore - This contains
7) Bsummary.txt - This contains the target model (either regular model or permuted model) summary output, when the outcome is binary.
8) Individual_risk_values.txt - This contains all the calculated individual risk scores using the target model (either regular model or permuted model), when the outcome is binary.

## Tasks of each function
1) GWAS_binary - This performs GWAS and outputs the file B_trd.sum which contains GWAS summary statistics of all additive SNP effects
2) GWEIS_binary - This performs GWEIS and outputs the files B_add.sum and B_gxe.sum which contain GWEIS summary statistics of all additive and interaction SNP effects
3) PRS_binary - This computes polygenic risk scores for each individual in the target dataset and outputs the files B_trd.sscore, B_trd.sscore and B_trd.sscore  
4) summary_regular_binary - This outputs the file Bsummary.txt which gives the summary of the fitted **regular** model for **binary** outcome
5) summary_permuted_binary - This outputs the file Bsummary.txt which gives the summary of the fitted **permuted** model for **binary** outcome
6) results_regular_binary - This outputs the file Individual_risk_values.txt containing all the calculated individual risk scores using the **regular** genomic prediction model for **binary** outcome
7) results_permuted_binary - This outputs the file Individual_risk_values.txt containing all the calculated individual risk scores using the **permuted** genomic prediction model for **binary** outcome


# Contact 
Please contact Hong Lee (hong.lee@unisa.edu.au) or Dovini Jayasinghe (dovinij@gmail.com) for queries.
