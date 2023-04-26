# GxEprs_Dummy
Step 1: Download and install plink2.exe in your working directory. Link: https://www.cog-genomics.org/plink/2.0/

Step 2: Download all the data files into the working directory:
1) DummyData.fam
2) DummyData.bim
3) DummyData.bed
4) Bcov_discovery.txt
5) Bphe_discovery.txt
6) Bcov_target.txt
7) Bphe_target.txt

Step 3: Download the file containing source code, i.e. GxEprs.R

Step 4: Open GxEprs.R in from working directory. It contains the following R functions (line 1 to line 125) and test functions (line 130 to line 151):
1) GWAS_binary 
2) GWEIS_binary
3) PRS_binary
4) results_binary

Step 5: See the output files generated in the same working directory

#File formats
1) DummyData <- "DummyData" #binary file name (DummyData.fam, DummyData.bim, DummyData.bed)
2) Bphe_discovery <- "Bphe_discovery.txt" #File format: FID, IID, PHE (1 = controls, 2 = cases)
3) Bcov_discovery <- "Bcov_discovery.txt" #File format: FID, IID, scaled_COV, scaled_COV_sq, all the confounders as separate columns
4) Bphe_target <- "Bphe_target.txt" #File format: FID, IID, PHE (0 = controls, 1 = cases)
5) Bcov_target <- "Bcov_target.txt" #File format: FID, IID, scaled_COV, scaled_COV_sq, all the confounders as separate columns

#Tasks of each function
1) GWAS_binary: This performs GWAS and outputs the file B_trd.sum which contains GWAS summary statistics of all additive SNP effects
2) GWEIS_binary: This performs GWEIS and outputs the files B_add.sum and G_gxe.sum which contain GWEIS summary statistics of all additive and interaction SNP effects
3) PRS_binary: This computes polygenic risk scores for each individual in the target dataset and outputs the files B_trd.sscore, B_trd.sscore and B_trd.sscore  
4) results_binary: This outputs the file Individual_risk_values_without_confounders.txt or Individual_risk_values.txt (depending on the number of confounders given as argument)containing all the calculated individual risk scores using the relevant genomic prediction model



