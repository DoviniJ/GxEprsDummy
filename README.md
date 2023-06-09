---
title: "GxEprs"
author: "Dovini Jayasinghe, Md Moksedul Momin and Hong Lee"
date: "01/05/2023" #to be edited
output: pdf_document
---

# GxEprs
The 'GxEprs' is an R package to detect and estimate GxE. It uses a novel PRS model that can enhance the prediction accuracy by utilising GxE effects. Firstly it performs Genome Wide Association Studies (GWAS)  and Genome Wide Environment Interaction Studies (GWEIS) using the discovery dataset (see functions ```GWAS_binary()```,```GWAS_quantitative()```, ```GWEIS_binary()```, ```GWEIS_quantitative()```). See the section $\color{red}{IMPORTANT}$ for the discovery models used. Secondly, it uses the GWAS and GWEIS summary statistics generated from the fucntions above to obtain polygenic risk scores (PRSs) (see functions ```PRS_binary()``` and ```PRS_quantitative()```) for the target sample. Finally it predicts the risk values of each individual in the target sample (see functions ```summary_regular_binary()``` and ```summary_regular_quantitative()```). Note that the users can fit 4 different models when the outcome is a quantitative trait, and 5 different models when the outcome is a binary disease trait. See the section $\color{red} {IMPORTANT}$ for the target models used. Finally, it is recommended to check the p-value from permutations using Model 4 (see function ```summary_permuted_quantitative()```), and Model 5 (see function```summary_permuted_binary()```), to make sure that the significance of GxE is not spurious due to model misspecification (see references).

# Data preparation

## File formats
### Input files
1) mydata.fam - This is a file associated with the PLINK binary format file which contains the following columns in order. The example dataset has 10,000 individuals. Note that the file has no column headings. This follows the PLINK .fam file format.
* family ID (FID) 
* individual ID (IID) 
* father's ID 
* mother's ID 
* sex 
* phenotype value

```
1001019 1001019 0 0 2 -9
1001022 1001022 0 0 2 -9
1001035 1001035 0 0 1 -9
1001054 1001054 0 0 2 -9
1001078 1001078 0 0 2 -9
```
  
2) mydata.bim - This is is a file associated with the PLINK binary format file which contains the following columns in order. The example dataset has 10,000 SNPs. Note that the file has no column headings. This follows the PLINK .bim file format.
* chromosome code 
* SNP ID 
* position of centimorgans 
* base-pair coordinate 
* minor allele  
* reference allele 

```
1	  snp_53131969	0	754182	A	G
1	  snp_52286139	0	761732	C	T
1	  snp_512562034	0	768448	A	G
1	  snp_54040617	0	779322	G	A
1	  snp_52980300	0	785989	T	C
```

3) mydata.bed - This is the PLINK binary format file which includes genotype information. This follows the PLINK .bed file format.
4) Bpd.txt - This is a .txt file which contains the following columns in order. The discovery dataset has 7916 individuals. Note that the file has no column headings.
* FID 
* IID  
* binary phenotype (1=controls, 2=cases) of the discovery sample

```
1050405 1050405 1
1036224 1036224 2
1168718 1168718 1
1033226 1033226 1
1056980 1056980 1
```

5) Bcd.txt - This is a .txt file which contains the following columns in order. The discovery dataset has 7916 individuals. Note that the file has no column headings.    
* FID 
* IID 
* standardized covariate 
* square of the standardized covariate (Note: This column is required to address model misspecification, if present. See reference for details.) 
* 14 confounders of the discovery sample

```
1050405 1050405 0.787403812314451 0.620004763647331 -3.04026 45 -12.048 2.17634 -0.940322 -0.446351 -5.45685 -2.53161 -2.13435 -1.95623 -3.82792 -0.380636 0 10
1036224 1036224 -0.119722138532781 0.0143333904548625 -5.1054 64 -14.5169 6.01889 -3.85803 3.62625 5.10717 -3.54574 0.393994 3.64275 4.42975 -2.26704 1 19
1168718 1168718 0.173372721351375 0.0300581005087816 -1.91044 59 -12.7462 5.95244 0.0738914 1.80523 4.76284 0.130369 -1.05615 0.316777 0.988783 -1.76502 1 7
1033226 1033226 -0.699321184695051 0.48905011936329 -1.83526 68 -10.3349 4.71264 -1.84521 -0.524855 -3.80275 0.837965 0.265233 2.10903 -0.210259 0.71504 0 20
1056980 1056980 3.69300366651739 13.6382760809109 -3.15649 69 -8.56737 4.78248 -1.49547 -7.49413 -5.39887 1.85316 4.07476 1.05351 0.825942 -2.09669 1 20
```

6) Bpt.txt - This is a .txt file which contains the following columns in order. The target dataset has 1939 individuals who are independent from the discovery dataset. Note that the file has no column headings.   
* FID 
* IID  
* binary phenotype (0=controls, 1=cases) of the target sample

```
1001035 1001035 0
1001105 1001105 0
1001129 1001129 0
1001210 1001210 0
1001323 1001323 0
```

7) Bct.txt - This is a .txt file which contains the following columns in order. The target dataset has 1939 individuals who are independent from the discovery dataset. Note that the file has no column headings.   
* FID 
* IID 
* standardized covariate 
* square of the standardized covariate  (Note: This column is required to address model misspecification, if present. See reference for details.)
* 14 confounders of the target sample

```
1001035 1001035 -0.420822976931972 0.177091977913887 -4.12263 57 -13.5185 5.40198 -4.81994 0.664494 -4.92217 -0.451329 3.14677 0.42704 0.821306 -2.77705 1 7
1001105 1001105 -0.0805583280660987 0.00648964422080519 2.92534 56 -13.6236 3.21643 -0.856048 0.750187 -2.01798 -0.350832 5.10141 2.1807 -6.04343 1.78928 1 19
1001129 1001129 -1.32752644870189 1.76232647200306 -3.09118 61 -9.94475 3.60562 -0.917639 0.905664 -5.09843 -1.16329 -1.88102 -1.24154 0.699574 2.2442 0 20
1001210 1001210 0.698007239555549 0.487214106471958 4.58829 49 -12.5471 4.09467 -2.58951 6.06898 12.9822 -0.704179 2.90357 -0.334968 5.04274 0.66175 0 10
1001323 1001323 -0.657981606980219 0.432939795124272 -3.53948 56 -12.795 2.91524 -2.72794 3.61555 3.92957 -2.93899 -0.454737 2.31013 2.51783 -4.15592 0 7
```

8) Qpd.txt - This is a .txt file which contains the following columns in order. The discovery dataset has 6426 individuals. Note that the file has no column headings.
* FID 
* IID  
* quantitative phenotype of the discovery sample

```
1109732 1109732 31.6534
1158310 1158310 25.5035
1156676 1156676 26.7391
1096610 1096610 25.5271
1006946 1006946 26.7165
```

9) Qcd.txt - This is a .txt file which contains the following columns in order. The discovery dataset has 6426 individuals. Note that the file has no column headings.    
* FID 
* IID 
* standardized covariate 
* square of the standardized covariate  (Note: Although squared covariate term is not involved when analyzing quantitative traits, this file can/may be used as the covariate file when modulating with a binary outcome. In that event, this column is required to address model misspecification, if present. See reference for details.)
* 14 confounders of the discovery sample

```
1109732 1109732 -0.644020461494502 0.414762354823591 -3.83142 64 -14.0364 5.51742 0.0714337 5.66263 0.865562 -2.26957 -0.0965859 -2.35497 1.05889 0.195302 0 7
1158310 1158310 -0.0278698075809014 0.00077672617459647 0.614044 66 -10.8505 2.11998 -0.882883 -0.441662 -2.64177 2.78944 0.524586 2.67134 -2.63724 -0.998764 1 20
1156676 1156676 2.1286574811167 4.53118267191409 -0.237792 55 -9.75369 3.18343 -2.09793 6.87345 11.3777 2.96961 -1.11879 0.873649 3.35523 -4.57831 1 10
1096610 1096610 2.1286574811167 4.53118267191409 6.69866 47 -9.07045 0.956878 -2.48407 1.06359 -3.13247 2.1232 -0.00976751 0.820582 0.0305345 1.6303 1 20
1006946 1006946 -0.952095788451302 0.906486390386706 -1.61423 59 -12.9379 1.29461 -1.79973 1.44404 -6.82898 -2.96795 -2.91577 -1.82881 7.15892 2.10916 1 20
```

10) Qpt.txt - This is a .txt file which contains the following columns in order. The target dataset has 1579 individuals who are independent from the discovery dataset. Note that the file has no column headings.   
* FID 
* IID  
* quantitative phenotype of the target sample

```
1001054 1001054 26.5723
1001090 1001090 20.2632
1001203 1001203 27.7365
1001450 1001450 18.75
1001516 1001516 23.3025
```

11) Qct.txt - This is a .txt file which contains the following columns in order. The target dataset has 1579 individuals who are independent from the discovery dataset. Note that the file has no column headings.   
* FID 
* IID 
* standardized covariate 
* square of the standardized covariate (Note: Although squared covariate term is not involved when analyzing quantitative traits, this file can/may be used as the covariate file when modulating with a binary outcome. In that event, this column is required to address model misspecification, if present. See reference for details.) 
* 14 confounders of the target sample

```
1001054 1001054 -0.644020461494502 0.414762354823591 -3.82659 69 -13.8514 3.9608 -1.78805 0.0692473 -6.32556 2.85359 1.08516 -1.30304 3.41659 1.41577 0 7
1001090 1001090 -0.952095788451302 0.906486390386706 2.06515 60 -12.2438 4.04169 -0.905739 5.9656 8.35545 -1.43576 -0.618153 0.746918 5.11019 -0.207188 1 19
1001203 1001203 -0.0278698075809014 0.00077672617459647 -0.795863 62 -10.9195 6.91985 -2.92088 1.26019 -5.56624 -0.552624 -0.0756095 -0.910047 -1.33896 1.72636 0 7
1001450 1001450 -0.644020461494502 0.414762354823591 -2.62088 67 -9.9271 4.1096 -2.35454 0.719021 -1.82806 -1.82107 1.21574 -3.56693 -7.91232 2.71011 0 10
1001516 1001516 0.280205519375899 0.0785151330887172 -3.33164 67 -11.8637 5.88272 1.07288 2.74488 -7.32776 -2.39477 -3.07983 -1.43625 2.08822 1.42939 1 15
```

# Package installation
The current GitHub version of **GxEprs** can be installed via:
```
library(devtools)
install_github("DoviniJ/GxEprsDummy") 
```
# Load the library
```
library(GxEprsDummy)
```

# Quick start
We will illustrate the usage of **GxEprs** using a few example datasets downloaded and white labelled from UK biobank. Follow the step wise process:

##### Step 1: Download and install plink2 in your machine. The package supports both Linux and Windows version.
Link: https://www.cog-genomics.org/plink/2.0/

##### Step 2: Obtain the path to the executable plink application <plink_path>

##### Step 3.1: Run the following code (functions should be run in the given order) to obtain the risk scores of individuals in the target dataset:

###### Step 3.1.1 Give the path where plink executable file is located
```
plink_path <- "<plink_path>/plink2" 
```
###### Step 3.1.2 It is always recommended to check how the files look like before using them in functions, for better understanding. You may directly use the data files embedded in the package as a trial. Note that, for convenience, we have used identical names for the embedded data object, and for the corresponding function argument. You can check the top proportion of each data file using the following code:
```
head(Bphe_discovery) #phenotype file of the discovery sample when the outcome is binary
head(Bcov_discovery) #covariate file of the discovery sample when the outcome is binary
head(Bphe_target) #phenotype file of the target sample when the outcome is binary
head(Bcov_target) #covariate file of the target sample when the outcome is binary
head(Qphe_discovery) #phenotype file of the discovery sample when the outcome is quantitative
head(Qcov_discovery) #covariate file of the discovery sample when the outcome is quantitative
head(Qphe_target) #phenotype file of the target sample when the outcome is quantitative
head(Qcov_target) #covariate file of the target sample when the outcome is quantitative
```

###### Step 3.1.3 To call the data files saved in "inst" directory, you can follow the following code to obtain the path of each data file. 
```
inst_path <- system.file(package = "GxEprsDummy") 
DummyData <- paste0(inst_path, "/DummyData") #this contains all .fam, .bed and .bim files. They can be accessed by a direct call of prefix "DummyData"
Bphe_discovery <- paste0(inst_path, "/Bphe_discovery.txt")
Bcov_discovery <- paste0(inst_path, "/Bcov_discovery.txt")
Bphe_target <- paste0(inst_path, "/Bphe_target.txt")
Bcov_target <- paste0(inst_path, "/Bcov_target.txt")
Qphe_discovery <- paste0(inst_path, "/Qphe_discovery.txt")
Qcov_discovery <- paste0(inst_path, "/Qcov_discovery.txt")
Qphe_target <- paste0(inst_path, "/Qphe_target.txt")
Qcov_target <- paste0(inst_path, "/Qcov_target.txt")
```
Note that the step 3.1.3 described above is to call the embedded data files in this package itself. However, when users have to call their own data, they can follow the same approach. It is more convenient if the users can store all their data files in the same working directory. For example, assume that the file names are as follows 
(Refer to 'File formats' section of this document to view the formatting details of each of the following input file): 
* binary files: **mydata.fam**, **mydata.bim** and **mydata.bed**
* phenotype file of discovery sample (binary outcome): **Bpd.txt**
* covariate file of discovery sample (binary outcome): **Bcd.txt**
* phenotype file of target sample (binary outcome): **Bpt.txt**
* covariate file of the target sample (binary outcome): **Bct.txt**
* phenotype file of discovery sample (quantitative outcome): **Qpd.txt** 
* covariate file of discovery sample (quantitative outcome): **Qcd.txt** 
* phenotype file of target sample (quantitative outcome): **Qpt.txt**
* covariate file of the target sample (quantitative outcome): **Qct.txt**

###### Additional note:
_Note that, all these files can be placed in a separate location. It is always upto the users choice. In that case remember to give the full path to the file location since R identifies files by name, only when they are in the same directory._
```
b_file <- "<path>/mydata"
Bphe_discovery <- "<path>/Bpd.txt"
Bcov_discovery <- "<path>/Bcd.txt"
Bphe_target <- "<path>/Bpt.txt"
Bcov_target <- "<path>/Bct.txt"
Qphe_discovery <- "<path>/Qpd.txt"
Qcov_discovery <- "<path>/Qcd.txt"
Qphe_target <- "<path>/Qpt.txt"
Qcov_target <- "<path>/Qct.txt"
```

###### Step 3.1.4 Set the working directory and run the following R functions in the given order
```
setwd("<path to working directory>") #set the working directory where you need to save the output files
```
$\color{red}{NOTE:}$ Read **manual.pdf** document for descriptions of arguments passed for each function.

###### When the outcome variable is binary
**Command**
```
GWAS_binary(plink_path, "mydata", "Bpd.txt", "Bcd.txt", thread = 20, summary_output = "B_out")
```
As explained above, “mydata” is the prefix of the PLINK format files, “Bpd.txt” is binary phenotype file of the discovery sample, "Bcd.txt" is the covariate file of the discovery sample, thread indicates the number of CPUs used to run the command which can be optionally specified by the user (default is 20). Here summary_output is an output file name (also an optional argument) which can be defined by the user (default is "B_out" and the function will generate the suffix ".trd.sum") This command performs GWAS using a logistic regression, and outputs GWAS summary statistics of all additive SNP effects stored in the file "B_out.trd.sum".

**Output**
```
V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12 V13 V14
1 754182 snp_53131969 G A A N ADD 7867 0.10380130071374 0.112336 0.923983 0.355495 .
1 761732 snp_52286139 T C C N ADD 7583 0.119914096665338 0.111597 1.07452 0.282589 .
1 768448 snp_512562034 G A A N ADD 7916 0.223183550514231 0.11595 1.92486 0.0542468 .
1 779322 snp_54040617 A G G N ADD 7916 0.0978342643483546 0.111966 0.873763 0.382247 .
```
B_out.trd.sum - This contains GWAS summary statistics of all additive SNP effects, when the outcome is binary. V1 to V14 denote the following columns in order. Note that all .sum files follow the same structure.
* chromosome 
* base pair position 
* SNP ID 
* reference allele 
* alternate allele 
* counted allele A1 (in regression) 
* firth regression status 
* test identifier 
* number of samples in regression 
* odds ratio for A1 allele 
* standard error of log odds 
* test statistic 
* p-value  
* error code 

In addition to the output file, users can assign the function to an object and call each component in the output file separately. See the topic GWAS_binary (page 5) in manual.pdf for examples.

**Command**
```
GWEIS_binary(plink_path, "mydata", "Bpd.txt", "Bcd.txt", thread = 20, summary_output = "B_out")
```
This performs GWEIS using a logistic regression, and outputs GWEIS summary statistics of all additive and interaction SNP effects in the files "B_out.add.sum" and "B_out.gxe.sum", respectively. 

**Output**
```
V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12 V13 V14
1 754182 snp_53131969 G A A N ADD 7867 0.0957646219843793 0.143743 0.66621 0.505277 .
1 761732 snp_52286139 T C C N ADD 7583 0.0992931463425387 0.142772 0.695495 0.486745 .
1 768448 snp_512562034 G A A N ADD 7916 0.14841138439146 0.149001 0.996065 0.319219 .
1 779322 snp_54040617 A G G N ADD 7916 0.0952647242257811 0.142304 0.669477 0.503192 .
```
B_out.add.sum - This contains GWEIS summary statistics of all additive SNP effects, when the outcome is binary. 

```
V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12 V13 V14
1 754182 snp_53131969 G A A N ADDxCOVAR1 7867 0.0291703772997799 0.0941336 0.30992 0.756622 .
1 761732 snp_52286139 T C C N ADDxCOVAR1 7583 0.0235993322503244 0.0939593 0.25121 0.801651 .
1 768448 snp_512562034 G A A N ADDxCOVAR1 7916 0.0201849071590975 0.0863014 0.233834 0.815114 .
1 779322 snp_54040617 A G G N ADDxCOVAR1 7916 0.0132222001691214 0.0938807 0.140788 0.888037 .
```
B_out.gxe.sum - This contains GWEIS summary statistics of all interaction SNP effects, when the outcome is binary. 

In addition to the output files, users can assign the function to an object and call each component in the output files separately. See the topic GWEIS_binary (page 7) in manual.pdf for examples.

**Command**
```
PRS_binary(plink_path, "mydata", summary_input = "B_out.trd.sum", summary_output = "B_trd")
PRS_binary(plink_path, "mydata", summary_input = "B_out.add.sum", summary_output = "B_add")
PRS_binary(plink_path, "mydata", summary_input = "B_out.gxe.sum", summary_output = "B_gxe")
```
As explained above, “mydata” is the prefix of the PLINK format files, “B_out.trd.sum” is summary statistics generated from previous functions and used as an input for this function to construct PRS. Finally, "B_trd" is an output file name that can be defined by user. These commands compute polygenic risk scores for each individual in the target dataset and outputs the files B_trd.sscore, B_add.sscore and B_gxe.sscore respectively.

**Output**
```
#FID	IID	ALLELE_CT	NAMED_ALLELE_DOSAGE_SUM	SCORE1_AVG
1001019	1001019	19814	4772	-0.000330369
1001022	1001022	19822	4441	-0.0012421
1001035	1001035	19856	4648	-0.00277161
1001054	1001054	19868	4552	-0.00163487
```
B_trd.sscore - This contains the following columns in order.
* FID 
* IID 
* number of alleles across scored variants (ALLELE_CT)
* sum of named allele dosages (NAMED_ALLELE_DOSAGE_SUM)
* SCORE1_AVG (polygenic risk scores (PRSs), computed from the additive effects of GWAS summary statistics), of the full dataset


```
#FID	IID	ALLELE_CT	NAMED_ALLELE_DOSAGE_SUM	SCORE1_AVG
1001019	1001019	19814	4772	-0.000325805
1001022	1001022	19822	4441	-0.000530173
1001035	1001035	19856	4648	-0.00417207
1001054	1001054	19868	4552	-0.00294334
```
B_add.sscore - This contains the the following columns in order.
* FID 
* IID 
* ALLELE_CT
* NAMED_ALLELE_DOSAGE_SUM
* SCORE1_AVG (polygenic risk scores (PRSs), computed from the additive effects of GWEIS summary statistics), of the full dataset 


```
#FID	IID	ALLELE_CT	NAMED_ALLELE_DOSAGE_SUM	SCORE1_AVG
1001019	1001019	19814	4772	-0.000316507
1001022	1001022	19822	4441	-0.000688814
1001035	1001035	19856	4648	0.00188325
1001054	1001054	19868	4552	0.00115625
```
B_gxe.sscore - This contains the the following columns in order.
* FID 
* IID 
* ALLELE_CT
* NAMED_ALLELE_DOSAGE_SUM
* SCORE1_AVG (polygenic risk scores (PRSs), computed from the interaction effects of GWEIS summary statistics), of the full dataset

In addition to the output file, users can assign the function to an object and call each component in the output file separately. See the topic PRS_binary (page 11) in manual.pdf for examples.

**Command**
```
summary_regular_binary("Bpt.txt", "Bct.txt", trd_score = "B_trd.sscore", Model = 1, summary_output = "Bsummary.txt", risk_output = "Individual_risk_values.txt")
summary_regular_binary("Bpt.txt", "Bct.txt", add_score = "B_add.sscore", Model = 2, summary_output = "Bsummary.txt", risk_output = "Individual_risk_values.txt")
summary_regular_binary("Bpt.txt", "Bct.txt", add_score = "B_add.sscore", gxe_score = "B_gxe.sscore", Model = 3, summary_output = "Bsummary.txt", risk_output = "Individual_risk_values.txt")
summary_regular_binary("Bpt.txt", "Bct.txt", add_score = "B_add.sscore", gxe_score = "B_gxe.sscore", Model = 4, summary_output = "Bsummary.txt", risk_output = "Individual_risk_values.txt")
summary_regular_binary("Bpt.txt", "Bct.txt", add_score = "B_add.sscore", gxe_score = "B_gxe.sscore", Model = 5, summary_output = "Bsummary.txt", risk_output = "Individual_risk_values.txt")
```
“Bpt.txt” is binary phenotype file of the target sample, "Bct.txt" is the covariate file of the target sample, "B_trd.sscore", "B_add.sscore" and "B_gxe.sscore" are the summary statistics output files generated at GWAS and GWEIS steps. Depending on the model used, the input file should be varied. (See section $\color{red}{IMPORTANT}$ for the target models.) This function outputs 2 files. The first file, Bsummary.txt gives the summary of the fitted **regular** model for **binary** outcome. The second file, Individual_risk_values.txt contains all the calculated individual risk scores of the fitted **regular** model for **binary** outcome. 

##### Refer to the section $\color{red}{IMPORTANT}$ at the end of this document for details about models fitted at this step.

**Output**
```
Call:
glm(formula = out ~ ., family = binomial(link = logit), data = df_new)

Deviance Residuals: 
    Min       1Q   Median       3Q      Max  
-1.0629  -0.3430  -0.2324  -0.1492   3.2762  

Coefficients:
                Estimate Std. Error z value Pr(>|z|)    
(Intercept)   -6.4270500  1.3285932  -4.837 1.31e-06 ***
E              1.1326786  0.1943609   5.828 5.62e-09 ***
`E squared`   -0.3047671  0.1228849  -2.480 0.013134 *  
PRS_add        0.1121382  0.1310665   0.856 0.392229    
PRS_gxe       -0.0289305  0.1537238  -0.188 0.850722    
`PRS_gxe x E` -0.0339760  0.1254533  -0.271 0.786525    
V7             0.1264787  0.0330392   3.828 0.000129 ***
V8             0.0396320  0.0143769   2.757 0.005840 ** 
V9            -0.0198777  0.0722361  -0.275 0.783180    
V10           -0.0237545  0.0729273  -0.326 0.744630    
V11           -0.0043770  0.0719826  -0.061 0.951513    
V12            0.0269025  0.0547924   0.491 0.623435    
V13           -0.0126348  0.0244999  -0.516 0.606060    
V14           -0.0766658  0.0654078  -1.172 0.241149    
V15            0.0004197  0.0630981   0.007 0.994692    
V16           -0.0918403  0.0597228  -1.538 0.124103    
V17           -0.0280825  0.0229614  -1.223 0.221318    
V18           -0.0926273  0.0532256  -1.740 0.081810 .  
V19            0.3436252  0.2209192   1.555 0.119843    
V20            0.0349171  0.0214295   1.629 0.103230    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for binomial family taken to be 1)

    Null deviance: 793.52  on 1938  degrees of freedom
Residual deviance: 682.31  on 1919  degrees of freedom
AIC: 722.31

Number of Fisher Scoring iterations: 7
```
Bsummary.txt - This contains the target regular model summary output, when the outcome is binary. 

```
1001035 1001035 0.0144851977697558
1001105 1001105 0.0732542373614076
1001129 1001129 0.00529574201513106
1001210 1001210 0.0570084654126558
1001323 1001323 0.0106232171799868
```
Individual_risk_values.txt - This contains all the calculated individual risk scores using the target dataset (e.g. Model 5), when the outcome is binary. The columns denote the following in order.
* FID
* IID 
* estimated risk value

Note: It is recommended to fit both regular and permuted models and obtain the summary of both fitted models (using ```summary_regular_binary("Bpt.txt", "Bct.txt", add_score = "B_add.sscore", gxe_score = "B_gxe.sscore", Model = 5)``` and ```summary_permuted_binary("Bpt.txt", "Bct.txt", iterations = 1000, add_score = "B_add.sscore", gxe_score = "B_gxe.sscore")```), if you choose to fit 'PRS_gxe x E' interaction component (i.e. novel proposed model, Model 5) when generating risk scores. If the 'PRS_gxe x E' term is significant in Model 5, and insignificant in Model 5* (permuted p value), consider that the 'PRS_gxe x E' interaction component is actually insignificant (always give priority to the p value obtained from the permuted model). 

In addition to the output files, users can assign the function to an object and call each component in the output files separately. See the topic summary_regular_binary (page 18) in manual.pdf for examples.

**Command**
```
summary_permuted_binary("Bpt.txt", "Bct.txt", iterations = 1000, add_score = "B_add.sscore", gxe_score = "B_gxe.sscore")
```
This outputs the p value of the fitted **permuted** model for **binary** outcome.




###### When the outcome variable is quantitative
**Command**
```
GWAS_quantitative(plink_path, "mydata", "Qpd.txt", "Qcd.txt", thread = 20, summary_output = "Q_out")
```
As explained above, “mydata” is the prefix of the PLINK format files, “Qpd.txt” is quantitative phenotype file of the discovery sample, "Qcd.txt" is the covariate file of the discovery sample, thread indicates the number of CPUs used to run the command which can be optionally specified by the user (default is 20). Here summary_output is an output file name (also an optional argument) which can be defined by the user (default is "Q_out" and the function will generate the suffix ".trd.sum") This command performs GWAS using a linear regression, and outputs GWAS summary statistics of all additive SNP effects stored in the file "Q_out.trd.sum".

**Output**
```
V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12 V13
1 754182 snp_53131969 G A A ADD 6385 -0.108261 0.124203 -0.87164 0.383438 .
1 761732 snp_52286139 T C C ADD 6152 -0.0577641 0.124443 -0.46418 0.642535 .
1 768448 snp_512562034 G A A ADD 6426 0.328349 0.131185 2.50295 0.0123411 .
1 779322 snp_54040617 A G G ADD 6426 -0.0420902 0.123417 -0.341039 0.733085 .
```
Q_out.trd.sum - This contains GWAS summary statistics of all additive SNP effects, when the outcome is quantitative. V1 to V13 denote the following columns in order. Note that all .sum files follow the same structure.
* chromosome 
* base pair position 
* SNP ID 
* reference allele 
* alternate allele 
* counted allele A1 (in regression) 
* test identifier 
* number of samples in regression 
* odds ratio for A1 allele 
* standard error of log odds 
* test statistic 
* p-value  
* error code 

In addition to the output file, users can assign the function to an object and call each component in the output file separately. See the topic GWAS_quantitative (page 6) in manual.pdf for examples.

**Command**
```
GWEIS_quantitative(plink_path, "mydata", "Qpd.txt", "Qcd.txt", thread = 20, summary_output = "Q_out")
```
This performs GWEIS using a linear regression, and outputs GWEIS summary statistics of all additive and interaction SNP effects in the files "Q_out.add.sum" and "Q_out.gxe.sum", respectively. 

**Output**
```
V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12 V13
1 754182 snp_53131969 G A A ADD 6385 -0.108204 0.124284 -0.870618 0.383995 .
1 761732 snp_52286139 T C C ADD 6152 -0.0579637 0.124507 -0.465544 0.641558 .
1 768448 snp_512562034 G A A ADD 6426 0.32818 0.131207 2.50124 0.0124006 .
1 779322 snp_54040617 A G G ADD 6426 -0.0426272 0.123464 -0.345261 0.729909 .
```
Q_out.add.sum - This contains GWEIS summary statistics of all additive SNP effects, when the outcome is quantitative. 

```
V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12 V13
1 754182 snp_53131969 G A A ADDxCOVAR1 6385 -0.0362558 0.124498 -0.291216 0.770895 .
1 761732 snp_52286139 T C C ADDxCOVAR1 6152 -0.030858 0.125198 -0.246473 0.805324 .
1 768448 snp_512562034 G A A ADDxCOVAR1 6426 -0.00407146 0.129867 -0.031351 0.974991 .
1 779322 snp_54040617 A G G ADDxCOVAR1 6426 -0.00553396 0.123911 -0.0446607 0.964379 .
```
Q_out.gxe.sum - This contains GWEIS summary statistics of all interaction SNP effects, when the outcome is quantitative. 

In addition to the output files, users can assign the function to an object and call each component in the output files separately. See the topic GWEIS_quantitative (page 9) in manual.pdf for examples.

**Command**
```
PRS_quantitative(plink_path, "mydata", summary_input = "Q_out.trd.sum", summary_output = "Q_trd")
PRS_quantitative(plink_path, "mydata", summary_input = "Q_out.add.sum", summary_output = "Q_add")
PRS_quantitative(plink_path, "mydata", summary_input = "Q_out.gxe.sum", summary_output = "Q_gxe")
```
As explained above, “mydata” is the prefix of the PLINK format files, “Q_out.trd.sum” is summary statistics generated from previous functions and used as an input for this function to construct PRS. Finally, "Q_trd" is an output file name that can be defined by user. These commands compute polygenic risk scores for each individual in the target dataset and outputs the files Q_trd.sscore, Q_add.sscore and Q_gxe.sscore respectively.

**Output**
```
#FID	IID	ALLELE_CT	NAMED_ALLELE_DOSAGE_SUM	SCORE1_AVG
1001019	1001019	19814	4772	0.000688852
1001022	1001022	19822	4441	0.00142371
1001035	1001035	19856	4648	-0.000863821
1001054	1001054	19868	4552	0.000226732
```
Q_trd.sscore - This contains the following columns in order.
* FID 
* IID 
* number of alleles across scored variants (ALLELE_CT)
* sum of named allele dosages (NAMED_ALLELE_DOSAGE_SUM)
* SCORE1_AVG (polygenic risk scores (PRSs), computed from the additive effects of GWAS summary statistics), of the full dataset

```
#FID	IID	ALLELE_CT	NAMED_ALLELE_DOSAGE_SUM	SCORE1_AVG
1001019	1001019	19814	4772	0.000629705
1001022	1001022	19822	4441	0.00142161
1001035	1001035	19856	4648	-0.000898764
1001054	1001054	19868	4552	0.000256159
```
Q_add.sscore - This contains the the following columns in order.
* FID 
* IID 
* ALLELE_CT
* NAMED_ALLELE_DOSAGE_SUM
* SCORE1_AVG (polygenic risk scores (PRSs), computed from the additive effects of GWEIS summary statistics), of the full dataset 

```
#FID	IID	ALLELE_CT	NAMED_ALLELE_DOSAGE_SUM	SCORE1_AVG
1001019	1001019	19814	4772	0.00121497
1001022	1001022	19822	4441	0.00156228
1001035	1001035	19856	4648	0.000543463
1001054	1001054	19868	4552	0.000137235
```
Q_gxe.sscore - This contains the the following columns in order.
* FID 
* IID 
* ALLELE_CT
* NAMED_ALLELE_DOSAGE_SUM
* SCORE1_AVG (polygenic risk scores (PRSs), computed from the interaction effects of GWEIS summary statistics), of the full dataset


In addition to the output file, users can assign the function to an object and call each component in the output file separately. See the topic PRS_quantitative (page 12) in manual.pdf for examples.

**Command**
```
summary_regular_quantitative("Qpt.txt", "Qct.txt", trd_score = "Q_trd.sscore", Model = 1, summary_output = "Qsummary.txt", risk_output = "Individual_risk_values.txt")
summary_regular_quantitative("Qpt.txt", "Qct.txt", add_score = "Q_add.sscore", Model = 2, summary_output = "Qsummary.txt", risk_output = "Individual_risk_values.txt")
summary_regular_quantitative("Qpt.txt", "Qct.txt", add_score = "Q_add.sscore", gxe_score = "Q_gxe.sscore", Model = 3, summary_output = "Qsummary.txt", risk_output = "Individual_risk_values.txt")
summary_regular_quantitative("Qpt.txt", "Qct.txt", add_score = "Q_add.sscore", gxe_score = "Q_gxe.sscore", Model = 4, summary_output = "Qsummary.txt", risk_output = "Individual_risk_values.txt")
```
“Qpt.txt” is quantitative phenotype file of the target sample, "Qct.txt" is the covariate file of the target sample, "Q_trd.sscore", "Q_add.sscore" and "Q_gxe.sscore" are the summary statistics output files generated at GWAS and GWEIS steps. Depending on the model used, the input file should be varied. (See section $\color{red}{IMPORTANT}$ for the target models.) This function outputs 2 files. The first file, Qsummary.txt gives the summary of the fitted **regular** model for **quantitative** outcome. The second file, Individual_risk_values.txt contains all the calculated individual risk scores of the fitted **regular** model for **quantitative** outcome. 

##### Refer to the section $\color{red}{IMPORTANT}$ at the end of this document for details about models fitted at this step.

**Output**
```
Call:
lm(formula = out ~ ., data = df_new)

Residuals:
    Min      1Q  Median      3Q     Max 
-2.1863 -0.6602 -0.1432  0.4995  6.0425 

Coefficients:
               Estimate Std. Error t value Pr(>|t|)    
(Intercept)   -0.074259   0.283690  -0.262 0.793540    
E              0.001652   0.028315   0.058 0.953492    
PRS_add       -0.034537   0.025300  -1.365 0.172414    
PRS_gxe        0.026966   0.025168   1.071 0.284140    
`PRS_gxe x E`  0.027421   0.027941   0.981 0.326565    
V6             0.012552   0.008906   1.409 0.158920    
V7             0.003991   0.003202   1.246 0.212841    
V8            -0.005173   0.016407  -0.315 0.752575    
V9             0.027789   0.017174   1.618 0.105835    
V10           -0.022333   0.016265  -1.373 0.169920    
V11           -0.013725   0.012155  -1.129 0.259012    
V12            0.011525   0.005596   2.060 0.039603 *  
V13            0.005250   0.015875   0.331 0.740885    
V14           -0.035308   0.014171  -2.492 0.012820 *  
V15           -0.030871   0.014077  -2.193 0.028455 *  
V16           -0.005749   0.005695  -1.009 0.312960    
V17           -0.007022   0.012450  -0.564 0.572838    
V18            0.170381   0.050142   3.398 0.000696 ***
V19           -0.026555   0.005044  -5.265  1.6e-07 ***
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Residual standard error: 0.9842 on 1560 degrees of freedom
Multiple R-squared:  0.04241,	Adjusted R-squared:  0.03136 
F-statistic: 3.838 on 18 and 1560 DF,  p-value: 1.004e-07
```
Qsummary.txt - This contains the target regular model summary output, when the outcome is quantitative. 

```
1001054 1001054 0.110000856152333
1001090 1001090 -0.024816480821851
1001203 1001203 0.234383736636126
1001450 1001450 0.117064079395998
1001516 1001516 0.224683433022139
```
Individual_risk_values.txt - This contains all the calculated individual risk scores using the target dataset (e.g. Model 4), when the outcome is quantitative. The columns denote the following in order.
* FID
* IID 
* estimated risk value

Note: It is recommended to fit both regular and permuted models and obtain the summary of both fitted models (using ```summary_regular_quantitative("Qpt.txt", "Qct.txt", add_score = "Q_add.sscore", gxe_score = "Q_gxe.sscore", Model = 4)``` and ```summary_permuted_quantitative("qpt.txt", "qct.txt", iterations = 1000, add_score = "Q_add.sscore", gxe_score = "Q_gxe.sscore")```), if you choose to fit 'PRS_gxe x E' interaction component (i.e. novel proposed model, Model 4) when generating risk scores. If the 'PRS_gxe x E' term is significant in Model 4, and insignificant in Model 4* (permuted p value), consider that the 'PRS_gxe x E' interaction component is actually insignificant (always give priority to the p value obtained from the permuted model). 

In addition to the output files, users can assign the function to an object and call each component in the output files separately. See the topic summary_regular_quantitative (page 17) in manual.pdf for examples.

**Command**
```
summary_permuted_quantitative("Qpt.txt", "Qct.txt", iterations = 1000, add_score = "Q_add.sscore", gxe_score = "Q_gxe.sscore")
```
This outputs the p value of the fitted **permuted** model for **quantitative** outcome.



## $$\color{red}{IMPORTANT}$$
The discovery model used in ```GWAS_binary(plink_path, "mydata", "Bpd.txt", "Bcd.txt", thread = 20, summary_output = "B_trd.sum")``` or ```GWAS_quantitative(plink_path, "mydata", "Qpd.txt", "Qcd.txt", thread = 20, summary_output = "Q_trd.sum")``` is as follows:
* y = b_trd.W + error
 where y is the outcome variable, b_trd is the estimated SNP effect and W is the SNP genotype.
 
The discovery model used in ```GWEIS_binary(plink_path, "mydata", "Bpd.txt", "Bcd.txt", thread = 20, summary_output = "B_trd.sum")``` is as follows:
 * y = b_add.W + b_cov.E + b_cov2.E^2 + b_gxe.(WxE) + error
where y is the outcome variable, b_add is the estimated additive SNP effect, E is the covariate, W is the SNP genotype, b_cov is the estimated effect of the covariate, b_cov2 is the estimated effect of the squared covariate and b_gxe is the estimated effect of the Hadamard product of WxE.

The discovery model used in ```GWEIS_quantitative(plink_path, "mydata", "Qpd.txt", "Qcd.txt", thread = 20, summary_output = "Q_trd.sum")``` is as follows:
 * y = b_add.W + b_cov.E + b_gxe.(WxE) + error
where y is the outcome variable, b_add is the estimated additive SNP effect, E is the covariate, W is the SNP genotype, b_cov is the estimated effect of the covariate and b_gxe is the estimated effect of the Hadamard product of WxE.


The fitted (target) models in ```summary_regular_binary("Bpt.txt", "Bct.txt", trd_score = "B_trd.sscore", add_score = "B_add.sscore", gxe_score = "B_gxe.sscore", Model)``` or ```summary_regular_quantitative("Qpt.txt", "Qct.txt", trd_score = "Q_trd.sscore", add_score = "Q_add.sscore", gxe_score = "Q_gxe.sscore", Model)``` are as follows:

* Model 1: y = PRS_trd + E + PRS_trd x E + confounders + error
* Model 2: y = PRS_add + E + PRS_add x E + confounders + error
* Model 3: y = PRS_add + E + PRS_gxe x E + confounders + error
* Model 4: y = PRS_add + E + PRS_gxe + PRS_gxe x E + confounders + error
* Model 4*: permuted Model 4
* Model 5: y = PRS_add + E + E^2 + PRS_gxe + PRS_gxe x E + confounders + error
* Model 5*: permuted Model 5

where y is the outcome variable, E is the covariate of interest, PRS_trd and PRS_add are the polygenic risk scores computed using additive SNP effects of GWAS and GWEIS summary statistics respectively, and PRS_gxe is the polygenic risk scores computed using GxE interaction SNP effects of GWEIS summary statistics.

When deciding on the number of permutations to be used, always get an idea from the p value obtained from either Model 4 or 5 (accordingly). If that p value is insignificant, you can use any number of permututations (e.g. 1000), but if that p value is highly significant (p value is very small and very close to zero), it is recommended to increase the number of permutations for better results. For example assume that the p value obtained from Model 4 is 0.000154, then it is recommended to use at least 10,000 iterations in the permutation model.


# Contact 
Please contact Dovini Jayasinghe (dovini.jayasinghe@mymail.unisa.edu.au) or Md Moksedul Momin (md_moksedul.momin@mymail.unisa.edu.au) or Hong Lee (hong.lee@unisa.edu.au) for queries.
