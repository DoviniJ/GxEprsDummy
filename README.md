---
title: "GxEprs"
author: "Dovini Jayasinghe, Md Moksedul Momin and Hong Lee"
date: "01/05/2023" #to be edited
output: pdf_document
---

# GxEprs
The 'GxEprs' is an R package for **genomic prediction** that uses a sophisticated method that has been enhanced for its prediction accuracy. It performs Genome Wide Association Studies (GWAS) and Genome Wide Environment Interaction Studies (GWEIS) using plink2 which is easily accessible from within R.

# Data preparation

## File formats
### Input files
1) mydata.fam - This is one of the binary files which contains the following columns in order. The example dataset has 10,000 individuals. Note that the file has no column headings.
* family ID (FID) 
* individual ID (IID) 
* father's ID 
* mother's ID 
* sex 
* phenotype value

<!--- ![image](https://user-images.githubusercontent.com/131835334/236634478-99a10ec5-2e05-4259-981e-d67562b1a06a.png) -->
```
1001019 1001019 0 0 2 -9
1001022 1001022 0 0 2 -9
1001035 1001035 0 0 1 -9
1001054 1001054 0 0 2 -9
1001078 1001078 0 0 2 -9
```
  
2) mydata.bim - This is one of the binary files which contains the following columns in order. The example dataset has 10,000 SNPs. Note that the file has no column headings. 
* chromosome code 
* SNP ID 
* position of centimorgans 
* base-pair coordinate 
* minor allele  
* reference allele 

<!--- ![image](https://user-images.githubusercontent.com/131835334/236634694-5dbe6a29-5ae0-44c9-b076-b80fcabb7144.png) -->
```
1	  snp_53131969	0	754182	A	G
1	  snp_52286139	0	761732	C	T
1	  snp_512562034	0	768448	A	G
1	  snp_54040617	0	779322	G	A
1	  snp_52980300	0	785989	T	C
```

3) mydata.bed - This is also a binary file which cannot be read by humans.
4) Bpd.txt - This is a .txt file which contains the following columns in order. The discovery dataset has 7916 individuals. Note that the file has no column headings.
* FID 
* IID  
* binary phenotype (1=controls, 2=cases) of the discovery sample

<!--- ![image](https://user-images.githubusercontent.com/131835334/236635016-88560176-a22a-4863-b200-4ddca8ca6980.png) -->
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
* square of the standardized covariate  
* 14 confounders of the discovery sample

<!--- ![image](https://user-images.githubusercontent.com/131835334/236635276-7e1c6d92-3a84-4f9e-b68b-a171d9684da3.png) -->
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

<!--- ![image](https://user-images.githubusercontent.com/131835334/236635388-53c3ff05-ae8b-498c-8354-0e1419aaf56f.png) -->
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
* square of the standardized covariate  
* 14 confounders of the target sample

<!--- ![image](https://user-images.githubusercontent.com/131835334/236635508-d08cec8f-38bb-4008-8bc5-18ad87c8eece.png) -->
```
1001035 1001035 -0.420822976931972 0.177091977913887 -4.12263 57 -13.5185 5.40198 -4.81994 0.664494 -4.92217 -0.451329 3.14677 0.42704 0.821306 -2.77705 1 7
1001105 1001105 -0.0805583280660987 0.00648964422080519 2.92534 56 -13.6236 3.21643 -0.856048 0.750187 -2.01798 -0.350832 5.10141 2.1807 -6.04343 1.78928 1 19
1001129 1001129 -1.32752644870189 1.76232647200306 -3.09118 61 -9.94475 3.60562 -0.917639 0.905664 -5.09843 -1.16329 -1.88102 -1.24154 0.699574 2.2442 0 20
1001210 1001210 0.698007239555549 0.487214106471958 4.58829 49 -12.5471 4.09467 -2.58951 6.06898 12.9822 -0.704179 2.90357 -0.334968 5.04274 0.66175 0 10
1001323 1001323 -0.657981606980219 0.432939795124272 -3.53948 56 -12.795 2.91524 -2.72794 3.61555 3.92957 -2.93899 -0.454737 2.31013 2.51783 -4.15592 0 7
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
headQBcov_target) #covariate file of the target sample when the outcome is quantitative
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
Note that the step 3.1.3 described above is to call the embedded data files in this package itself. However, when the users have to call their own data, they can follow the same approach. It is more convenient if the users can store all their data files in the same working directory. For example, assume that the file names are as follows 
(Refer to 'File Formats' section of this document to view the format details of each of the following input file): 
* binary files: mydata.fam, mydata.bim and mydata.bed
* phenotype file of discovery sample (binary outcome): Bpd.txt 
* covariate file of discovery sample (binary outcome): Bcd.txt 
* phenotype file of target sample (binary outcome): Bpt.txt
* covariate file of the target sample (binary outcome): Bct.txt
* phenotype file of discovery sample (quantitative outcome): Qpd.txt 
* covariate file of discovery sample (quantitative outcome): Qcd.txt 
* phenotype file of target sample (quantitative outcome): Qpt.txt
* covariate file of the target sample (quantitative outcome): Qct.txt

Note that, all these files can be placed in a separate location. It is always upto the users choice. In that case remember to give the full path to the file location since R identifies files by name, only when they are in the same directory.
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

###### Step 3.1.4 Set the number of confounders as 14 and number of threads (CPUs) as 20
```
n_confounders = 14 #this is the number of confounders in the covariate files of the embedded example datasets. Users can change this value according to the number of confounders they use in their covariate files.
thread = 20 #this is the number of threads specified for this example. Users can change this value according to their preference.
```
###### Step 3.1.5 Set the working directory and run the following R functions in the given order
```
setwd("<path to working directory>") #set the working directory where you need to save the output files
```
When the outcome variable is binary
```
GWAS_binary(plink_path, mydata, "Bpd.txt", "Bcd.txt", n_confounders, thread)
GWEIS_binary(plink_path, mydata, "Bpd.txt", "Bcd.txt", n_confounders, thread)
PRS_binary(plink_path, mydata)
summary_regular_binary("Bpt.txt", "Bct.txt", n_confounders)
```
When the outcome variable is quantitative (NOT ADDED TO THE PACKAGE YET)
```
GWAS_quantitative(mydata, "Qpd.txt", "Qcd.txt", n_confounders, thread)
GWEIS_quantitative(mydata, "Qpd.txt", "Qcd.txt", n_confounders, thread)
PRS_quantitative(plink_path, mydata)
summary_regular_quantitative("Qpt.txt", "Qct.txt", n_confounders)
```

Here, the fitted models in ```summary_regular_binary(Bphe_target, Bcov_target, n_confounders)``` or ```summary_regular_quantitative(Qphe_target, Qcov_target, n_confounders)``` are as follows:

* Model 1: y = PRS_trd + E + PRS_trd x E + confounders
* Model 2: y = PRS_add + E + PRS_add x E + confounders
* Model 3: y = PRS_add + E + PRS_gxe x E + confounders
* Model 4: y = PRS_add + E + PRS_gxe + PRS_gxe x E + confounders
* Model 4*: permuted Model 4
* Model 5: y = PRS_add + E + E^2 + PRS_gxe + PRS_gxe x E + confounders
* Model 5*: permuted Model 5

where y is the outcome variable, E is the covariate of interest, PRS_trd and PRS_add are the polygenic risk scores computed using additive SNP effects of GWAS and GWEIS summary statistics respectively, and PRS_gxe is the polygenic risk scores computed using GxE interaction SNP effects of GWEIS summary statistics.


##### Step 3.2: Run the following code line to check the significance of the interaction term ('PRS_gxe x E'):

When the outcome variable is binary:
```
summary_permuted_binary(Bphe_target, Bcov_target, n_confounders)
```
Note: It is recommended to fit both regular and permuted models and obtain the summary of both fitted models (using ```summary_regular_binary(Bphe_target, Bcov_target, n_confounders)``` and ```summary_permuted_binary(Bphe_target, Bcov_target, n_confounders)```), if you choose to fit 'PRS_gxe x E' interaction component (i.e. novel proposed model, Model 5) when generating risk scores. If the 'PRS_gxe x E' term is significant in Model 5, and insignificant in Model 5* (permuted p value), consider that the 'PRS_gxe x E' interaction component is actually insignificant (always give priority to the p value obtained from the permuted model). 


When the outcome variable is quantitative:
```
summary_permuted_quantitative(Qphe_target, Qcov_target, n_confounders)
```

Note: It is recommended to fit both regular and permuted models and obtain the summary of both fitted models (using ```summary_regular_quantitative(Qphe_target, Qcov_target, n_confounders)``` and ```summary_permuted_quantitative(Qphe_target, Qcov_target, n_confounders)```), if you choose to fit 'PRS_gxe x E' interaction component (i.e. novel proposed model, Model 5) when generating risk scores. If the 'PRS_gxe x E' term is significant in Model 4, and insignificant in Model 4* (permuted p value), consider that the 'PRS_gxe x E' interaction component is actually insignificant (always give priority to the p value obtained from the permuted model). 





### Output files
1) B_trd.sum - This contains GWAS summary statistics of all additive SNP effects, when the outcome is binary. V1 to V14 denotes the following columns in order. Note that all .sum files follow the same structure.
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

<!--- ![image](https://user-images.githubusercontent.com/131835334/236993782-75e2d666-5b83-4348-80fc-994801c406a4.png) -->
```
V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12 V13 V14
1 754182 snp_53131969 G A A N ADD 7867 0.10380130071374 0.112336 0.923983 0.355495 .
1 761732 snp_52286139 T C C N ADD 7583 0.119914096665338 0.111597 1.07452 0.282589 .
1 768448 snp_512562034 G A A N ADD 7916 0.223183550514231 0.11595 1.92486 0.0542468 .
1 779322 snp_54040617 A G G N ADD 7916 0.0978342643483546 0.111966 0.873763 0.382247 .
```

2) B_add.sum - This contains GWEIS summary statistics of all additive SNP effects, when the outcome is binary. 
<!--- ![image](https://user-images.githubusercontent.com/131835334/236993906-1ea97b12-af9e-4693-96a2-aff2128d1eb7.png) -->
```
V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12 V13 V14
1 754182 snp_53131969 G A A N ADD 7867 0.0957646219843793 0.143743 0.66621 0.505277 .
1 761732 snp_52286139 T C C N ADD 7583 0.0992931463425387 0.142772 0.695495 0.486745 .
1 768448 snp_512562034 G A A N ADD 7916 0.14841138439146 0.149001 0.996065 0.319219 .
1 779322 snp_54040617 A G G N ADD 7916 0.0952647242257811 0.142304 0.669477 0.503192 .
```

3) B_gxe.sum - This contains GWEIS summary statistics of all interaction SNP effects, when the outcome is binary. 
<!--- ![image](https://user-images.githubusercontent.com/131835334/236993968-f07a1493-4d11-494e-b7f5-8b4747641207.png) -->
```
V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12 V13 V14
1 754182 snp_53131969 G A A N ADDxCOVAR1 7867 0.0291703772997799 0.0941336 0.30992 0.756622 .
1 761732 snp_52286139 T C C N ADDxCOVAR1 7583 0.0235993322503244 0.0939593 0.25121 0.801651 .
1 768448 snp_512562034 G A A N ADDxCOVAR1 7916 0.0201849071590975 0.0863014 0.233834 0.815114 .
1 779322 snp_54040617 A G G N ADDxCOVAR1 7916 0.0132222001691214 0.0938807 0.140788 0.888037 .
```

4) B_trd.sscore - This contains the following columns in order.
* FID 
* IID 
* number of alleles across scored variants (ALLELE_CT)  
* polygenic risk scores (PRSs), computed from the additive effects of GWAS summary statistics, of the full dataset

<!--- ![image](https://user-images.githubusercontent.com/131835334/236994019-1ef3609a-4142-4fda-a89b-e05c81fc6d32.png) -->
```
#FID	IID	ALLELE_CT	NAMED_ALLELE_DOSAGE_SUM	SCORE1_AVG
1001019	1001019	19814	4772	-0.000330369
1001022	1001022	19822	4441	-0.0012421
1001035	1001035	19856	4648	-0.00277161
1001054	1001054	19868	4552	-0.00163487
```

5) B_add.sscore - This contains the the following columns in order.
* FID 
* IID 
* ALLELE_CT  
* polygenic risk scores (PRSs), computed from the additive effects of GWEIS summary statistics, of the full dataset 

<!--- ![image](https://user-images.githubusercontent.com/131835334/236994081-d346ae48-d22a-4a35-a608-e4ed7535ec6c.png) -->
```
#FID	IID	ALLELE_CT	NAMED_ALLELE_DOSAGE_SUM	SCORE1_AVG
1001019	1001019	19814	4772	-0.000325805
1001022	1001022	19822	4441	-0.000530173
1001035	1001035	19856	4648	-0.00417207
1001054	1001054	19868	4552	-0.00294334
```

6) B_gxe.sscore - This contains the the following columns in order.
* FID 
* IID 
* ALLELE_CT  
* polygenic risk scores (PRSs), computed from the interaction effects of GWEIS summary statistics, of the full dataset

<!--- ![image](https://user-images.githubusercontent.com/131835334/236994128-e8e58f0e-6e0c-4494-ad58-5dc9f765f6e9.png) -->
```
#FID	IID	ALLELE_CT	NAMED_ALLELE_DOSAGE_SUM	SCORE1_AVG
1001019	1001019	19814	4772	-0.000316507
1001022	1001022	19822	4441	-0.000688814
1001035	1001035	19856	4648	0.00188325
1001054	1001054	19868	4552	0.00115625
```

7) Bsummary.txt - This contains the target regular model summary output, when the outcome is binary. 
<!--- ![image](https://user-images.githubusercontent.com/131835334/236994166-c7abfafc-51e2-40c0-a240-6715aa04a457.png) -->
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

8) B_permuted_p.txt - This contains the p-value of the permuted model, when the outcome is binary. Since the permutation is random at each time, we did not include an example here.
9) Individual_risk_values.txt - This contains all the calculated individual risk scores using the target model (either regular model or permuted model), when the outcome is binary. The columns denote the following in order.
* FID
* IID 
* estimated risk value
 
For demonstration, we include the results of regular model below. 

<!--- ![image](https://user-images.githubusercontent.com/131835334/236994221-e61723a3-23f5-4e75-9144-252f3950795b.png) -->
```
1001035 1001035 0.014485064706852
1001105 1001105 0.0732548491904503
1001129 1001129 0.00529571917986906
1001210 1001210 0.057008389028729
1001323 1001323 0.0106231802372828
```


## Tasks of each function
1) GWAS_binary - This performs GWAS and outputs the file B_trd.sum which contains GWAS summary statistics of all additive SNP effects
2) GWEIS_binary - This performs GWEIS and outputs the files B_add.sum and B_gxe.sum which contain GWEIS summary statistics of all additive and interaction SNP effects
3) PRS_binary - This computes polygenic risk scores for each individual in the target dataset and outputs the files B_trd.sscore, B_trd.sscore and B_trd.sscore  
4) summary_regular_binary - This outputs the file Bsummary.txt which gives the summary of the fitted **regular** model for **binary** outcome
5) summary_permuted_binary - This outputs the file Bsummary.txt which gives the summary of the fitted **permuted** model for **binary** outcome
6) results_regular_binary - This outputs the file Individual_risk_values.txt containing all the calculated individual risk scores using the **regular** genomic prediction model for **binary** outcome
7) results_permuted_binary - This outputs the file Individual_risk_values.txt containing all the calculated individual risk scores using the **permuted** genomic prediction model for **binary** outcome


# Contact 
Please contact Dovini Jayasinghe (dovinij@gmail.com) or Moksedul Momin (md_moksedul.momin@mymail.unisa.edu.au) or Hong Lee (hong.lee@unisa.edu.au) or for queries.
