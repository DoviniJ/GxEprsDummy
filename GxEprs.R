GWAS_binary <- function(b_file, pheno_file, covar_file, n_confounders, thread){
  
  #Pre-define fixed arguments
  parameters <- c(1, (1:n_confounders)+3) #1  4  5  6  7  8  9 10 11 12 13 14 15 16 17 when n_confounders = 14
  param_vec <- paste0(parameters, collapse = ", ")
  
  # Run plink analysis
  system(paste0("./plink2 --bfile ", b_file, 
                " --glm --pheno ", 
                pheno_file, 
                " --covar ", covar_file, 
                " --parameters ", param_vec, 
                " --allow-no-sex --threads ", 
                thread,
                " --out B_gwas"))
  
  # Filter plink output
  plink_output <- read.table("B_gwas.PHENO1.glm.logistic.hybrid", header = FALSE)
  filtered_output <- plink_output[(plink_output$V8=="ADD"),]
  filtered_output$V10 = log(filtered_output$V10)
  
  #Write B_trd.sum file
  sink("B_trd.sum")
  write.table(filtered_output, sep = " ", row.names = FALSE, quote = FALSE)
  sink()
}

GWEIS_binary <- function(b_file, pheno_file, covar_file, n_confounders, thread){
  
  #Pre-define fixed arguments
  parameters <- c(1, 2, 3, (1:(n_confounders+1))+3) #1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 when n_confounders = 14
  param_vec <- paste0(parameters, collapse = ", ")
  
  # Run plink analysis
  system(paste0("./plink2 --bfile ", b_file, 
                " --glm interaction --pheno ", 
                pheno_file, 
                " --covar ", covar_file, 
                " --parameters ", param_vec, 
                " --allow-no-sex --threads ", 
                thread,
                " --out B_gweis"))
  
  # Filter plink output
  plink_output <- read.table("B_gweis.PHENO1.glm.logistic.hybrid", header = FALSE)
  filtered_output <- as.data.frame(plink_output[(plink_output$V8=="ADD"),])
  filtered_output$V10 = log(filtered_output$V10)
  filtered_output2 <- plink_output[(plink_output$V8=="ADDxCOVAR1"),]
  filtered_output2$V10 <- log(filtered_output2$V10)
  
  #Write B_add.sum file
  sink("B_add.sum")
  write.table(filtered_output, sep = " ", row.names = FALSE, quote = FALSE)
  sink()
  
  #Write B_gxe.sum file
  sink("B_gxe.sum")
  write.table(filtered_output2, sep = " ", row.names = FALSE, quote = FALSE)
  sink()
}

PRS_binary <- function(b_file){
  # Run plink code to obtain polygenic risk scores
  system(paste0("./plink2 --bfile ", b_file, 
                " --score B_trd.sum 3 6 10 header --out B_trd"))
  system(paste0("./plink2 --bfile ", b_file, 
                " --score B_add.sum 3 6 10 header --out B_add"))
  system(paste0("./plink2 --bfile ", b_file, 
                " --score B_gxe.sum 3 6 10 header --out B_gxe"))
}

results_binary <- function(n_confounders){
  fam=read.table("Bphe_target.txt",header=F) 
  colnames(fam) <- c("FID", "IID", "PHENOTYPE")
  dat=read.table("Bcov_target.txt",header=F)
  colnames(dat)[1] <- "FID"
  colnames(dat)[2] <- "IID"
  prs0_all=read.table("B_trd.sscore")
  colnames(prs0_all)[1] <- "FID"
  colnames(prs0_all)[2] <- "IID"
  prs0=merge(fam, prs0_all, by = "FID")
  prs1_all=read.table("B_add.sscore")
  colnames(prs1_all)[1] <- "FID"
  colnames(prs1_all)[2] <- "IID"
  prs1=merge(fam, prs1_all, by = "FID")
  prs2_all=read.table("B_gxe.sscore")
  colnames(prs2_all)[1] <- "FID"
  colnames(prs2_all)[2] <- "IID"
  prs2=merge(fam, prs2_all, by = "FID")
  
  m1 <- match(prs0$IID.x, dat$IID)
  
  out = fam$PHENOTYPE[m1]
  
  cov=dat$V3[m1]
  ps0=scale(prs0$V5)
  ps1=scale(prs1$V5)
  ps2=scale(prs2$V5)
  xv0=scale(prs0$V5*cov)
  xv1=scale(prs1$V5*cov)
  xv2=scale(prs2$V5*cov)
  
  cov2=scale(cov^2)
  
  if(n_confounders == 0){
    m = glm(out ~ cov + cov2 + ps1 + ps2 + xv2, family = binomial(link = logit))
    m_fit <- fitted.values(m)
    sink("Individual_risk_values_without_confounders.txt")
    write.table(cbind(fam$FID[m1], fam$IID[m1], m_fit), row.names = F, col.names = F, quote = F)
    sink()
  }else{
    conf_var <- matrix(ncol = n_confounders, nrow = nrow(dat))
    for (k in 1:n_confounders) {
      conf_var[, k] <- as.numeric(dat[, k+4])
    }
    conf_var <- conf_var[m1,]
    df_new <- as.data.frame(cbind(out, cov, cov2, ps1, ps2, xv2, conf_var))
    m = glm(out ~., data = df_new, family = binomial(link = logit))
    m_fit <- fitted.values(m)
    sink("Individual_risk_values.txt")
    write.table(cbind(fam$FID[m1], fam$IID[m1], m_fit), row.names = F, col.names = F, quote = F)
    sink()
  }
  
}




######################################################################################## 
########################################################################################

DummyData <- "DummyData" #binary file name (DummyData.fam, DummyData.bim, DummyData.bed)
Bphe_discovery <- "Bphe_discovery.txt" #File format: FID, IID, PHE (1 = controls, 2 = cases)
Bcov_discovery <- "Bcov_discovery.txt" #File format: FID, IID, scaled_COV, scaled_COV_sq, all the confounders as separate columns
Bphe_target <- "Bphe_target.txt" #File format: FID, IID, PHE (0 = controls, 1 = cases)
Bcov_target <- "Bcov_target.txt" #File format: FID, IID, scaled_COV, scaled_COV_sq, all the confounders as separate columns

#Test GWAS_binary
GWAS_binary(DummyData, Bphe_discovery, Bcov_discovery, 14, 20)

#Test GWEIS_binary
GWEIS_binary(DummyData, Bphe_discovery, Bcov_discovery, 14, 20)

#Test PRS_binary
PRS_binary(DummyData)


#Test results_binary
results_binary(0)
results_binary(14)



