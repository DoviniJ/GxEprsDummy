#' GWEIS_binary function
#' 
#' This function performs GWEIS using plink2 and outputs the GWEIS summary statistics files with additive SNP effects named B_add.sum and interaction SNP effects named B_gxe.sum
#' 
#' @param b_file Prefix of the binary files, where all .fam, .bed and .bim files have a common prefix
#' @param pheno_file Name (with file extension) of the phenotype file containing family ID, individual ID and phenotype of the discovery dataset as columns, without heading
#' @param covar_file Name (with file extension) of the covariate file containing family ID, individual ID, standardized covariate, square of standardized covariate, and/or confounders of the discovery dataset as columns, without heading
#' @param n_confounders Number of confounding variables in the discovery dataset
#' @param thread Number of threads used
#' 
#' @keywords gwies, interaction, gxe
#' 
#' @export 
#' 
#' @importFrom 
#' 
#' @return This function will perform GWEIS and output
#' \item{B_add.sum} GWEIS summary statistics file with additive SNP effects
#' \item{B_gxe.sum} GWEIS summary statistics file with interaction SNP effects
#' 
#' @example 
#' GWEIS_binary(DummyData, Bphe_discovery.txt, Bcov_discovery.txt, 14, 20)
#' GWEIS_binary(DummyData, Bphe_discovery.txt, Bcov_discovery.txt, 0, 20)


GWEIS_binary <- function(b_file, pheno_file, covar_file, n_confounders, thread){
  
  parameters <- c(1, 2, 3, (1:(n_confounders+1))+3)
  param_vec <- paste0(parameters, collapse = ", ")
  
  system(paste0("./plink2 --bfile ", b_file, 
                " --glm interaction --pheno ", 
                pheno_file, 
                " --covar ", covar_file, 
                " --parameters ", param_vec, 
                " --allow-no-sex --threads ", 
                thread,
                " --out B_gweis"))
  
  plink_output <- read.table("B_gweis.PHENO1.glm.logistic.hybrid", header = FALSE)
  filtered_output <- as.data.frame(plink_output[(plink_output$V8=="ADD"),])
  filtered_output$V10 = log(filtered_output$V10)
  filtered_output2 <- plink_output[(plink_output$V8=="ADDxCOVAR1"),]
  filtered_output2$V10 <- log(filtered_output2$V10)
  
  sink("B_add.sum")
  write.table(filtered_output, sep = " ", row.names = FALSE, quote = FALSE)
  sink()

  sink("B_gxe.sum")
  write.table(filtered_output2, sep = " ", row.names = FALSE, quote = FALSE)
  sink()
  
}