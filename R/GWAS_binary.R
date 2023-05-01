#' GWAS_binary function
#' 
#' This function performs GWAS using plink2 and outputs the GWAS summary statistics file with additive SNP effects named B_trd.sum
#' 
#' @param b_file Prefix of the binary files, where all .fam, .bed and .bim files have a common prefix
#' @param pheno_file Name (with file extension) of the phenotype file containing family ID, individual ID and phenotype of the discovery dataset as columns, without heading
#' @param covar_file Name (with file extension) of the covariate file containing family ID, individual ID, standardized covariate, square of standardized covariate, and/or confounders of the discovery dataset as columns, without heading
#' @param n_confounders Number of confounding variables in the discovery dataset
#' @param thread Number of threads used
#' 
#' @keywords gwas
#' 
#' @export 
#' 
#' @importFrom 
#' 
#' @return This function will perform GWAS and output
#' \item{B_trd.sum} GWAS summary statistics file with additive SNP effects
#' 
#' @example 
#' GWAS_binary(DummyData, Bphe_discovery.txt, Bcov_discovery.txt, 14, 20)
#' GWAS_binary(DummyData, Bphe_discovery.txt, Bcov_discovery.txt, 0, 20)


GWAS_binary <- function(b_file, pheno_file, covar_file, n_confounders, thread){
  
  parameters <- c(1, (1:n_confounders)+3)
  param_vec <- paste0(parameters, collapse = ", ")

  system(paste0("./plink2 --bfile ", b_file, 
                " --glm --pheno ", 
                pheno_file, 
                " --covar ", covar_file, 
                " --parameters ", param_vec, 
                " --allow-no-sex --threads ", 
                thread,
                " --out B_gwas"))
  
  plink_output <- read.table("B_gwas.PHENO1.glm.logistic.hybrid", header = FALSE)
  filtered_output <- plink_output[(plink_output$V8=="ADD"),]
  filtered_output$V10 = log(filtered_output$V10)
  
  sink("B_trd.sum")
  write.table(filtered_output, sep = " ", row.names = FALSE, quote = FALSE)
  sink()
  
}
