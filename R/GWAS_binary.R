#' GWAS_binary function
#' 
#' This function performs GWAS using plink2 and outputs the GWAS summary statistics file with additive SNP effects named B_trd.sum
#' 
#' @param plink_path Path to the PLINK executable application
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
#' 
#' @return This function will perform GWAS and output
#' \item{B_trd.sum} GWAS summary statistics file with additive SNP effects
#' 
#' @example GWAS_binary(plink_path, DummyData, Bphe_discovery, Bcov_discovery, 14, 20)


GWAS_binary <- function(plink_path, b_file, pheno_file, covar_file, n_confounders, thread){
  if(n_confounders > 0){
    parameters <- c(1, (1:n_confounders)+3)
  }
  else{
    parameters <- 1
  }
  param_vec <- paste0(parameters, collapse = ", ")
  runPLINK <- function(PLINKoptions = "") system(paste(plink_path, PLINKoptions))
  runPLINK(paste0(" --bfile ", b_file, 
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

