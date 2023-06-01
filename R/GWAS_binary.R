#' GWAS_binary function
#' This function performs GWAS using plink2 and outputs the GWAS summary statistics file with additive SNP effects named B_trd.sum
#' @param plink_path Path to the PLINK executable application
#' @param b_file Prefix of the binary files, where all .fam, .bed and .bim files have a common prefix
#' @param Bphe_discovery Name (with file extension) of the phenotype file containing family ID, individual ID and phenotype of the discovery dataset as columns, without heading
#' @param Bcov_discovery Name (with file extension) of the covariate file containing family ID, individual ID, standardized covariate, square of standardized covariate, and/or confounders of the discovery dataset as columns, without heading
#' @param thread Number of threads used
#' @param summary_output Name of the SNP effects of the GWAS summary statistics file specified by the user
#' @keywords gwas
#' @export 
#' @importFrom stats D cor dnorm
#' @return This function will perform GWAS and output
#' \item{B_trd.sum} GWAS summary statistics file with additive SNP effects
#' @example x <- GWAS_binary(plink_path, DummyData, Bphe_discovery, Bcov_discovery)
#' @example head(x)
GWAS_binary <- function(plink_path, b_file, Bphe_discovery, Bcov_discovery, thread = 20, summary_output = "B_trd.sum"){  
  cov_file <- read.table(Bcov_discovery)
  n_confounders = ncol(cov_file) - 4
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
                  Bphe_discovery, 
                  " --covar ", Bcov_discovery, 
                  " --parameters ", param_vec, 
                  " --allow-no-sex --threads ", 
                  thread,
                  " --out B_gwas"))
  plink_output <- read.table("B_gwas.PHENO1.glm.logistic.hybrid", header = FALSE)
  filtered_output <- plink_output[(plink_output$V8=="ADD"),]
  filtered_output$V10 = log(filtered_output$V10)
  sink(summary_output)
  write.table(filtered_output, sep = " ", row.names = FALSE, quote = FALSE)
  sink()
  return(filtered_output)
}
