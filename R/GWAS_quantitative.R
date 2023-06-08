#' GWAS_quantitative function
#' This function performs GWAS using plink2 and outputs the GWAS summary statistics file with additive SNP effects named Q_out.trd.sum
#' @param plink_path Path to the PLINK executable application
#' @param b_file Prefix of the binary files, where all .fam, .bed and .bim files have a common prefix
#' @param Qphe_discovery Name (with file extension) of the phenotype file containing family ID, individual ID and phenotype of the discovery dataset as columns, without heading
#' @param Qcov_discovery Name (with file extension) of the covariate file containing family ID, individual ID, standardized covariate, square of standardized covariate, and/or confounders of the discovery dataset as columns, without heading
#' @param thread Number of threads used
#' @param summary_output Name of the SNP effects of the GWAS summary statistics file specified by the user
#' @keywords gwas
#' @export 
#' @importFrom stats D cor dnorm
#' @return This function will perform GWAS and output
#' \item{Q_out.trd.sum} GWAS summary statistics file with additive SNP effects
#' @example x <- GWAS_quantitative(plink_path, DummyData, Qphe_discovery, Qcov_discovery)
#' @example head(x)
GWAS_quantitative <- function(plink_path, b_file, Qphe_discovery, Qcov_discovery, thread = 20, summary_output = "Q_out.trd.sum"){  
  cov_file <- read.table(Qcov_discovery)
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
                  " --linear interaction --pheno ", 
                  Qphe_discovery, 
                  " --covar ", Qcov_discovery, 
                  " --parameters ", param_vec, 
                  " --allow-no-sex --threads ", 
                  thread,
                  " --out Q_gwas"))
  plink_output <- read.table("Q_gwas.PHENO1.glm.linear", header = FALSE)
  filtered_output <- plink_output[(plink_output$V7=="ADD"),]
  sink(summary_output)
  write.table(filtered_output, sep = " ", row.names = FALSE, quote = FALSE)
  sink()
  return(filtered_output)
}
