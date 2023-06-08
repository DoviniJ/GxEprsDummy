#' GWEIS_binary function
#' This function performs GWEIS using plink2 and outputs the GWEIS summary statistics files with additive SNP effects named B_out.add.sum and interaction SNP effects named B_out.gxe.sum
#' @param plink_path Path to the PLINK executable application
#' @param b_file Prefix of the binary files, where all .fam, .bed and .bim files have a common prefix
#' @param Bphe_discovery Phenotype file containing family ID, individual ID and phenotype of the discovery dataset as columns, without heading
#' @param Bcov_discovery Covariate file containing family ID, individual ID, standardized covariate, square of standardized covariate, and/or confounders of the discovery dataset as columns, without heading
#' @param thread Number of threads used
#' @param summary_output Name (prefix) of the additive or interaction SNP effects of the GWEIS summary statistics file specified by the user
#' @keywords gwies, interaction, gxe
#' @export 
#' @importFrom stats D cor dnorm
#' @return This function will perform GWEIS and output
#' \item{B_out.add.sum} GWEIS summary statistics file with additive SNP effects
#' \item{B_out.gxe.sum} GWEIS summary statistics file with interaction SNP effects
#' @example x <- GWEIS_binary(plink_path, DummyData, Bphe_discovery, Bcov_discovery)
#' @example head(x[[1]])
#' @example head(x[[2]])
GWEIS_binary <- function(plink_path, b_file, Bphe_discovery, Bcov_discovery, thread = 20, summary_output = "B_out"){
  cov_file <- read.table(Bcov_discovery)
  n_confounders = ncol(cov_file) - 4
  if(n_confounders > 0){
    parameters <- c(1, 2, 3, (1:(n_confounders+1))+3)
  }
  else{
    parameters <- c(1, 2, 3, 4)
  }
  param_vec <- paste0(parameters, collapse = ", ")
  runPLINK <- function(PLINKoptions = "") system(paste(plink_path, PLINKoptions))
  runPLINK(paste0(" --bfile ", b_file,
                " --glm interaction --pheno ", 
                Bphe_discovery, 
                " --covar ", Bcov_discovery, 
                " --parameters ", param_vec, 
                " --allow-no-sex --threads ", 
                thread,
                " --out B_gweis"))
  plink_output <- read.table("B_gweis.PHENO1.glm.logistic.hybrid", header = FALSE)
  filtered_output <- as.data.frame(plink_output[(plink_output$V8=="ADD"),])
  filtered_output$V10 = log(filtered_output$V10)
  filtered_output2 <- plink_output[(plink_output$V8=="ADDxCOVAR1"),]
  filtered_output2$V10 <- log(filtered_output2$V10)
  summary_output1 <- paste0(noquote(summary_output), ".add.sum")
  summary_output2 <- paste0(noquote(summary_output), ".gxe.sum")
  sink(summary_output1)
  write.table(filtered_output, sep = " ", row.names = FALSE, quote = FALSE)
  sink()
  sink(summary_output2)
  write.table(filtered_output2, sep = " ", row.names = FALSE, quote = FALSE)
  sink()
 out <- list(filtered_output, filtered_output2)
 return(out)
 }

