#' PRS_binary function
#' This function uses plink2 and outputs PRSs of each individual in the target dataset, using pre-generated GWAS and GWEIS summary statistics files named B_trd.sum, B_add.sum and B_gxe.sum
#' @param plink_path Path to the PLINK executable application
#' @param b_file Prefix of the binary files, where all .fam, .bed and .bim files have a common prefix
#' @keywords prs, profile scores
#' @export 
#' @importFrom stats D cor dnorm
#' @return This function will output
#' \item{B_trd.sscore} PRSs for each target individual using GWAS additive effects
#' \item{B_add.sscore} PRSs for each target individual using GWEIS additive effects
#' \item{B_gxe.sscore} PRSs for each target individual using GWEIS interaction effects
#' @example x <- PRS_binary(plink_path, DummyData)
#' @example head(x[[1]])
#' @example head(x[[2]])
#' @example head(x[[3]])
PRS_binary <- function(plink_path, b_file, summary_input1 = "B_trd.sum", summary_input2 = "B_add.sum", summary_input3 = "B_gxe.sum", summary_output1 = NULL, summary_output2 = NULL, summary_output3 = NULL){
  runPLINK <- function(PLINKoptions = "") system(paste(plink_path, PLINKoptions))
  if(missing(summary_output1)){
    runPLINK(paste0(" --bfile ", b_file, 
                    " --score ", noquote(summary_input1), " 3 6 10 header --out B_trd"))
    s1 = read.table("B_trd.sscore", header = F)
  }
  else{
    runPLINK(paste0(" --bfile ", b_file, 
                    " --score ", noquote(summary_input1), " 3 6 10 header --out ", summary_output1))
    summary_output1_name <- paste0(noquote(summary_output1), ".sscore")
    summary_output1 <- as.character(summary_output1_name)
    s1 = read.table(summary_output1, header = F)
  }
  if(missing(summary_output2)){
    runPLINK(paste0(" --bfile ", b_file, 
                    " --score ", noquote(summary_input2), " 3 6 10 header --out B_add"))
    s2 = read.table("B_add.sscore", header = F)
  }
  else{
    runPLINK(paste0(" --bfile ", b_file, 
                    " --score ", noquote(summary_input2), " 3 6 10 header --out ", summary_output2))
    summary_output2_name <- paste0(noquote(summary_output2), ".sscore")
    summary_output2 <- as.character(summary_output2_name)
    s2 = read.table(summary_output2, header = F)
  }
  if(missing(summary_output3)){
    runPLINK(paste0(" --bfile ", b_file, 
                    " --score ", noquote(summary_input3), " 3 6 10 header --out B_gxe"))
    s3 = read.table("B_gxe.sscore", header = F)
  }
  else{
    runPLINK(paste0(" --bfile ", b_file, 
                    " --score ", noquote(summary_input3), " 3 6 10 header --out ", summary_output3))
    summary_output3_name <- paste0(noquote(summary_output3), ".sscore")
    summary_output3 <- as.character(summary_output3_name)
    s3 = read.table(summary_output3, header = F)
  }
  out <- list(s1, s2, s3)
  return(out)
}
