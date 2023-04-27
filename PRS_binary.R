#' PRS_binary function
#' 
#' This function uses plink2 and outputs PRSs of each individual in the target dataset, using pre-generated GWAS and GWEIS summary statistics files named B_trd.sum, B_add.sum and B_gxe.sum
#' 
#' @param b_file Prefix of the binary files, where all .fam, .bed and .bim files have a common prefix
#' 
#' @keywords prs, profile scores
#' 
#' @export 
#' 
#' @importFrom 
#' 
#' @return This function will output
#' \item{B_trd.sscore} PRSs for each target individual using GWAS additive effects
#' \item{B_add.sscore} PRSs for each target individual using GWEIS additive effects
#' \item{B_gxe.sscore} PRSs for each target individual using GWEIS interaction effects
#' 
#' @example 
#' PRS_binary(DummyData)


PRS_binary <- function(b_file){

  system(paste0("plink2 --bfile ", b_file, 
                " --score B_trd.sum 3 6 10 header --out B_trd"))
  system(paste0("plink2 --bfile ", b_file, 
                " --score B_add.sum 3 6 10 header --out B_add"))
  system(paste0("plink2 --bfile ", b_file, 
                " --score B_gxe.sum 3 6 10 header --out B_gxe"))
  
}