#' PRS_binary function
#' This function uses plink2 and outputs PRSs of each individual in the target dataset, using pre-generated GWAS and GWEIS summary statistics files named Q_out.trd.sum, Q_out.add.sum or Q_out.gxe.sum
#' @param plink_path Path to the PLINK executable application
#' @param b_file Prefix of the binary files, where all .fam, .bed and .bim files have a common prefix
#' @param summary_input Name of the summary statistics file specified by the user
#' @param summary_output Name of the PRS file generated using provided summary statistics file specified by the user
#' @keywords prs, profile scores
#' @export 
#' @importFrom stats D cor dnorm
#' @return This function will output
#' \item{Q_trd.sscore} PRSs for each target individual using GWAS additive effects
#' @examples \dontrun{ 
#' x <- PRS_quantitative(plink_path, DummyData)
#' head(x[[1]])
#' head(x[[2]])
#' head(x[[3]])
#' }
PRS_quantitative <- function(plink_path, b_file, summary_input = "Q_out.trd.sum", summary_output = "Q_trd"){               
  runPLINK <- function(PLINKoptions = "") system(paste(plink_path, PLINKoptions))
    runPLINK(paste0(" --bfile ", b_file, 
                    " --score ", noquote(summary_input), " 3 6 9 header --out ", summary_output))
    summary_output_name <- paste0(noquote(summary_output), ".sscore")
    summary_output <- as.character(summary_output_name)
    out = read.table(summary_output, header = F)
  return(out)
}