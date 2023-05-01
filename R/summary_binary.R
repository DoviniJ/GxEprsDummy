#' summary_binary function
#' 
#' This function uses plink2 and outputs the summary of model in the target dataset, using pre-generated GWAS and GWEIS summary statistics files named B_trd.sum, B_add.sum and B_gxe.sum
#' 
#' @param n_confounders Number of confounding variables in the target dataset
#' 
#' @keywords risk scores, risk values, individual risk, disease risk
#' 
#' @export 
#' 
#' @importFrom 
#' 
#' @return This function will output
#' \item{Bsummary.txt} the summary of the fitted model
#' 
#' @example 
#' summary_binary(0)
#' summary_binary(14)


summary_binary <- function(n_confounders){
  
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
  
  cov=scale(dat$V3[m1])
  ps0=scale(prs0$V5)
  ps1=scale(prs1$V5)
  ps2=scale(prs2$V5)
  xv0=scale(prs0$V5*cov)
  xv1=scale(prs1$V5*cov)
  xv2=scale(prs2$V5*cov)
  
  cov2=scale(cov^2)
  
  if(n_confounders == 0){
    df_new <- as.data.frame(cbind(out, cov, cov2, ps1, ps2, xv2))
    colnames(df_new)[2] <- "E"
    colnames(df_new)[3] <- "E squared"
    colnames(df_new)[4] <- "PRS_add"
    colnames(df_new)[5] <- "PRS_gxe"
    colnames(df_new)[6] <- "PRS_gxe x E"
    m = glm(out ~., data = df_new, family = binomial(link = logit))
    sink("Bsummary.txt")
    print(summary(m))
    sink()
  }else{
    conf_var <- matrix(ncol = n_confounders, nrow = nrow(dat))
    for (k in 1:n_confounders) {
      conf_var[, k] <- as.numeric(dat[, k+4])
    }
    conf_var <- conf_var[m1,]
    df_new <- as.data.frame(cbind(out, cov, cov2, ps1, ps2, xv2, conf_var))
    colnames(df_new)[2] <- "E"
    colnames(df_new)[3] <- "E squared"
    colnames(df_new)[4] <- "PRS_add"
    colnames(df_new)[5] <- "PRS_gxe"
    colnames(df_new)[6] <- "PRS_gxe x E"
    m = glm(out ~., data = df_new, family = binomial(link = logit))
    sink("Bsummary.txt")
    print(summary(m))
    sink()
  }
  
}
