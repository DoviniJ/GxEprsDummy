#' summary_permuted_binary function
#' 
#' This function uses plink2 and outputs the p value of permuted model in the target dataset, using pre-generated GWAS and GWEIS summary statistics files named B_trd.sum, B_add.sum and B_gxe.sum
#' 
#' @param Bphe_target Phenotype file containing family ID, individual ID and phenotype of the target dataset as columns, without heading
#' @param Bcov_target Covariate file containing family ID, individual ID, standardized covariate, square of standardized covariate, and/or confounders of the target dataset as columns, without heading
#' @param n_confounders Number of confounding variables in the target dataset
#' 
#' @keywords regression summary
#' 
#' @export 
#' 
#' 
#' @return This function will output
#' \item{B_permuted_p.txt} the p value of the permuted model
#' 
#' @example summary_permuted_binary(Bphe_target, Bcov_target, 14)


summary_permuted_binary <- function(Bphe_target, Bcov_target, n_confounders){
  fam=read.table(Bphe_target, header=F) 
  colnames(fam) <- c("FID", "IID", "PHENOTYPE")
  dat=read.table(Bcov_target, header=F)
  colnames(dat)[1] <- "FID"
  colnames(dat)[2] <- "IID"
  prs0_all=read.table("B_trd.sscore")
  colnames(prs0_all)[1] <- "FID"
  colnames(prs0_all)[2] <- "IID"
  prs0=merge(fam, prs0_all, by = "FID")
  prs1_all=read.table("B_add.sscore", header=F)
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
    pn=1000; pp_gxe_x_E=0
    df_regular_new <- as.data.frame(cbind(out, cov, cov2, ps1, ps2, xv2))
    regular_m = glm(out ~., data = df_regular_new, family = binomial(link = logit))
    regular_p = summary(regular_m)$coefficients[6,4]
    for(i in 1:pn){
      sv2=sample(seq(1, length(out)))
      xv2=scale(prs2$V5[sv2]*cov)
    df_new <- as.data.frame(cbind(out, cov, cov2, ps1, ps2, xv2))
    colnames(df_new)[2] <- "E"
    colnames(df_new)[3] <- "E squared"
    colnames(df_new)[4] <- "PRS_add"
    colnames(df_new)[5] <- "PRS_gxe"
    colnames(df_new)[6] <- "PRS_gxe x E"
    m = glm(out ~., data = df_new, family = binomial(link = logit))
    if (regular_p < summary(m)$coefficients[6,4]) pp_gxe_x_E=pp_gxe_x_E+1
    }
    sink("B_permuted_p.txt")
    print(pp_gxe_x_E)
    sink()}else{
    pn=1000; pp_gxe_x_E=0
    df_regular_new <- as.data.frame(cbind(out, cov, cov2, ps1, ps2, xv2))
    regular_m = glm(out ~., data = df_regular_new, family = binomial(link = logit))
    regular_p = summary(regular_m)$coefficients[6,4]
    for(i in 1:pn){
      sv2=sample(seq(1, length(out)))
      xv2=scale(prs2$V5[sv2]*cov)
    }
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
    if (regular_p < summary(m)$coefficients[6,4]) pp_gxe_x_E=pp_gxe_x_E+1
  }
  sink("B_permuted_p.txt")
  print(pp_gxe_x_E)
  sink()
}

