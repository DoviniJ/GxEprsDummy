#' summary_regular_binary function
#' This function uses plink2 and outputs the summary of regular model in the target dataset, using pre-generated GWAS and GWEIS summary statistics files named B_trd.sum, B_add.sum and B_gxe.sum
#' @param Bphe_target Phenotype file containing family ID, individual ID and phenotype of the target dataset as columns, without heading
#' @param Bcov_target Covariate file containing family ID, individual ID, standardized covariate, square of standardized covariate, and/or confounders of the target dataset as columns, without heading
#' @param input_score1 The .sscore file generated using additive SNP effects of GWAS summary statistics
#' @param input_score2 The .sscore file generated using additive SNP effects of GWEIS summary statistics
#' @param input_score3 The .sscore file generated using interaction SNP effects of GWEIS summary statistics
#' @param Model Specify the model number (1: y = PRS_trd + E + PRS_trd x E + confounders, 2: y = PRS_add + E + PRS_add x E + confounders, 3: y = PRS_add + E + PRS_gxe x E + confounders, 4: y = PRS_add + E + PRS_gxe + PRS_gxe x E + confounders, 5: y = PRS_add + E + E^2 + PRS_gxe + PRS_gxe x E + confounders, where y is the outcome variable, E is the covariate of interest, PRS_trd and PRS_add are the polygenic risk scores computed using additive SNP effects of GWAS and GWEIS summary statistics respectively, and PRS_gxe is the polygenic risk scores computed using GxE interaction SNP effects of GWEIS summary statistics.)
#' @keywords regression summary
#' @export 
#' @importFrom stats D cor dnorm
#' @return This function will output
#' \item{Bsummary.txt} the summary of the fitted model
#' \item{Individual_risk_values.txt} the estimated risk values of individuals in the target sample
#' @example x <- summary_regular_binary(Bphe_target, Bcov_target, Model = 5)
#' @example x[[1]][[1]]
#' @example x[[1]][[2]]
#' @example x[[1]][[3]]
#' @example x[[1]][[4]]
#' @example x[[1]][[5]]
#' @example x[[1]][[6]]
#' @example x[[1]][[7]]
#' @example x[[1]][[8]]
#' @example x[[1]][[9]]
#' @example x[[1]][[10]]
#' @example x[[1]][[11]]
#' @example x[[1]][[12]]
#' @example x[[1]][[13]]
#' @example x[[1]][[14]]
#' @example x[[1]][[15]]
#' @example x[[1]][[16]]
#' @example x[[1]][[17]]
#' @example x[[1]][[18]]
#' @example head(x[[2]])
summary_regular_binary <- function(Bphe_target, Bcov_target, input_score1 = "B_trd.sscore", input_score2 = "B_add.sscore", input_score3 = "B_gxe.sscore", Model){
  cov_file <- read.table(Bcov_target)
  n_confounders = ncol(cov_file) - 4
  fam=read.table(Bphe_target, header=F) 
  colnames(fam) <- c("FID", "IID", "PHENOTYPE")
  dat=read.table(Bcov_target ,header=F)
  colnames(dat)[1] <- "FID"
  colnames(dat)[2] <- "IID"
  if(file.exists(input_score1)){
    prs0_all=read.table(input_score1)
    colnames(prs0_all)[1] <- "FID"
    colnames(prs0_all)[2] <- "IID"
    prs0=merge(fam, prs0_all, by = "FID")
    m1 <- match(prs0$IID.x, dat$IID)
    ps0=scale(prs0$V5)
    out = fam$PHENOTYPE[m1]
    cov=scale(dat$V3[m1])
    xv0=scale(prs0$V5*cov)
    cov2=scale(cov^2)
  }
  if(file.exists(input_score2)){
    prs1_all=read.table(input_score2)
    colnames(prs1_all)[1] <- "FID"
    colnames(prs1_all)[2] <- "IID"
    prs1=merge(fam, prs1_all, by = "FID")
    m1 <- match(prs1$IID.x, dat$IID)
    ps1=scale(prs1$V5)
    out = fam$PHENOTYPE[m1]
    cov=scale(dat$V3[m1])
    xv1=scale(prs1$V5*cov)
    cov2=scale(cov^2)
  }
  if(file.exists(input_score3)){
    prs2_all=read.table(input_score3)
    colnames(prs2_all)[1] <- "FID"
    colnames(prs2_all)[2] <- "IID"
    prs2=merge(fam, prs2_all, by = "FID")
    m1 <- match(prs2$IID.x, dat$IID)
    ps2=scale(prs2$V5)
    out = fam$PHENOTYPE[m1]
    cov=scale(dat$V3[m1])
    xv2=scale(prs2$V5*cov)
    cov2=scale(cov^2)
  }
  if(Model == 1){
    if(n_confounders == 0){
      df_new <- as.data.frame(cbind(out, cov, ps0, xv0))
      colnames(df_new)[2] <- "E"
      colnames(df_new)[3] <- "PRS_trd"
      colnames(df_new)[4] <- "PRS_trd x E"
      m = glm(out ~., data = df_new, family = binomial(link = logit))
      sink("Bsummary.txt")
      print(summary(m))
      sink()
      m_fit <- fitted.values(m)
      sink("Individual_risk_values.txt")
      write.table(cbind(fam$FID[m1], fam$IID[m1], m_fit), row.names = F, col.names = F, quote = F)
      sink()
    }else{
      conf_var <- matrix(ncol = n_confounders, nrow = nrow(dat))
      for (k in 1:n_confounders) {
        conf_var[, k] <- as.numeric(dat[, k+4])
      }
      conf_var <- conf_var[m1,]
      df_new <- as.data.frame(cbind(out, cov, ps0, xv0, conf_var))
      colnames(df_new)[2] <- "E"
      colnames(df_new)[3] <- "PRS_trd"
      colnames(df_new)[4] <- "PRS_trd x E"
      m = glm(out ~., data = df_new, family = binomial(link = logit))
      sink("Bsummary.txt")
      print(summary(m))
      sink()
      m_fit <- fitted.values(m)
      sink("Individual_risk_values.txt")
      write.table(cbind(fam$FID[m1], fam$IID[m1], m_fit), row.names = F, col.names = F, quote = F)
      sink()
    }
    s <- summary(m)
    out1 <- list(s, s$call, s$terms, s$family, s$deviance, s$aic, s$contrasts, 
                 s$df.residual, s$null.deviance, s$df.null, s$iter, s$deviance.resid, 
                 s$coefficients, s$aliesed, s$dispersion, s$df, s$cov.unscaled, s$cov.scaled)
    out2 <- cbind(fam$FID[m1], fam$IID[m1], m_fit)
    out_all <- list(out1, out2)
  }
  if(Model == 2){
    if(n_confounders == 0){
      df_new <- as.data.frame(cbind(out, cov, ps1, xv1))
      colnames(df_new)[2] <- "E"
      colnames(df_new)[3] <- "PRS_add"
      colnames(df_new)[4] <- "PRS_add x E"
      m = glm(out ~., data = df_new, family = binomial(link = logit))
      sink("Bsummary.txt")
      print(summary(m))
      sink()
      m_fit <- fitted.values(m)
      sink("Individual_risk_values.txt")
      write.table(cbind(fam$FID[m1], fam$IID[m1], m_fit), row.names = F, col.names = F, quote = F)
      sink()
    }else{
      conf_var <- matrix(ncol = n_confounders, nrow = nrow(dat))
      for (k in 1:n_confounders) {
        conf_var[, k] <- as.numeric(dat[, k+4])
      }
      conf_var <- conf_var[m1,]
      df_new <- as.data.frame(cbind(out, cov, ps1, xv1, conf_var))
      colnames(df_new)[2] <- "E"
      colnames(df_new)[3] <- "PRS_add"
      colnames(df_new)[4] <- "PRS_add x E"
      m = glm(out ~., data = df_new, family = binomial(link = logit))
      sink("Bsummary.txt")
      print(summary(m))
      sink()
      m_fit <- fitted.values(m)
      sink("Individual_risk_values.txt")
      write.table(cbind(fam$FID[m1], fam$IID[m1], m_fit), row.names = F, col.names = F, quote = F)
      sink()
    }
    s <- summary(m)
    out1 <- list(s, s$call, s$terms, s$family, s$deviance, s$aic, s$contrasts, 
                 s$df.residual, s$null.deviance, s$df.null, s$iter, s$deviance.resid, 
                 s$coefficients, s$aliesed, s$dispersion, s$df, s$cov.unscaled, s$cov.scaled)
    out2 <- cbind(fam$FID[m1], fam$IID[m1], m_fit)
    out_all <- list(out1, out2)
  }
  if(Model == 3){
    if(n_confounders == 0){
      df_new <- as.data.frame(cbind(out, cov, ps1, xv2))
      colnames(df_new)[2] <- "E"
      colnames(df_new)[3] <- "PRS_add"
      colnames(df_new)[4] <- "PRS_gxe x E"
      m = glm(out ~., data = df_new, family = binomial(link = logit))
      sink("Bsummary.txt")
      print(summary(m))
      sink()
      m_fit <- fitted.values(m)
      sink("Individual_risk_values.txt")
      write.table(cbind(fam$FID[m1], fam$IID[m1], m_fit), row.names = F, col.names = F, quote = F)
      sink()
    }else{
      conf_var <- matrix(ncol = n_confounders, nrow = nrow(dat))
      for (k in 1:n_confounders) {
        conf_var[, k] <- as.numeric(dat[, k+4])
      }
      conf_var <- conf_var[m1,]
      df_new <- as.data.frame(cbind(out, cov, ps1, xv2, conf_var))
      colnames(df_new)[2] <- "E"
      colnames(df_new)[3] <- "PRS_add"
      colnames(df_new)[4] <- "PRS_gxe x E"
      m = glm(out ~., data = df_new, family = binomial(link = logit))
      sink("Bsummary.txt")
      print(summary(m))
      sink()
      m_fit <- fitted.values(m)
      sink("Individual_risk_values.txt")
      write.table(cbind(fam$FID[m1], fam$IID[m1], m_fit), row.names = F, col.names = F, quote = F)
      sink()
    }
    s <- summary(m)
    out1 <- list(s, s$call, s$terms, s$family, s$deviance, s$aic, s$contrasts, 
                 s$df.residual, s$null.deviance, s$df.null, s$iter, s$deviance.resid, 
                 s$coefficients, s$aliesed, s$dispersion, s$df, s$cov.unscaled, s$cov.scaled)
    out2 <- cbind(fam$FID[m1], fam$IID[m1], m_fit)
    out_all <- list(out1, out2)
  }
  if(Model == 4){
    if(n_confounders == 0){
      df_new <- as.data.frame(cbind(out, cov, ps1, ps2, xv2))
      colnames(df_new)[2] <- "E"
      colnames(df_new)[3] <- "PRS_add"
      colnames(df_new)[4] <- "PRS_gxe"
      colnames(df_new)[5] <- "PRS_gxe x E"
      m = glm(out ~., data = df_new, family = binomial(link = logit))
      sink("Bsummary.txt")
      print(summary(m))
      sink()
      m_fit <- fitted.values(m)
      sink("Individual_risk_values.txt")
      write.table(cbind(fam$FID[m1], fam$IID[m1], m_fit), row.names = F, col.names = F, quote = F)
      sink()
    }else{
      conf_var <- matrix(ncol = n_confounders, nrow = nrow(dat))
      for (k in 1:n_confounders) {
        conf_var[, k] <- as.numeric(dat[, k+4])
      }
      conf_var <- conf_var[m1,]
      df_new <- as.data.frame(cbind(out, cov, ps1, ps2, xv2, conf_var))
      colnames(df_new)[2] <- "E"
      colnames(df_new)[3] <- "PRS_add"
      colnames(df_new)[4] <- "PRS_gxe"
      colnames(df_new)[5] <- "PRS_gxe x E"
      m = glm(out ~., data = df_new, family = binomial(link = logit))
      sink("Bsummary.txt")
      print(summary(m))
      sink()
      m_fit <- fitted.values(m)
      sink("Individual_risk_values.txt")
      write.table(cbind(fam$FID[m1], fam$IID[m1], m_fit), row.names = F, col.names = F, quote = F)
      sink()
    }
    s <- summary(m)
    out1 <- list(s, s$call, s$terms, s$family, s$deviance, s$aic, s$contrasts, 
                 s$df.residual, s$null.deviance, s$df.null, s$iter, s$deviance.resid, 
                 s$coefficients, s$aliesed, s$dispersion, s$df, s$cov.unscaled, s$cov.scaled)
    out2 <- cbind(fam$FID[m1], fam$IID[m1], m_fit)
    out_all <- list(out1, out2)
  }
  if(Model == 5){
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
      m_fit <- fitted.values(m)
      sink("Individual_risk_values.txt")
      write.table(cbind(fam$FID[m1], fam$IID[m1], m_fit), row.names = F, col.names = F, quote = F)
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
      m_fit <- fitted.values(m)
      sink("Individual_risk_values.txt")
      write.table(cbind(fam$FID[m1], fam$IID[m1], m_fit), row.names = F, col.names = F, quote = F)
      sink()
    }
    s <- summary(m)
    out1 <- list(s, s$call, s$terms, s$family, s$deviance, s$aic, s$contrasts, 
                 s$df.residual, s$null.deviance, s$df.null, s$iter, s$deviance.resid, 
                 s$coefficients, s$aliesed, s$dispersion, s$df, s$cov.unscaled, s$cov.scaled)
    out2 <- cbind(fam$FID[m1], fam$IID[m1], m_fit)
    out_all <- list(out1, out2)
  }
  return(out_all)
}
