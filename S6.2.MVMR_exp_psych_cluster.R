rm(list=ls())

library(TwoSampleMR)
library(data.table)
library(ieugwasr)
library(MVMR)

setwd("")
set.seed(4477)

# load scz GWAS
scz <- data.frame(fread("scz_EUR.gwaslab.tsv.gz"))
scz$Phenotype <- "scz"

# load mdd GWAS
mdd <- data.frame(fread("mdd_EUR_new.gwaslab.tsv.gz"))
mdd$Phenotype <- "mdd"

# load phenotypic correlation
pheno_cor <- read.csv("./mdd_new/pheno_cor.csv", header = T)
row.names(pheno_cor) <- pheno_cor$X
pheno_cor <- pheno_cor[-1]

############# psych as outcome

pheno <- c("bmi", "smkinit", "cigday")
covar <- c("bmi", "smkinit", "cigday", "canud", "ea", "income")

mv_res_exp_scz <- list()
mv_res_exp_mdd <- list()

for (x in 1:length(pheno)){
  
  mv_res_exp_scz[[x]] <- list()
  mv_res_exp_mdd[[x]] <- list()
  
  covar_new <- covar[-x]
  
  for (i in 1:length(covar_new)){
    
    print(c(pheno[x], covar_new[i]))
    
    if (covar_new[i]=="canud"){
      mv_exposure_dat <- mv_extract_exposures_local(filenames_exposure = c(paste0("", pheno[x], "_EUR.gwaslab.tsv.gz"), 
                                                                           paste0("",covar_new[i],"_EUR.gwaslab.tsv.gz")),
                                                    sep = "\t",
                                                    snp_col = "rsID",
                                                    beta_col = "BETA",
                                                    se_col = "SE",
                                                    eaf_col = c("EAF","rsid"), # canud does not have eaf, so we load rsID as a placeholder
                                                    effect_allele_col = "EA",
                                                    other_allele_col = "NEA",
                                                    pval_col = "P",
                                                    pval_threshold = 5e-08,
                                                    plink_bin = "plink",
                                                    bfile = "1kg/EUR")
    } else {
      mv_exposure_dat <- mv_extract_exposures_local(filenames_exposure = c(paste0("", pheno[x], "_EUR.gwaslab.tsv.gz"), 
                                                                           paste0("",covar_new[i],"_EUR.gwaslab.tsv.gz")),
                                                    sep = "\t",
                                                    snp_col = "rsID",
                                                    beta_col = "BETA",
                                                    se_col = "SE",
                                                    eaf_col = "EAF",
                                                    effect_allele_col = "EA",
                                                    other_allele_col = "NEA",
                                                    pval_col = "P",
                                                    pval_threshold = 5e-08,
                                                    plink_bin = "plink",
                                                    bfile = "1kg/EUR")
    }
    
    mv_exposure_dat$exposure <- ifelse(mv_exposure_dat$exposure=="exposure1", pheno[x], covar_new[i])
    
    mv_exposure_dat$id.exposure <- ifelse(mv_exposure_dat$exposure==pheno[x], 1, 2)
    
    ###### for scz
    mv_outcome_dat <- scz[scz$rsID %in% mv_exposure_dat$SNP,]
    
    mv_outcome_dat <- format_data(mv_outcome_dat, 
                                  type = "outcome",
                                  phenotype_col = "Phenotype",
                                  snp_col = "rsID",
                                  beta_col = "BETA",
                                  se_col = "SE",
                                  eaf_col = "EAF",
                                  effect_allele_col = "EA",
                                  other_allele_col = "NEA",
                                  pval_col = "P")
    
    # mv harmonise
    mvdat <- mv_harmonise_data(mv_exposure_dat, mv_outcome_dat)
    
    # get number of instruments
    nsnp <- as.numeric(lapply(data.frame(mvdat$exposure_pval), function(x) sum(x<5e-08)))

    # run MVMR using MVMR package
    bx <- as.matrix(mvdat$exposure_beta)
    bxse <- as.matrix(mvdat$exposure_se)
    F.data <- MendelianRandomization::mr_mvinput(bx = bx,
                                                 bxse = bxse,
                                                 by = mvdat$outcome_beta,
                                                 byse = mvdat$outcome_se,
                                                 snps= row.names(bx))
    F.data <- mrmvinput_to_mvmr_format(F.data)

    # Estimating pairwise covariances
    mvmrcovmatrix <- as.matrix(pheno_cor[row.names(pheno_cor) %in% c(pheno[x],covar_new[i]),
                                         names(pheno_cor) %in% c(pheno[x],covar_new[i])])
    Xcovmat <- phenocov_mvmr(mvmrcovmatrix,F.data[,c("sebetaX1","sebetaX2")])
    
    # test weak instrument
    sres <- strength_mvmr(r_input = F.data, gencov = Xcovmat)
    
    # test for horizontal pleiotropy
    pres <- pleiotropy_mvmr(r_input = F.data, gencov = Xcovmat)
    
    # estimate causal effects
    mvmr_res <- data.frame(ivw_mvmr(r_input = F.data, gencov = Xcovmat))
    
    # combine
    mvmr_res$nsnp <- nsnp 
    mvmr_res$F <- as.numeric(sres[1,])
    mvmr_res$Qstat <- pres$Qstat
    mvmr_res$Qpval <- pres$Qpval
    
    # Robust causal effect estimation
    mvmr_res_1 <- qhet_mvmr(F.data, mvmrcovmatrix, CI = F, iterations = 100)
    mvmr_res$Estimates_robust <- as.numeric(unlist(mvmr_res_1))
    
    # add exp and out info
    mvmr_res$exposure <- c(pheno[x], covar_new[i])
    mvmr_res$outcome <- "scz"
    
    # add to list
    mv_res_exp_scz[[x]][[i]] <- mvmr_res
    
    
    ###### for mdd
    mv_outcome_dat <- mdd[mdd$rsID %in% mv_exposure_dat$SNP,]
    
    mv_outcome_dat <- format_data(mv_outcome_dat, 
                                  type = "outcome",
                                  phenotype_col = "Phenotype",
                                  snp_col = "rsID",
                                  beta_col = "BETA",
                                  se_col = "SE",
                                  eaf_col = "EAF",
                                  effect_allele_col = "EA",
                                  other_allele_col = "NEA",
                                  pval_col = "P")
    
    # mv harmonise
    mvdat <- mv_harmonise_data(mv_exposure_dat, mv_outcome_dat)
    
    # get number of instruments
    nsnp <- as.numeric(lapply(data.frame(mvdat$exposure_pval), function(x) sum(x<5e-08)))
    
    # run MVMR using MVMR package
    bx <- as.matrix(mvdat$exposure_beta)
    bxse <- as.matrix(mvdat$exposure_se)
    F.data <- MendelianRandomization::mr_mvinput(bx = bx,
                                                 bxse = bxse,
                                                 by = mvdat$outcome_beta,
                                                 byse = mvdat$outcome_se,
                                                 snps= row.names(bx))
    F.data <- mrmvinput_to_mvmr_format(F.data)
    
    # Estimating pairwise covariances
    mvmrcovmatrix <- as.matrix(pheno_cor[row.names(pheno_cor) %in% c(pheno[x],covar_new[i]),
                                         names(pheno_cor) %in% c(pheno[x],covar_new[i])])
    Xcovmat <- phenocov_mvmr(mvmrcovmatrix,F.data[,c("sebetaX1","sebetaX2")])
    
    # test weak instrument
    sres <- strength_mvmr(r_input = F.data, gencov = Xcovmat)
    
    # test for horizontal pleiotropy
    pres <- pleiotropy_mvmr(r_input = F.data, gencov = Xcovmat)
    
    # estimate causal effects
    mvmr_res <- data.frame(ivw_mvmr(r_input = F.data, gencov = Xcovmat))
    
    # combine
    mvmr_res$nsnp <- nsnp 
    mvmr_res$F <- as.numeric(sres[1,])
    mvmr_res$Qstat <- pres$Qstat
    mvmr_res$Qpval <- pres$Qpval
    
    # Robust causal effect estimation
    mvmr_res_1 <- qhet_mvmr(F.data, mvmrcovmatrix, CI = F, iterations = 100)
    mvmr_res$Estimates_robust <- as.numeric(unlist(mvmr_res_1))
    
    # add exp and out info
    mvmr_res$exposure <- c(pheno[x], covar_new[i])
    mvmr_res$outcome <- "mdd"
    
    # add to list
    mv_res_exp_mdd[[x]][[i]] <- mvmr_res
    
    }
  
}

save(mv_res_exp_scz, file = "mv_res_exp_scz.RData")
save(mv_res_exp_mdd, file = "mv_res_exp_mdd.RData")
