rm(list=ls())

library(TwoSampleMR)
library(ieugwasr)
library(data.table)

setwd("")
set.seed(4477)

############################### EAS

######################## exp to scz

# load bmi GWAS
bmi_EAS <- data.frame(fread("bmi_EAS.gwaslab.tsv.gz"))

bmi_EAS$Phenotype <- "bmi"

bmi_EAS_clump <- ld_clump(
  dplyr::tibble(rsid=bmi_EAS$rsID, pval=bmi_EAS$P, id=bmi_EAS$Phenotype),
  plink_bin = "plink",
  bfile = "1kg/EAS")

bmi_EAS_clumped <- bmi_EAS[bmi_EAS$rsID %in% bmi_EAS_clump$rsid,]

bmi_EAS_clumped <- format_data(bmi_EAS_clumped,
                               type = "exposure",
                               phenotype_col = "Phenotype",
                               snp_col = "rsID",
                               beta_col = "BETA",
                               se_col = "SE",
                               eaf_col = "EAF",
                               effect_allele_col = "EA",
                               other_allele_col = "NEA",
                               pval_col = "P")

bmi_EAS_exp_data <- bmi_EAS_clumped[bmi_EAS_clumped$pval.exposure < 1e-5,]
head(bmi_EAS_exp_data)


# load smkinit GWAS
smkinit_EAS <- data.frame(fread("smkinit_EAS.gwaslab.tsv.gz"))

smkinit_EAS$Phenotype <- "smkinit"

smkinit_EAS_clump <- ld_clump(
  dplyr::tibble(rsid=smkinit_EAS$rsID, pval=smkinit_EAS$P, id=smkinit_EAS$Phenotype),
  plink_bin = "plink",
  bfile = "1kg/EAS")

smkinit_EAS_clumped <- smkinit_EAS[smkinit_EAS$rsID %in% smkinit_EAS_clump$rsid,]

smkinit_EAS_clumped <- format_data(smkinit_EAS_clumped,
                                   type = "exposure",
                                   phenotype_col = "Phenotype",
                                   snp_col = "rsID",
                                   beta_col = "BETA",
                                   se_col = "SE",
                                   eaf_col = "EAF",
                                   effect_allele_col = "EA",
                                   other_allele_col = "NEA",
                                   pval_col = "P")

smkinit_EAS_exp_data <- smkinit_EAS_clumped[smkinit_EAS_clumped$pval.exposure < 1e-5,]
head(smkinit_EAS_exp_data)


# load cigday GWAS
cigday_EAS <- data.frame(fread("cigday_EAS.gwaslab.tsv.gz"))

cigday_EAS$Phenotype <- "cigday"

cigday_EAS_clump <- ld_clump(
  dplyr::tibble(rsid=cigday_EAS$rsID, pval=cigday_EAS$P, id=cigday_EAS$Phenotype),
  plink_bin = "plink",
  bfile = "1kg/EAS")

cigday_EAS_clumped <- cigday_EAS[cigday_EAS$rsID %in% cigday_EAS_clump$rsid,]

cigday_EAS_clumped <- format_data(cigday_EAS_clumped,
                                  type = "exposure",
                                  phenotype_col = "Phenotype",
                                  snp_col = "rsID",
                                  beta_col = "BETA",
                                  se_col = "SE",
                                  eaf_col = "EAF",
                                  effect_allele_col = "EA",
                                  other_allele_col = "NEA",
                                  pval_col = "P")

cigday_EAS_exp_data <- cigday_EAS_clumped[cigday_EAS_clumped$pval.exposure < 1e-5,]
head(cigday_EAS_exp_data)


# combine
EAS_exp_data <- rbind(bmi_EAS_exp_data, smkinit_EAS_exp_data, cigday_EAS_exp_data)


# load scz_EAS GWAS
scz_EAS <- data.frame(fread("scz_EAS.gwaslab.tsv.gz"))
scz_EAS$Phenotype <- "scz"

scz_EAS_out_data <- scz_EAS[scz_EAS$rsID %in% EAS_exp_data$SNP,]

scz_EAS_out_data <- format_data(scz_EAS_out_data, 
                                type = "outcome",
                                phenotype_col = "Phenotype",
                                snp_col = "rsID",
                                beta_col = "BETA",
                                se_col = "SE",
                                eaf_col = "EAF",
                                effect_allele_col = "EA",
                                other_allele_col = "NEA",
                                pval_col = "P")

dat_exp_scz_EAS <- harmonise_data(EAS_exp_data, scz_EAS_out_data, action = 2)

mr_exp_scz_EAS<- mr(dat_exp_scz_EAS, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))

# caluclate F-statistics
mr_exp_scz_EAS$F <- NA

F_bmi_scz_EAS <- bmi_EAS_exp_data[bmi_EAS_exp_data$SNP %in% scz_EAS_out_data$SNP,]
mr_exp_scz_EAS$F[mr_exp_scz_EAS$exposure=="bmi"] <- mean(F_bmi_scz_EAS$beta.exposure^2/F_bmi_scz_EAS$se.exposure^2)

F_smkinit_scz_EAS <- smkinit_EAS_exp_data[smkinit_EAS_exp_data$SNP %in% scz_EAS_out_data$SNP,]
mr_exp_scz_EAS$F[mr_exp_scz_EAS$exposure=="smkinit"]<- mean(F_smkinit_scz_EAS$beta.exposure^2/F_smkinit_scz_EAS$se.exposure^2)

F_cigday_scz_EAS <- cigday_EAS_exp_data[cigday_EAS_exp_data$SNP %in% scz_EAS_out_data$SNP,]
mr_exp_scz_EAS$F[mr_exp_scz_EAS$exposure=="cigday"] <- mean(F_cigday_scz_EAS$beta.exposure^2/F_cigday_scz_EAS$se.exposure^2)


################ scz to out

# scz
scz_EAS_clump <- ld_clump(
  dplyr::tibble(rsid=scz_EAS$rsID, pval=scz_EAS$P, id=scz_EAS$Phenotype),
  plink_bin = "plink",
  bfile = "1kg/EAS")

scz_EAS_clumped <- scz_EAS[scz_EAS$rsID %in% scz_EAS_clump$rsid,]

scz_EAS_clumped <- format_data(scz_EAS_clumped, 
                               type = "exposure",
                               phenotype_col = "Phenotype",
                               snp_col = "rsID",
                               beta_col = "BETA",
                               se_col = "SE",
                               eaf_col = "EAF",
                               effect_allele_col = "EA",
                               other_allele_col = "NEA",
                               pval_col = "P")

scz_EAS_exp_data <- scz_EAS_clumped[scz_EAS_clumped$pval.exposure < 1e-5,]
head(scz_EAS_exp_data)

# bmi_EAS
bmi_EAS_out_data <- bmi_EAS[bmi_EAS$rsID %in% scz_EAS_exp_data$SNP,]

bmi_EAS_out_data <- format_data(bmi_EAS_out_data, 
                                type = "outcome",
                                phenotype_col = "Phenotype",
                                snp_col = "rsID",
                                beta_col = "BETA",
                                se_col = "SE",
                                eaf_col = "EAF",
                                effect_allele_col = "EA",
                                other_allele_col = "NEA",
                                pval_col = "P")

# smkinit_EAS
smkinit_EAS_out_data <- smkinit_EAS[smkinit_EAS$rsID %in% scz_EAS_exp_data$SNP,]

smkinit_EAS_out_data <- format_data(smkinit_EAS_out_data, 
                                    type = "outcome",
                                    phenotype_col = "Phenotype",
                                    snp_col = "rsID",
                                    beta_col = "BETA",
                                    se_col = "SE",
                                    eaf_col = "EAF",
                                    effect_allele_col = "EA",
                                    other_allele_col = "NEA",
                                    pval_col = "P")


# cigday_EAS
cigday_EAS_out_data <- cigday_EAS[cigday_EAS$rsID %in% scz_EAS_exp_data$SNP,]

cigday_EAS_out_data <- format_data(cigday_EAS_out_data, 
                                   type = "outcome",
                                   phenotype_col = "Phenotype",
                                   snp_col = "rsID",
                                   beta_col = "BETA",
                                   se_col = "SE",
                                   eaf_col = "EAF",
                                   effect_allele_col = "EA",
                                   other_allele_col = "NEA",
                                   pval_col = "P")

EAS_out_data <- rbind(bmi_EAS_out_data, smkinit_EAS_out_data, cigday_EAS_out_data)

dat_scz_out_EAS <- harmonise_data(scz_EAS_exp_data, EAS_out_data, action = 2)

mr_scz_out_EAS <- mr(dat_scz_out_EAS, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))


# caluclate F-statistics
mr_scz_out_EAS$F <- NA

F_scz_bmi_EAS <- scz_EAS_exp_data[scz_EAS_exp_data$SNP %in% bmi_EAS_out_data$SNP,]
mr_scz_out_EAS$F[mr_scz_out_EAS$outcome=="bmi"] <- mean(F_scz_bmi_EAS$beta.exposure^2/F_scz_bmi_EAS$se.exposure^2)

F_scz_smkinit_EAS <- scz_EAS_exp_data[scz_EAS_exp_data$SNP %in% smkinit_EAS_out_data$SNP,]
mr_scz_out_EAS$F[mr_scz_out_EAS$outcome=="smkinit"] <- mean(F_scz_smkinit_EAS$beta.exposure^2/F_scz_smkinit_EAS$se.exposure^2)

F_scz_cigday_EAS <- scz_EAS_exp_data[scz_EAS_exp_data$SNP %in% cigday_EAS_out_data$SNP,]
mr_scz_out_EAS$F[mr_scz_out_EAS$outcome=="cigday"] <- mean(F_scz_cigday_EAS$beta.exposure^2/F_scz_cigday_EAS$se.exposure^2)

mr_EAS <- rbind(mr_exp_scz_EAS, mr_scz_out_EAS)

mr_EAS$pop <- "EAS"





############################### EUR

######################## exp to scz

# load bmi GWAS
bmi_EUR <- data.frame(fread("bmi_EUR.gwaslab.tsv.gz"))

bmi_EUR$Phenotype <- "bmi"

bmi_EUR_clump <- ld_clump(
  dplyr::tibble(rsid=bmi_EUR$rsID, pval=bmi_EUR$P, id=bmi_EUR$Phenotype),
  plink_bin = "plink",
  bfile = "1kg/EUR")

bmi_EUR_clumped <- bmi_EUR[bmi_EUR$rsID %in% bmi_EUR_clump$rsid,]

bmi_EUR_clumped <- format_data(bmi_EUR_clumped,
                               type = "exposure",
                               phenotype_col = "Phenotype",
                               snp_col = "rsID",
                               beta_col = "BETA",
                               se_col = "SE",
                               eaf_col = "EAF",
                               effect_allele_col = "EA",
                               other_allele_col = "NEA",
                               pval_col = "P")

bmi_EUR_exp_data <- bmi_EUR_clumped[bmi_EUR_clumped$pval.exposure < 5e-8,]
head(bmi_EUR_exp_data)


# load smkinit GWAS
smkinit_EUR <- data.frame(fread("smkinit_EUR.gwaslab.tsv.gz"))

smkinit_EUR$Phenotype <- "smkinit"

smkinit_EUR_clump <- ld_clump(
  dplyr::tibble(rsid=smkinit_EUR$rsID, pval=smkinit_EUR$P, id=smkinit_EUR$Phenotype),
  plink_bin = "plink",
  bfile = "1kg/EUR")

smkinit_EUR_clumped <- smkinit_EUR[smkinit_EUR$rsID %in% smkinit_EUR_clump$rsid,]

smkinit_EUR_clumped <- format_data(smkinit_EUR_clumped,
                                   type = "exposure",
                                   phenotype_col = "Phenotype",
                                   snp_col = "rsID",
                                   beta_col = "BETA",
                                   se_col = "SE",
                                   eaf_col = "EAF",
                                   effect_allele_col = "EA",
                                   other_allele_col = "NEA",
                                   pval_col = "P")

smkinit_EUR_exp_data <- smkinit_EUR_clumped[smkinit_EUR_clumped$pval.exposure < 5e-8,]
head(smkinit_EUR_exp_data)


# load cigday GWAS
cigday_EUR <- data.frame(fread("cigday_EUR.gwaslab.tsv.gz"))

cigday_EUR$Phenotype <- "cigday"

cigday_EUR_clump <- ld_clump(
  dplyr::tibble(rsid=cigday_EUR$rsID, pval=cigday_EUR$P, id=cigday_EUR$Phenotype),
  plink_bin = "plink",
  bfile = "1kg/EUR")

cigday_EUR_clumped <- cigday_EUR[cigday_EUR$rsID %in% cigday_EUR_clump$rsid,]

cigday_EUR_clumped <- format_data(cigday_EUR_clumped,
                                  type = "exposure",
                                  phenotype_col = "Phenotype",
                                  snp_col = "rsID",
                                  beta_col = "BETA",
                                  se_col = "SE",
                                  eaf_col = "EAF",
                                  effect_allele_col = "EA",
                                  other_allele_col = "NEA",
                                  pval_col = "P")

cigday_EUR_exp_data <- cigday_EUR_clumped[cigday_EUR_clumped$pval.exposure < 5e-8,]
head(cigday_EUR_exp_data)


# combine
EUR_exp_data <- rbind(bmi_EUR_exp_data, smkinit_EUR_exp_data, cigday_EUR_exp_data)

# load scz_EUR GWAS
scz_EUR <- data.frame(fread("scz_EUR.gwaslab.tsv.gz"))
scz_EUR$Phenotype <- "scz"

scz_EUR_out_data <- scz_EUR[scz_EUR$rsID %in% EUR_exp_data$SNP,]

scz_EUR_out_data <- format_data(scz_EUR_out_data, 
                                type = "outcome",
                                phenotype_col = "Phenotype",
                                snp_col = "rsID",
                                beta_col = "BETA",
                                se_col = "SE",
                                eaf_col = "EAF",
                                effect_allele_col = "EA",
                                other_allele_col = "NEA",
                                pval_col = "P")

dat_exp_scz_EUR <- harmonise_data(EUR_exp_data, scz_EUR_out_data, action = 2)

mr_exp_scz_EUR<- mr(dat_exp_scz_EUR, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))

# caluclate F-statistics
mr_exp_scz_EUR$F <- NA

F_bmi_scz_EUR <- bmi_EUR_exp_data[bmi_EUR_exp_data$SNP %in% scz_EUR_out_data$SNP,]
mr_exp_scz_EUR$F[mr_exp_scz_EUR$exposure=="bmi"] <- mean(F_bmi_scz_EUR$beta.exposure^2/F_bmi_scz_EUR$se.exposure^2)

F_smkinit_scz_EUR <- smkinit_EUR_exp_data[smkinit_EUR_exp_data$SNP %in% scz_EUR_out_data$SNP,]
mr_exp_scz_EUR$F[mr_exp_scz_EUR$exposure=="smkinit"]<- mean(F_smkinit_scz_EUR$beta.exposure^2/F_smkinit_scz_EUR$se.exposure^2)

F_cigday_scz_EUR <- cigday_EUR_exp_data[cigday_EUR_exp_data$SNP %in% scz_EUR_out_data$SNP,]
mr_exp_scz_EUR$F[mr_exp_scz_EUR$exposure=="cigday"] <- mean(F_cigday_scz_EUR$beta.exposure^2/F_cigday_scz_EUR$se.exposure^2)



################ scz to out

# scz
scz_EUR_clump <- ld_clump(
  dplyr::tibble(rsid=scz_EUR$rsID, pval=scz_EUR$P, id=scz_EUR$Phenotype),
  plink_bin = "plink",
  bfile = "1kg/EUR")

scz_EUR_clumped <- scz_EUR[scz_EUR$rsID %in% scz_EUR_clump$rsid,]

scz_EUR_clumped <- format_data(scz_EUR_clumped, 
                               type = "exposure",
                               phenotype_col = "Phenotype",
                               snp_col = "rsID",
                               beta_col = "BETA",
                               se_col = "SE",
                               eaf_col = "EAF",
                               effect_allele_col = "EA",
                               other_allele_col = "NEA",
                               pval_col = "P")

scz_EUR_exp_data <- scz_EUR_clumped[scz_EUR_clumped$pval.exposure < 5e-8,]
head(scz_EUR_exp_data)

# bmi_EUR
bmi_EUR_out_data <- bmi_EUR[bmi_EUR$rsID %in% scz_EUR_exp_data$SNP,]

bmi_EUR_out_data <- format_data(bmi_EUR_out_data, 
                                type = "outcome",
                                phenotype_col = "Phenotype",
                                snp_col = "rsID",
                                beta_col = "BETA",
                                se_col = "SE",
                                eaf_col = "EAF",
                                effect_allele_col = "EA",
                                other_allele_col = "NEA",
                                pval_col = "P")

# smkinit_EUR
smkinit_EUR_out_data <- smkinit_EUR[smkinit_EUR$rsID %in% scz_EUR_exp_data$SNP,]

smkinit_EUR_out_data <- format_data(smkinit_EUR_out_data, 
                                    type = "outcome",
                                    phenotype_col = "Phenotype",
                                    snp_col = "rsID",
                                    beta_col = "BETA",
                                    se_col = "SE",
                                    eaf_col = "EAF",
                                    effect_allele_col = "EA",
                                    other_allele_col = "NEA",
                                    pval_col = "P")


# cigday_EUR
cigday_EUR_out_data <- cigday_EUR[cigday_EUR$rsID %in% scz_EUR_exp_data$SNP,]

cigday_EUR_out_data <- format_data(cigday_EUR_out_data, 
                                   type = "outcome",
                                   phenotype_col = "Phenotype",
                                   snp_col = "rsID",
                                   beta_col = "BETA",
                                   se_col = "SE",
                                   eaf_col = "EAF",
                                   effect_allele_col = "EA",
                                   other_allele_col = "NEA",
                                   pval_col = "P")

EUR_out_data <- rbind(bmi_EUR_out_data, smkinit_EUR_out_data, cigday_EUR_out_data)

dat_scz_out_EUR <- harmonise_data(scz_EUR_exp_data, EUR_out_data, action = 2)

mr_scz_out_EUR <- mr(dat_scz_out_EUR, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))

# caluclate F-statistics
mr_scz_out_EUR$F <- NA

F_scz_bmi_EUR <- scz_EUR_exp_data[scz_EUR_exp_data$SNP %in% bmi_EUR_out_data$SNP,]
mr_scz_out_EUR$F[mr_scz_out_EUR$outcome=="bmi"] <- mean(F_scz_bmi_EUR$beta.exposure^2/F_scz_bmi_EUR$se.exposure^2)

F_scz_smkinit_EUR <- scz_EUR_exp_data[scz_EUR_exp_data$SNP %in% smkinit_EUR_out_data$SNP,]
mr_scz_out_EUR$F[mr_scz_out_EUR$outcome=="smkinit"] <- mean(F_scz_smkinit_EUR$beta.exposure^2/F_scz_smkinit_EUR$se.exposure^2)

F_scz_cigday_EUR <- scz_EUR_exp_data[scz_EUR_exp_data$SNP %in% cigday_EUR_out_data$SNP,]
mr_scz_out_EUR$F[mr_scz_out_EUR$outcome=="cigday"] <- mean(F_scz_cigday_EUR$beta.exposure^2/F_scz_cigday_EUR$se.exposure^2)

mr_EUR <- rbind(mr_exp_scz_EUR, mr_scz_out_EUR)

mr_EUR$pop <- "EUR"


mr_scz <- rbind(mr_EAS,mr_EUR)

write.csv(mr_scz, "mr_scz.csv", row.names=F)

save.image("mr_scz.RData")

