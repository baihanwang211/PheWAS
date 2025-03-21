rm(list = ls())

# load packages 
library(tidyverse)
library(plyr)
library(parallel)
library(PheWAS)

options(future.globals.maxSize = 8000 * 1024^2)
setwd("")

# load phecode data
load("processed_phecodes_parent.RData")

# Save the SNPs to a vector:
snp_list <- c("pgs_scz_EAS",,"pgs_scz_EUR_EAS","pgs_mdd_EAS",,"pgs_mdd_EUR_EAS")

## phewas with exclusion

# Function for applying the log model:
apply_log_model <- function(phe, df, snp){
  
  # Find exclusions:
  phe_exclusions <- PheWAS::phecode_exclude 
  
  # convert to ckb format for matching
  phe_exclusions$ckb_phecode <- gsub(".", "_", phe_exclusions$code, fixed = TRUE)
  phe_exclusions$ckb_phecode <- paste0("phe_0", phe_exclusions$ckb_phecode)
  phe_exclusions$ckb_exclusion <- gsub(".", "_", phe_exclusions$exclusion_criteria, fixed = TRUE)
  phe_exclusions$ckb_exclusion <- paste0("phe_0", phe_exclusions$ckb_exclusion)
  
  exclude <- phe_exclusions$ckb_exclusion[phe_exclusions$ckb_phecode == phe]
  exclude <- exclude[exclude %in% colnames(df)] # keep only those that exist in the df
  
  ## select cases from whole cohort, and controls from within population subset only
  tab <- df[df$is_in_gwas_population_subset == 1 &
              !is.na(df$is_in_gwas_population_subset) |
              df[[phe]] == 1, ]
  
  # split cases and controls
  tab_cases <- tab[tab[[phe]]==1,]
  tab_controls <- tab[tab[[phe]]==0,]
  
  # exclude any controls that have an event of the exclusion phecode - this needs to be an if statement otherwise there is an error if there are no controls to exclude
  if (sum(tab_controls[exclude]) > 0) {
    tab_controls <- tab_controls %>%
      filter_at(vars(exclude), all_vars(. == 0))
  }
  
  # recombine
  tab <- rbind(tab_cases,tab_controls)
  
  # Both sexes are being analysed together, need to adjust for sex:
  adjust <- c("sex", "age", "age2", "region_code", paste0("national_pc", 1:11))
  
  # Fit logistic model:
  model <- summary(glm(as.formula(paste0(phe, " ~ ", snp, " + ",
                                         paste0(adjust, collapse = " + "))),
                       data = tab, 
                       family = "binomial"))
  # Store results:
  beta <- model$coef[names(model$aliased) == snp, 1]
  se <- model$coef[names(model$aliased) == snp, 2]
  p <- model$coef[names(model$aliased) == snp, 4]
  phe_table <- gsub("_", ".", gsub("phe_0", "", phe))
  phecode_output <- data.frame(phe_table, snp, beta, se, exp(beta), p, nrow(tab), 
                               table(tab[[phe]])[2], table(tab[[phe]])[1])
  names(phecode_output) <- c("phenotype", "snp", "beta", "SE", "OR", "p", "n_total", 
                             "n_cases", "n_controls")
  return(phecode_output)
}


# Function for running it over each SNP:
run_phewas_over_snps <- function(snp, df){
  
  # Select only those with SNP data
  df  <- df[!is.na(df[[snp]]), ] 
  # get sex restriction info
  phe_sex_info <- PheWAS::sex_restriction
  # convert to ckb format for matching
  phe_sex_info$ckb_phecode <- gsub(".", "_", phe_sex_info$phecode, fixed = TRUE)
  phe_sex_info$ckb_phecode <- paste0("phe_0", phe_sex_info$ckb_phecode)
  
  female_only <- phe_sex_info$ckb_phecode[phe_sex_info$female_only=="TRUE"]
  male_only <- phe_sex_info$ckb_phecode[phe_sex_info$male_only=="TRUE"]
  
  phe_cols <- names(df)[grep("phe_", names(df))]
  
  combined_phe <- phe_cols[!(phe_cols %in% female_only | phe_cols %in% male_only)]
  
  # Apply log model on every phecode
  snp_result_list <-future.apply::future_lapply(combined_phe,
                                                apply_log_model,
                                                df,
                                                snp)
  
  snp_output_tab <- do.call(rbind, snp_result_list)
  
  return(snp_output_tab)
}

phewas_results_list_2 <- future.apply::future_lapply(snp_list,     # run over list of snps
                                                     run_phewas_over_snps, # run function
                                                     all_phe_100_cases)    # use the df containing only phecodes >= 100 cases           

saveRDS(phewas_results_list_2, "phewas_pgs_list_parent.rds")

phewas_results_list_2_annot <- phewas_results_list_2

for (i in 1:4){
  phewas_results_list_2_annot[[i]] <- addPhecodeInfo(phewas_results_list_2_annot[[i]], 
                                                     groupnums = T, # add colours for plotting 
                                                     groupcolors = T)
  write.csv(phewas_results_list_2_annot[[i]], paste0("phewas_", snp_list[i], "_parent.csv"), row.names = F)
}