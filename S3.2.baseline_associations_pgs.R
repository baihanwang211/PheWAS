rm(list = ls())

library(ggplot2)

setwd("")

# read

baseline_pgs <- readRDS("baseline_pgs.rds")

pgs <- c("pgs_scz_EAS","pgs_scz_EUR_EAS","pgs_mdd_EAS","pgs_mdd_EUR_EAS")


for (x in 1:4){
  
  # continuous variables
  
  variable_con <- c("ownership_index", "household_size", "children", "met", "alc_age_start", "smoking_age_start", "cig_equiv_day", "food_diversity_score",
                    "sleep_hours", "bmi", "standing_height", "sitting_height", "sbp", "dbp", "heart_rate", "random_glucose")

  assoc_pgs_baseline_con <- data.frame(variable=variable_con,beta=rep(NA,length(variable_con)),SE=rep(NA,length(variable_con)),OR=rep(NA,length(variable_con)),p=rep(NA,length(variable_con)),
                                       n_total=rep(NA,length(variable_con)),n_cases=rep(NA,length(variable_con)),n_controls=rep(NA,length(variable_con)))
  
  for (i in 1:nrow(assoc_pgs_baseline_con)){
    data <- baseline_pgs[baseline_pgs$is_in_gwas_population_subset==1,]
    model <- lm(paste0(variable_con[i], " ~ ", pgs[x], " + age + age2 + is_female + region_code +
                    national_pc01 + national_pc02 + national_pc03 + national_pc04 + national_pc05 + national_pc06 + national_pc07 + national_pc08 + national_pc09 + national_pc10 + national_pc11"), data)
    assoc_pgs_baseline_con$beta[i] <- model$coefficients[2]
    assoc_pgs_baseline_con$SE[i] <- summary(model)$coefficients[2,"Std. Error"]
    assoc_pgs_baseline_con$p[i] <- summary(model)$coefficients[2,"Pr(>|t|)"]  
    assoc_pgs_baseline_con$n_total[i] <- sum(!is.na(data[variable_con[i]]))
    print(i)
  }
  
  # categorical variables
  
  variable_cat <- c("region_is_urban", "high_school", "employed", "income", "married_now", "ever_reg_alcohol", "current_reg_alcohol", 
                     "problem_drinking", "ever_reg_smoking", "current_reg_smoking", "self_rated_poor_heath", "comparative_poor_health",
                    "life_unsatisfied", "stressful_life_events","snoring", "sleep_problems")
  
  assoc_pgs_baseline_cat <-data.frame(variable=variable_cat,beta=rep(NA,length(variable_cat)),SE=rep(NA,length(variable_cat)),OR=rep(NA,length(variable_cat)),p=rep(NA,length(variable_cat)),
                                      n_total=rep(NA,length(variable_cat)),n_cases=rep(NA,length(variable_cat)),n_controls=rep(NA,length(variable_cat)))
  
  for (i in 1:nrow(assoc_pgs_baseline_cat)){
    data <- baseline_pgs[baseline_pgs$is_in_gwas_population_subset==1,]
    if (variable_cat[i]=="region_is_urban"){
      model <- glm(paste0(variable_cat[i], " ~ ", pgs[x], " + age + age2 + is_female +
                    national_pc01 + national_pc02 + national_pc03 + national_pc04 + national_pc05 + national_pc06 + national_pc07 + national_pc08 + national_pc09 + national_pc10 + national_pc11"), family = "binomial", data)
      
    } else {
      model <- glm(paste0(variable_cat[i], " ~ ", pgs[x], " + age + age2 + is_female + region_code +
                          national_pc01 + national_pc02 + national_pc03 + national_pc04 + national_pc05 + national_pc06 + national_pc07 + national_pc08 + national_pc09 + national_pc10 + national_pc11"), family = "binomial", data)
    }
    assoc_pgs_baseline_cat$beta[i] <- model$coefficients[2]
    assoc_pgs_baseline_cat$SE[i] <- summary(model)$coefficients[2,"Std. Error"]
    assoc_pgs_baseline_cat$OR[i] <- exp(assoc_pgs_baseline_cat$beta[i])
    assoc_pgs_baseline_cat$p[i] <- summary(model)$coefficients[2,"Pr(>|z|)"]  
    assoc_pgs_baseline_cat$n_total[i] <- sum(!is.na(data[variable_cat[i]]))
    assoc_pgs_baseline_cat$n_cases[i] <- table(data[variable_cat[i]])[2]
    assoc_pgs_baseline_cat$n_controls[i] <- table(data[variable_cat[i]])[1]
    print(i)
  }
  
  # disease-related variables
  
  variable_dis <- c("asthma_diag", "emph_bronc_diag", "has_copd", "tb_diag", "kidney_dis_diag", 
                   "cirrhosis_hep_diag", "gall_diag", "has_diabetes", "cancer_diag", "fracture_diag", "head_injury_diag", "peptic_ulcer_diag", 
                   "rheum_arthritis_diag", "rheum_heart_dis_diag", "stroke_or_tia_diag", "chd_diag", 
                   "depression", "dep_smp", "anxiety", "continuous_anxiety", "continuous_pain", "panic_attacks", "phobic", "neurasthenia_diag", "psych_disorder_diag")
  
  assoc_pgs_baseline_dis <-data.frame(variable=variable_dis,beta=rep(NA,length(variable_dis)),SE=rep(NA,length(variable_dis)),OR=rep(NA,length(variable_dis)),p=rep(NA,length(variable_dis)),
                                      n_total=rep(NA,length(variable_dis)),n_cases=rep(NA,length(variable_dis)),n_controls=rep(NA,length(variable_dis)))
  
  for (i in 1:nrow(assoc_pgs_baseline_dis)){
    # keep all cases and controls in subcohort
    data_1 <- baseline_pgs[baseline_pgs[variable_dis[i]]==1,]
    data_2 <- baseline_pgs[baseline_pgs[variable_dis[i]]==0 & baseline_pgs$is_in_gwas_population_subset==1,]
    data <- rbind(data_1, data_2)
    model <- glm(paste0(variable_dis[i], " ~ ", pgs[x], " + age + age2 + is_female + region_code +
                          national_pc01 + national_pc02 + national_pc03 + national_pc04 + national_pc05 + national_pc06 + national_pc07 + national_pc08 + national_pc09 + national_pc10 + national_pc11"), family = "binomial", data)
    assoc_pgs_baseline_dis$beta[i] <- model$coefficients[2]
    assoc_pgs_baseline_dis$SE[i] <- summary(model)$coefficients[2,"Std. Error"]
    assoc_pgs_baseline_dis$OR[i] <- exp(assoc_pgs_baseline_dis$beta[i])
    assoc_pgs_baseline_dis$p[i] <- summary(model)$coefficients[2,"Pr(>|z|)"]  
    assoc_pgs_baseline_dis$n_total[i] <- sum(!is.na(data[variable_dis[i]]))
    assoc_pgs_baseline_dis$n_cases[i] <- table(data[variable_dis[i]])[2]
    assoc_pgs_baseline_dis$n_controls[i] <- table(data[variable_dis[i]])[1]
    print(i)
  }
  
  # combind all results and save
  
  assoc_pgs_baseline <- rbind(assoc_pgs_baseline_con,assoc_pgs_baseline_cat, assoc_pgs_baseline_dis)
  
  write.csv(assoc_pgs_baseline, paste0("assoc_baseline_", pgs[x], ".csv"), row.names=F)
  
}
