rm(list = ls())

library(tidyverse)
library(datawizard)

setwd("")

# load data

CKB_pheno <- readRDS("data_baseline_questionnaires.rds")

linkage <- readRDS("data_gwas_genetics.rds")
linkage <- linkage[,c("csid","national_pc01","national_pc02","national_pc03","national_pc04","national_pc05","national_pc06","national_pc07","national_pc08","national_pc09","national_pc10","national_pc11","is_in_gwas_population_subset")]
CKB_pheno <- merge(CKB_pheno,linkage,by="csid")

# merge pgs with dataset
pgs <- read.csv("pgs.csv",header = T)

CKB_pheno_pgs <- merge(CKB_pheno,pgs,by="csid")

# check follow-up time
CKB_pheno_pgs$study_date
CKB_pheno_pgs$censoring_date

CKB_pheno_pgs$follow_up_time <- as.numeric(difftime(CKB_pheno_pgs$censoring_date, CKB_pheno_pgs$study_date, units="days"))/365.25 
summary(CKB_pheno_pgs$follow_up_time)

# standardise based on population representative samples

for (i in c("pgs_scz_EAS","pgs_scz_EUR_EAS","pgs_mdd_EAS","pgs_mdd_EUR_EAS")){
  CKB_pheno_pgs[[i]] <- (CKB_pheno_pgs[[i]] -mean(CKB_pheno_pgs[[i]] [CKB_pheno_pgs$is_in_gwas_population_subset==1]))/sd(CKB_pheno_pgs[[i]] [CKB_pheno_pgs$is_in_gwas_population_subset==1])
}


# Section 1. Background information
#==================================

# Age
CKB_pheno_pgs$age <- as.numeric(CKB_pheno_pgs$age_at_study_date_x100/100)
summary(CKB_pheno_pgs$age)
CKB_pheno_pgs$age2 <- CKB_pheno_pgs$age*CKB_pheno_pgs$age

# Sex
CKB_pheno_pgs$is_female <- as.factor(CKB_pheno_pgs$is_female)
table(CKB_pheno_pgs$is_female, useNA="ifany")

# RCs
CKB_pheno_pgs$region_is_urban <- as.factor(CKB_pheno_pgs$region_is_urban)
table(CKB_pheno_pgs$region_is_urban, useNA="ifany")

CKB_pheno_pgs$region_code <- as.factor(CKB_pheno_pgs$region_code)
table(CKB_pheno_pgs$region_code, useNA="ifany")

# Education (>9 years of education)
table(CKB_pheno_pgs$highest_education)
CKB_pheno_pgs$high_school <- as.factor(ifelse(CKB_pheno_pgs$highest_education %in% c(3,4,5),1,0))
table(CKB_pheno_pgs$high_school, useNA="ifany")

# Occupation (employed vs. retired, house wife/husband, unemployed)
CKB_pheno_pgs$employed <- as.factor(ifelse(CKB_pheno_pgs$occupation %in% c(0:4,7),1,0))    
table(CKB_pheno_pgs$employed, useNA="ifany")

# Household income (>= 20,000)
CKB_pheno_pgs$income <- 0
CKB_pheno_pgs$income[CKB_pheno_pgs$household_income %in% c(4,5)] <- 1
CKB_pheno_pgs$income <- as.factor(CKB_pheno_pgs$income)
table(CKB_pheno_pgs$income, useNA="ifany")

# Household size
summary(CKB_pheno_pgs$household_size)

# number of children
summary(CKB_pheno_pgs$children)

# Ownership index 
CKB_pheno_pgs$ownership_index <- rowSums(CKB_pheno_pgs[, c("has_health_cover", "has_own_home", "has_private_toilet", "has_phone", "has_motor_vehicle", "had_recent_holiday")])
table(CKB_pheno_pgs$ownership_index, useNA="ifany")
summary(CKB_pheno_pgs$ownership_index, useNA="ifany")

# Marrital status
CKB_pheno_pgs$married_now <- as.factor(ifelse(CKB_pheno_pgs$marital_status==0,1,0))
table(CKB_pheno_pgs$married_now, useNA="ifany")

# Section 2. Lifestyle factors
#=============================
# physical activity
summary(CKB_pheno_pgs$met)

# Alcohol 
# Ever regular
CKB_pheno_pgs$ever_reg_alcohol <- as.factor(ifelse(CKB_pheno_pgs$alcohol_category %in% c(2,5,6),1,0))
table(CKB_pheno_pgs$ever_reg_alcohol, useNA="ifany")

# current regular
CKB_pheno_pgs$current_reg_alcohol <- as.factor(ifelse(CKB_pheno_pgs$alcohol_category==6,1,0))
table(CKB_pheno_pgs$current_reg_alcohol, useNA="ifany")

# age of alcohol initiation
summary(CKB_pheno_pgs$alc_age_start)

# Problem drinking among current drinkers
CKB_pheno_pgs$problem_drinking <- as.factor(ifelse(CKB_pheno_pgs$alcohol_category==6 & (CKB_pheno_pgs$alc_problem_cant_stop==1 | CKB_pheno_pgs$alc_problem_mental==1 | CKB_pheno_pgs$alc_problem_shakes==1 | CKB_pheno_pgs$alc_problem_work==1), 1, 0))
summary(CKB_pheno_pgs$problem_drinking)

# Smoking
# Ever regular
CKB_pheno_pgs$ever_reg_smoking <- as.factor(ifelse(CKB_pheno_pgs$smoking_category %in% c(3, 4), 1, 0))
table(CKB_pheno_pgs$ever_reg_smoking, useNA="ifany")

# Current regular
CKB_pheno_pgs$current_reg_smoking <- as.factor(ifelse(CKB_pheno_pgs$smoking_category==4, 1, 0))
table(CKB_pheno_pgs$current_reg_smoking, useNA="ifany")

# age of smoking initiation
# CKB_pheno_pgs$smoking_age_start[CKB_pheno_pgs$is_female==1] <- NA
summary(CKB_pheno_pgs$smoking_age_start)

# Amount of smoking
CKB_pheno_pgs$cig_equiv_day<-as.numeric(CKB_pheno_pgs$cig_equiv_day)
# CKB_pheno_pgs$cig_equiv_day[CKB_pheno_pgs$is_female==1] <- NA
summary(CKB_pheno_pgs$cig_equiv_day)

# Food diversity score
calculateDietVariables <- function(data) {
  food <- c("rice", "wheat", "other_staple", "meat", "poultry", "fish", "eggs", "dairy", "fresh_fruit", "fresh_veg", "preserved_veg",  "soybean")
  
  # Loop over food items
  for (ii in food) {
    freq_col <- paste0("diet_freq_", ii)
    diet_col <- paste0("diet_", ii)
    
    # Create diet variables
    data[[diet_col]] <- 0
    data[[diet_col]][data[[freq_col]] %in% c(2)] <- 1
    data[[diet_col]][data[[freq_col]] %in% c(0, 1)] <- 2
    
    # Display the table for the new diet variable
    print(table(data[[diet_col]]))
  }
  
  # Food diversity score
  data$food_diversity_score <- rowSums(data[, paste0("diet_", food)], na.rm = TRUE)
  
  return(data)
}
CKB_pheno_pgs <- calculateDietVariables(CKB_pheno_pgs)

# Section 3. Medical history and self-rated wellness
#=================================================== 

# Self-rated health  
CKB_pheno_pgs$self_rated_poor_heath<-as.factor(ifelse(CKB_pheno_pgs$self_rated_health==3,1,0))
table(CKB_pheno_pgs$self_rated_poor_heath, useNA="ifany")

# Comparative health      
CKB_pheno_pgs$comparative_poor_health <-as.factor(ifelse(CKB_pheno_pgs$comparative_health==2,1,0))
table(CKB_pheno_pgs$comparative_poor_health, useNA="ifany")

# asthma
CKB_pheno_pgs$asthma_diag <-as.factor(CKB_pheno_pgs$asthma_diag)

# emph
CKB_pheno_pgs$emph_bronc_diag <-as.factor(CKB_pheno_pgs$emph_bronc_diag)

# copd
CKB_pheno_pgs$has_copd <-as.factor(CKB_pheno_pgs$has_copd)

# tb 
CKB_pheno_pgs$tb_diag <-as.factor(CKB_pheno_pgs$tb_diag)

# kidney disease
CKB_pheno_pgs$kidney_dis_diag <-as.factor(CKB_pheno_pgs$kidney_dis_diag)

# cirrhosis_hep_diag
CKB_pheno_pgs$cirrhosis_hep_diag <-as.factor(CKB_pheno_pgs$cirrhosis_hep_diag)

# gall_diag
CKB_pheno_pgs$gall_diag <-as.factor(CKB_pheno_pgs$gall_diag)

# Diabetes 
CKB_pheno_pgs$has_diabetes<-as.factor(CKB_pheno_pgs$has_diabetes)
table(CKB_pheno_pgs$has_diabetes, useNA="ifany")

# Cancer 
CKB_pheno_pgs$cancer_diag<-as.factor(CKB_pheno_pgs$cancer_diag)  
table(CKB_pheno_pgs$cancer_diag, useNA="ifany")

# fracture
CKB_pheno_pgs$fracture_diag <-as.factor(CKB_pheno_pgs$fracture_diag)

# head injury
CKB_pheno_pgs$head_injury_diag <-as.factor(CKB_pheno_pgs$head_injury_diag)

# peptic_ulcer_diag
CKB_pheno_pgs$peptic_ulcer_diag <-as.factor(CKB_pheno_pgs$peptic_ulcer_diag)

# ra
CKB_pheno_pgs$rheum_arthritis_diag <-as.factor(CKB_pheno_pgs$rheum_arthritis_diag)

# rheum_heart_dis_diag
CKB_pheno_pgs$rheum_heart_dis_diag <-as.factor(CKB_pheno_pgs$rheum_heart_dis_diag)

# stroke
CKB_pheno_pgs$stroke_or_tia_diag <-as.factor(CKB_pheno_pgs$stroke_or_tia_diag)

# chd
CKB_pheno_pgs$chd_diag <-as.factor(CKB_pheno_pgs$chd_diag)




# Section 4. Mental health
#=========================

# Life satisfaction
CKB_pheno_pgs$life_unsatisfied<-as.factor(ifelse(CKB_pheno_pgs$satisfaction_level==3,1,0))
table(CKB_pheno_pgs$life_unsatisfied, useNA="ifany")    

# Stressful events
CKB_pheno_pgs$stressful_life_events <- factor(
  as.integer(CKB_pheno_pgs$sep_divorce == 1 | CKB_pheno_pgs$family_conflict == 1 |
               CKB_pheno_pgs$spouse_death_or_ill == 1 | CKB_pheno_pgs$other_family_death_or_ill == 1 |
               CKB_pheno_pgs$loss_of_income_debt == 1 | CKB_pheno_pgs$job_loss_retirement == 1 |
               CKB_pheno_pgs$bankruptcy == 1 | CKB_pheno_pgs$violence == 1 |
               CKB_pheno_pgs$injury_traffic_accident == 1 | CKB_pheno_pgs$natural_disaster == 1))
table(CKB_pheno_pgs$stressful_life_events, useNA="ifany")

# major depression
CKB_pheno_pgs$depression <- 0
CKB_pheno_pgs$depression[CKB_pheno_pgs$cidi_mdep==1] <- 1
CKB_pheno_pgs$depression <- as.factor(CKB_pheno_pgs$depression)

# depressive symptoms
CKB_pheno_pgs$dep_smp <- 0
CKB_pheno_pgs$dep_smp[CKB_pheno_pgs$feeling_worthless==1 |  CKB_pheno_pgs$loss_of_appetite==1 | CKB_pheno_pgs$loss_of_interest==1 | CKB_pheno_pgs$sad_depressed==1] <- 1
CKB_pheno_pgs$dep_smp <- factor(CKB_pheno_pgs$dep_smp)

# generalised anxiety disorder
CKB_pheno_pgs$anxiety <- 0
CKB_pheno_pgs$anxiety[CKB_pheno_pgs$cidi_gad==1] <- 1
CKB_pheno_pgs$anxiety <- factor(CKB_pheno_pgs$anxiety)

# anxiety symptoms
CKB_pheno_pgs$continuous_anxiety <- factor(CKB_pheno_pgs$continuous_anxiety)

# pain
CKB_pheno_pgs$continuous_pain <- factor(CKB_pheno_pgs$continuous_pain)

# panic attacks
CKB_pheno_pgs$panic_attacks <- factor(CKB_pheno_pgs$panic_attacks)

# phobia
CKB_pheno_pgs$phobic <- factor(CKB_pheno_pgs$phobic)

# neurasthenia
CKB_pheno_pgs$neurasthenia_diag <- as.factor(CKB_pheno_pgs$neurasthenia_diag)

# psychiatric disorder
CKB_pheno_pgs$psych_disorder_diag <- as.factor(CKB_pheno_pgs$psych_disorder_diag)

# Sleep duration
CKB_pheno_pgs$sleep_hours<-as.numeric(CKB_pheno_pgs$sleep_hours)
summary(CKB_pheno_pgs$sleep_hours, useNA="ifany")    

# Snoring
CKB_pheno_pgs$snoring<-as.factor(ifelse(CKB_pheno_pgs$sleep_snoring==0,1,0))
table(CKB_pheno_pgs$snoring,useNA="ifany")

# sleep problems
CKB_pheno_pgs$sleep_problems = if_else(CKB_pheno_pgs$sleep_affects_day==1 | CKB_pheno_pgs$sleep_delayed_fitful==1 | CKB_pheno_pgs$sleep_waking_too_early==1, 1, 0)
CKB_pheno_pgs$sleep_problems <- as.factor(CKB_pheno_pgs$sleep_problems)


# Section 5. Clinical measurements
#================================ 

# BMI
CKB_pheno_pgs$bmi <- CKB_pheno_pgs$bmi_calc
summary(CKB_pheno_pgs$bmi)

# Standing height  
CKB_pheno_pgs$standing_height <- CKB_pheno_pgs$standing_height_mm/10
summary(CKB_pheno_pgs$standing_height)

# Sitting height  
CKB_pheno_pgs$sitting_height <- CKB_pheno_pgs$sitting_height_mm/10
summary(CKB_pheno_pgs$sitting_height)

# SBP
CKB_pheno_pgs$sbp <- CKB_pheno_pgs$sbp_mean
summary(CKB_pheno_pgs$sbp)

# DBP
CKB_pheno_pgs$dbp <- CKB_pheno_pgs$dbp_mean
summary(CKB_pheno_pgs$dbp)

# Heart rate
CKB_pheno_pgs$heart_rate <- CKB_pheno_pgs$heart_rate_mean
summary(CKB_pheno_pgs$heart_rate)

# Random plasma glucose 
CKB_pheno_pgs$random_glucose <- CKB_pheno_pgs$random_glucose_x10/10
summary(CKB_pheno_pgs$random_glucose)

# Section 6. reproductive history

# age at menarche
summary(CKB_pheno_pgs$first_period_age)

# post-menopause
table(CKB_pheno_pgs$menopause_status)
CKB_pheno_pgs$post_menopause <- as.factor(ifelse(CKB_pheno_pgs$menopause_status==2,1,0))
summary(CKB_pheno_pgs$post_menopause)

# age at menopause
summary(CKB_pheno_pgs$menopause_age)

# had hysterectomy
CKB_pheno_pgs$had_hysterectomy <- as.factor(CKB_pheno_pgs$had_hysterectomy)
summary(CKB_pheno_pgs$had_hysterectomy)

# had breast lump removed
CKB_pheno_pgs$had_breast_lump_removed <- as.factor(CKB_pheno_pgs$had_breast_lump_removed)
summary(CKB_pheno_pgs$had_breast_lump_removed)

# had ovary removed
CKB_pheno_pgs$had_ovary_removed <- as.factor(CKB_pheno_pgs$had_ovary_removed)
summary(CKB_pheno_pgs$had_ovary_removed)

# number of pregnancies
summary(CKB_pheno_pgs$preg_count)

# # parity
# CKB_pheno_pgs$parity <- as.factor(ifelse(CKB_pheno_pgs$live_birth_count>0 | CKB_pheno_pgs$still_birth_count>0, 1, 0))
# summary(CKB_pheno_pgs$parity)

# number of live birth
summary(CKB_pheno_pgs$live_birth_count)

# age at first live birth
summary(CKB_pheno_pgs$live_birth_age01)

# number of still birth
summary(CKB_pheno_pgs$still_birth_count)


saveRDS(CKB_pheno_pgs, "baseline_pgs.rds")
