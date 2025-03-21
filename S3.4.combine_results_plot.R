rm(list = ls())

library(ggplot2)
library(ggrepel)
library(stringr)
library(dplyr)
library(ggpubr)

setwd("")

### combine results

dis <- c("scz", "mdd")

pop <- c("EAS","EUR_EAS")

### with exclusion

results <- list()

c <- 1

for (d in 1:2){
  
  for (p in 1:2){
    
    assoc_pgs_baseline <- read.csv(paste0("assoc_baseline_pgs_", dis[d], "_", pop[p], ".csv"))
    # get female-specific variables
    assoc_pgs_baseline_female <- read.csv(paste0("assoc_baseline_pgs_", dis[d], "_", pop[p], "_female.csv"))
    variable_baseline_female <- setdiff(assoc_pgs_baseline_female$variable, assoc_pgs_baseline$variable)
    assoc_pgs_baseline <- rbind(assoc_pgs_baseline, assoc_pgs_baseline_female[assoc_pgs_baseline_female$variable %in% variable_baseline_female,])
    assoc_pgs_baseline$dis <- dis[d]

    # name
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="region_is_urban"] <-"Urban area"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="life_unsatisfied"] <-"Unsatisfied with life"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="stressful_life_events"] <-"Stressful life events"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="depression"] <-"Major depression"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="dep_smp"] <-"Depressive symptoms"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="anxiety"] <-"Generalised anxiety disorder"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="continuous_anxiety"] <-"Continous anxiety"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="continuous_pain"] <-"Continuous pain"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="panic_attacks"] <-"Panic attacks"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="phobic"] <-"Phobia symptoms"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="neurasthenia_diag"] <-"Neurasthenia"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="psych_disorder_diag"] <-"Psychiatric disorder"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="sleep_hours"] <-"Sleep duration"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="snoring"] <-"Snoring"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="sleep_problems"] <-"Sleep problems"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="self_rated_poor_heath"] <-"Poor self-rated health"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="comparative_poor_health"] <-"Worse comparative health"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="asthma_diag"] <-"Asthma"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="emph_bronc_diag"] <-"Emphysema/bronchitis"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="has_copd"] <-"Chronic obstructive pulmonary disease"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="tb_diag"] <-"Tuberculosis"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="kidney_dis_diag"] <-"Kidney disease"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="cirrhosis_hep_diag"] <- "Cirrhosis/chronic hepatitis"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="gall_diag"] <- "Gallstone/gallbladder disease"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="has_diabetes"] <- "Diabetes"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="cancer_diag"] <- "Cancer"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="fracture_diag"] <- "Fracture"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="head_injury_diag"] <- "Head injury"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="peptic_ulcer_diag"] <- "Peptic ulcer"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="rheum_arthritis_diag"] <- "Rheumatoid arthritis"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="rheum_heart_dis_diag"] <- "Rheumatoid heart disease"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="stroke_or_tia_diag"] <- "Stroke/transient ischemic attack"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="chd_diag"] <- "Chronic heart disease"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="bmi"] <- "Body mass index"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="standing_height"] <- "Standing height"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="sitting_height"] <- "Sitting height"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="sbp"] <- "Systolic blood pressure"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="dbp"] <- "Dystolic blood pressure"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="heart_rate"] <- "Heart rate"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="random_glucose"] <- "Random glucose"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="high_school"] <- "> 9 years education"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="employed"] <- "Employed"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="income"] <- "Household income >= 20,000"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="married_now"] <- "Married"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="ownership_index"] <- "Ownership index"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="household_size"] <- "Household size"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="children"] <- "Number of children"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="met"] <- "Physical activity"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="ever_reg_alcohol"] <- "Ever regular alcohol drinking"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="current_reg_alcohol"] <- "Current regular alcohol drinking"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="alc_age_start"] <- "Age of alcohol initiation "
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="problem_drinking"] <- "Problem alcohol drinking"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="ever_reg_smoking"] <- "Ever regular smoking"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="current_reg_smoking"] <- "Current regular smoking"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="smoking_age_start"] <- "Age of smoking initation"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="cig_equiv_day"] <- "Cigarettes per day"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="food_diversity_score"] <- "Food diversity"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="first_period_age"] <- "Age at menarche"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="post_menopause"] <- "Post menopause"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="menopause_age"] <- "Age at menopause"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="preg_count"] <- "Number of pregnancies"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="live_birth_count"] <- "Number of live births"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="live_birth_age01"] <- "Age at first live birth"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="still_birth_count"] <- "Number of still births"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="had_hysterectomy"] <- "Had hysterectomy"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="had_breast_lump_removed"] <- "Had breast lump removed"
    assoc_pgs_baseline$description[assoc_pgs_baseline$variable=="had_ovary_removed"] <- "Had ovary removed"
    
    
    assoc_pgs_baseline$description <- factor(assoc_pgs_baseline$description, 
                                             levels=c("Urban area",
                                                      "Married", 
                                                      "Employed", 
                                                      "> 9 years education", 
                                                      "Household size",
                                                      "Number of children",
                                                      "Ownership index",
                                                      "Household income >= 20,000",
                                                      "Major depression",
                                                      "Depressive symptoms",
                                                      "Generalised anxiety disorder",
                                                      "Continous anxiety",
                                                      "Continuous pain",
                                                      "Panic attacks",
                                                      "Phobia symptoms",
                                                      "Neurasthenia",
                                                      "Psychiatric disorder",
                                                      "Unsatisfied with life",
                                                      "Stressful life events",
                                                      "Sleep duration",
                                                      "Sleep problems",
                                                      "Snoring",
                                                      "Poor self-rated health",
                                                      "Worse comparative health",
                                                      "Asthma",
                                                      "Chronic obstructive pulmonary disease",
                                                      "Emphysema/bronchitis",
                                                      "Tuberculosis",
                                                      "Rheumatoid arthritis",
                                                      "Rheumatoid heart disease",
                                                      "Chronic heart disease",
                                                      "Stroke/transient ischemic attack",
                                                      "Diabetes",
                                                      "Peptic ulcer",
                                                      "Cirrhosis/chronic hepatitis",
                                                      "Gallstone/gallbladder disease",
                                                      "Kidney disease",
                                                      "Head injury",
                                                      "Fracture",
                                                      "Cancer",
                                                      "Body mass index", 
                                                      "Standing height", 
                                                      "Sitting height", 
                                                      "Systolic blood pressure", 
                                                      "Dystolic blood pressure", 
                                                      "Heart rate", 
                                                      "Random glucose",
                                                      "Physical activity",
                                                      "Food diversity",
                                                      "Ever regular alcohol drinking", 
                                                      "Current regular alcohol drinking",
                                                      "Age of alcohol initiation ",
                                                      "Problem alcohol drinking",
                                                      "Ever regular smoking",
                                                      "Current regular smoking",
                                                      "Age of smoking initation", 
                                                      "Cigarettes per day",
                                                      "Age at menarche",
                                                      "Post menopause",
                                                      "Age at menopause",
                                                      "Number of pregnancies",
                                                      "Number of live births",
                                                      "Age at first live birth",
                                                      "Number of still births",
                                                      "Had hysterectomy",
                                                      "Had breast lump removed",
                                                      "Had ovary removed"))
    
    assoc_pgs_baseline$group[assoc_pgs_baseline$variable %in% c("life_unsatisfied", 
                                                                "stressful_life_events", 
                                                                "depression", 
                                                                "dep_smp",
                                                                "continuous_anxiety",
                                                                "continuous_pain",
                                                                "panic_attacks",
                                                                "phobic",
                                                                "anxiety", 
                                                                "neurasthenia_diag",
                                                                "psych_disorder_diag", 
                                                                "sleep_hours",
                                                                "snoring", 
                                                                "sleep_problems")] <- "Mental health and sleep"
    
    assoc_pgs_baseline$group[assoc_pgs_baseline$variable %in% c("self_rated_poor_heath",
                                                                "comparative_poor_health",
                                                                "asthma_diag",
                                                                "emph_bronc_diag",
                                                                "has_copd",
                                                                "tb_diag",
                                                                "kidney_dis_diag", 
                                                                "cirrhosis_hep_diag", 
                                                                "gall_diag",
                                                                "has_diabetes",
                                                                "cancer_diag", 
                                                                "fracture_diag", 
                                                                "head_injury_diag", 
                                                                "peptic_ulcer_diag",
                                                                "rheum_arthritis_diag", 
                                                                "rheum_heart_dis_diag", 
                                                                "stroke_or_tia_diag", 
                                                                "chd_diag")] <- "Physical health"
    
    assoc_pgs_baseline$group[assoc_pgs_baseline$variable %in% c("bmi", "standing_height", "sitting_height", 
                                                                "sbp", "dbp", "heart_rate", "random_glucose")] <- "Clinical measurements"
    
    assoc_pgs_baseline$group[assoc_pgs_baseline$variable %in% c("region_is_urban", "high_school", "employed", "income",
                                                                "married_now", "ownership_index", "household_size", "children")] <- "Socio-demographics"
    
    assoc_pgs_baseline$group[assoc_pgs_baseline$variable %in% c("met", "food_diversity_score", "ever_reg_alcohol", "current_reg_alcohol", "alc_age_start", "problem_drinking",
                                                                "ever_reg_smoking", "current_reg_smoking", "smoking_age_start", "cig_equiv_day")] <- "Lifestyle"
    
    assoc_pgs_baseline$group[assoc_pgs_baseline$variable %in% c("first_period_age",
                                                                "post_menopause",
                                                                "menopause_age",
                                                                "preg_count",
                                                                "live_birth_count",
                                                                "live_birth_age01",
                                                                "still_birth_count",                                                                "Age at first live birth",
                                                                "had_hysterectomy",
                                                                "had_breast_lump_removed",
                                                                "had_ovary_removed")] <- "Reproductive history"
    
    assoc_pgs_baseline$time <- "Baseline"
    
    ##### load phecode results and combine
    phewas_results <- read.csv(paste0("phewas_pgs_", dis[d], "_", pop[p], "_parent.csv"))
    
    # get male and female specific phecodes
    phewas_results_male <- read.csv(paste0("phewas_pgs_", dis[d], "_", pop[p], "_male_parent.csv"))
    variable_phecode_male <- setdiff(phewas_results_male$phenotype, phewas_results$phenotype)
    phewas_results <- rbind(phewas_results, phewas_results_male[phewas_results_male$phenotype %in% variable_phecode_male,])
    phewas_results_female <- read.csv(paste0("phewas_pgs_", dis[d], "_", pop[p], "_female_parent.csv"))
    variable_phecode_female <- setdiff(phewas_results_female$phenotype, phewas_results$phenotype)
    phewas_results <- rbind(phewas_results, phewas_results_female[phewas_results_female$phenotype %in% variable_phecode_female,])
    phewas_results$dis <- dis[d]
    
    phewas_results$time <- "Follow-up"
    
    phewas_results <- phewas_results[c("phenotype", "beta", "SE", "OR", "p", "n_total", "n_cases", "n_controls", "description", "group", "dis", "time")]
    
    baseline_results <- assoc_pgs_baseline[c("variable", "beta", "SE", "OR", "p", "n_total", "n_cases", "n_controls", "description", "group", "dis", "time")]
    names(phewas_results) <- names(baseline_results)
    
    results[[c]] <- rbind(baseline_results, phewas_results)
    results[[c]]$group <- str_to_sentence(results[[c]]$group)
    results[[p]]$group <- factor(results[[p]]$group, levels=c("Socio-demographics","Physical health","Mental health and sleep", "Clinical measurements", "Lifestyle", "Reproductive history",
                                                              "Infectious diseases","Neoplasms","Endocrine/metabolic","Hematopoietic","Mental disorders","Neurological","Sense organs","Circulatory system",
                                                              "Respiratory","Digestive","Genitourinary","Dermatologic","Musculoskeletal","Injuries & poisonings","Symptoms"))
    results[[c]]$description <- str_to_sentence(results[[c]]$description)
    
    results[[c]]$pop <- pop[p]
    
    c <- c+1
    
  }
  
}

results_combined <- rbind(results[[1]], results[[2]], results[[3]], results[[4]])

results_combined$p_fdr <- p.adjust(results_combined$p, method="fdr")

# save
write.csv(results_combined, "results_combined.csv", row.names=F)

##### plots

plots <- list()

# plot_title <- c("PGS-SCZ", "PGS-MD")

for (d in 1:2){
  
  plot_data <- results_combined[results_combined$dis==dis[d],]
  
  # Calculate the -log10 of the p-values
  plot_data$LogP = -log10(plot_data$p)
  
  # Assign a positive or negative sign to the LogP values based on the sign of the Estimate
  plot_data$SignedLogP = ifelse(plot_data$beta > 0, plot_data$LogP, -plot_data$LogP)
  
  # Recreate the SignGroup variable
  plot_data$SignGroup <- ifelse(plot_data$p_fdr < 0.05, ifelse(plot_data$beta > 0, "Positive", "Inverse"), "Non-significant")

  # Make sure it's a factor, especially if you want to control the order of legend items or have specific factor levels
  plot_data$SignGroup <- factor(plot_data$SignGroup, levels = c("Positive", "Inverse", "Non-significant"))
  
  # labels
  plot_data$label <- as.character(plot_data$description)
  plot_data$label[plot_data$p_fdr>=0.05] <- NA
  
  # only keep top 2 in each category
  for (i in levels(plot_data$group)){
    df <- plot_data[plot_data$group==i,]
    df$label[order(df$p)[3:nrow(df)]] <- NA
    plot_data$label[plot_data$group==i] <- df$label
  }
  
  plot_data$source <- ifelse(plot_data$pop=="EAS", "EAS", "Multi (EAS + EUR)")
  
  # plot
  options(ggrepel.max.overlaps = Inf)
  plots[[d]] <- ggplot(plot_data, aes(x = group, y = SignedLogP, shape=source, color = SignGroup,label=str_wrap(label, width = 35))) +
    geom_point(alpha = 0.5) +
    scale_y_continuous(labels = abs, limits = c(-8.5,13.5)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    scale_shape_manual(values=c(1, 16)) +
    scale_color_manual(values = c("Positive" = "red", "Inverse" = "blue", "Non-significant" = "gray")) +
    geom_text_repel(size = 3.5, color = "black", xlim = c(0,1000), force = 50) +
    labs(y = "-log10(p-value)", x = "Phenotype category") +
    theme_minimal() +
    # ggtitle(plot_title[d]) +
    labs(shape="GWAS source", color="Direction of association") +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line(color = "gray")) +
    facet_grid(~time, scales = "free_x", space = "free_x", switch = "x") +
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle=60, hjust = 1)) +
    coord_cartesian(clip = "off")
  
  if (d==1){
    plots[[d]] <- plots[[d]] +theme(axis.title.x=element_blank(), axis.text.x=element_blank())
  }
}

plots_combined <- ggarrange(plots[[1]], plots[[2]], ncol = 1, nrow = 2, legend = "right", common.legend = T, labels="auto", font.label = list(size = 16))

ggsave("mirror_manhattan_plot.jpg", plots_combined, width = 8, height = 12)
