
rm(list = ls())

library(rcompanion)
library(r2redux)
library(dplyr)
library(ckbplotr)
library(ggplot2)
library(ggpubr)

# load data
setwd("J:/psych/data/mdd_new")

baseline_pgs <- readRDS("baseline_pgs.rds")

# load endpoints
endpoints <- readRDS("endpoints.rds")
data <- merge(baseline_pgs,endpoints,by="csid")

# define mdd
data$mdd <- 0
data$mdd[data$depression==1 | data$ep_diy_p2130194614_combined_ep==1] <- 1

table(data$depression)
table(data$ep_diy_p2130194614_combined_ep)
table(data$mdd)
intersect(data$csid[data$depression==1], data$csid[data$ep_diy_p2130194614_combined_ep==1])

data$scz <- data$ep_Schizophrenia_combined_ep

dis <- c("scz", "mdd")
dis_label <- c("SCZ", "MD")
pop <- c("EAS","EUR_EAS")
pop_label <- c("EAS", "multi")
prevalence <- c(0.006, 0.018)

assoc_pgs <- data.frame(variable=rep(NA,4),beta=rep(NA,4),SE=rep(NA,4),p=rep(NA,4),
                        or=rep(NA,4),lci=rep(NA,4),uci=rep(NA,4),
                        r2_nagelkerke=rep(NA,4),r2_liability=rep(NA,4))

c <- 1
plot <- list()

for (d in 1:2){
  
  df_1 <- data[data[dis[d]]==1,]
  df_2 <- data[data[dis[d]]==0 & data$is_in_gwas_population_subset==1,]
  df <- rbind(df_1, df_2)
  
  n_total <- nrow(df)
  n_case <- table(df[dis[d]])[2]
  n_control <- table(df[dis[d]])[1]

  hist <- list()
  plot_quartile <- list()

  for (p in 1:2){
    
    pgs <- paste("pgs",dis[d],pop[p],sep="_")

    # histogram
    hist[[p]] <- ggplot(df, aes(x = .data[[pgs]], fill=factor(.data[[dis[d]]], levels=c(1,0)))) + 
      geom_density(alpha=0.5) +
      # ggtitle(paste0("Distributions of PGS"), subtitle ="") +
      xlab(paste0("Standardised PGS-", dis_label[d], "-", pop_label[p])) +
      ylab("Density") +
      theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
      theme(text = element_text(size = 12)) +
      scale_fill_discrete(name = "", labels = c(paste0("Cases (n = ", n_case, ")"), paste0("Controls (n = ", n_control, ")")))
    
    # regression
    model <- glm(paste(dis[d], " ~", pgs, "+ age + age2 + is_female + region_code + 
                       national_pc01 + national_pc02 + national_pc03 + national_pc04 + national_pc05 + national_pc06 + 
                       national_pc07 + national_pc08 + national_pc09 + national_pc10 + national_pc11"), 
                 family = binomial(link = "logit"), data = df)
    summary(model)
    
    assoc_pgs$variable[c] <- pgs
    assoc_pgs$beta[c] <- model$coefficients[2]
    assoc_pgs$SE[c] <- summary(model)$coefficients[2,"Std. Error"]
    assoc_pgs$or[c] <- exp(model$coefficients[pgs])
    assoc_pgs$lci[c] <- exp(confint(model, pgs))[1]
    assoc_pgs$uci[c] <- exp(confint(model, pgs))[2]
    assoc_pgs$p[c] <- summary(model)$coefficients[2,"Pr(>|z|)"]  
    
    # nagelkerke's pseudo r2
    null_model <- glm(paste(dis[d], " ~ age + age2 + is_female + region_code + 
                       national_pc01 + national_pc02 + national_pc03 + national_pc04 + national_pc05 + national_pc06 + 
                       national_pc07 + national_pc08 + national_pc09 + national_pc10 + national_pc11"),
                      family = binomial(link = "logit"), data = df)
    assoc_pgs$r2_nagelkerke[c] <- nagelkerke(model, null = null_model)$Pseudo.R.squared.for.model.vs.null[3]
    
    # R2 on liability scale
    nt=n_total # sample size
    ncase = n_case
    ncont = n_control
    K=prevalence[d] # population prevalence
    P=ncase/nt # proprotion of cases in the sample
    thd = -qnorm(K,0,1)
    zv = dnorm(thd) #z (normal density)
    mv = zv/K #mean liability for case
    lmv = lm(paste(dis[d], "~", pgs),data = df)# linear model
    #R2O : R2 on the observed scale
    R2O = var(lmv$fitted.values)/(ncase/nt*ncont/nt)
    # calculate correction factors
    theta = mv*(P-K)/(1-K)*(mv*(P-K)/(1-K)-thd) 
    cv = K*(1-K)/zv^2*K*(1-K)/(P*(1-P)) 
    # convert to R2 on the liability scale
    assoc_pgs$r2_liability[c] = R2O*cv/(1+R2O*theta*cv)
    
    # plot by quartiles
    df$pgs_quartile <- ntile(df[,pgs], 4)  
    assoc_pgs_quartile <- data.frame(quartile=c("Q1","Q2","Q3","Q4"),es=rep(NA,4),se=rep(NA,4),p=rep(NA,4),or=rep(NA,4),lci=rep(NA,4),hci=rep(NA,4))
    
    for (x in c(1:4)){
      if (x==1){
        assoc_pgs_quartile$es[x] <- 0
        assoc_pgs_quartile$se[x] <- 0
        assoc_pgs_quartile$p[x] <- NA
        
        assoc_pgs_quartile$or[x] <- 1
        assoc_pgs_quartile$lci[x] <- 1
        assoc_pgs_quartile$hci[x] <- 1 
        
        assoc_pgs_quartile$n[x] <- sum(df$pgs_quartile==x & df[dis[d]]==1)
      }
      
      else {
        df_quartile <- df
        df_quartile$pgs_quartile <- as.factor(ifelse(df_quartile$pgs_quartile==x,1, ifelse(df_quartile$pgs_quartile==1,0,NA)))
        summary(df_quartile$pgs_quartile)
        model <- glm(paste(dis[d], " ~ pgs_quartile + age + age2 + is_female + region_code + 
                       national_pc01 + national_pc02 + national_pc03 + national_pc04 + national_pc05 + national_pc06 + 
                       national_pc07 + national_pc08 + national_pc09 + national_pc10 + national_pc11"), 
                     family = binomial(link = "logit"), data = df_quartile)
        
        assoc_pgs_quartile$es[x] <- model$coefficients["pgs_quartile"]
        assoc_pgs_quartile$se[x] <- summary(model)$coefficients["pgs_quartile1","Std. Error"]
        assoc_pgs_quartile$p[x] <- summary(model)$coefficients["pgs_quartile1","Pr(>|z|)"]
        
        assoc_pgs_quartile$or[x] <- exp(model$coefficients["pgs_quartile1"])
        assoc_pgs_quartile$lci[x] <- exp(confint(model, "pgs_quartile1"))[1]
        assoc_pgs_quartile$hci[x] <- exp(confint(model, "pgs_quartile1"))[2]
        
        assoc_pgs_quartile$n[x] <- sum(df$pgs_quartile==x & df[dis[d]]==1)
      }
      print(x)
    }
    
    assoc_pgs_quartile$quartile <- paste0(assoc_pgs_quartile$quartile, " (", assoc_pgs_quartile$n, ")")

    plot_quartile[[p]] <- ggplot(assoc_pgs_quartile, aes(x=quartile, y=or)) +
      geom_point(size=2) +
      geom_errorbar(aes(ymin = lci, ymax = hci),width=0.2) +
      # ylim(0,8) +
      # ggtitle(paste0("Associations with ", dis_label[d]," by quartile"), subtitle = "") +
      xlab(paste("Quartile of", paste0("PGS-", dis_label[d], "-", pop_label[p]), "(cases)")) +
      ylab("Odds ratio") +
      theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
      theme(text = element_text(size = 12)) +
      if (d==1){
        ylim(0.5,6)
      } else {
        ylim(0.9,2.2)
      }
    
    print(p)
    
    c <- c+1
  }
  
  plot[[d]] <- ggarrange(hist[[1]], hist[[2]], plot_quartile[[1]], plot_quartile[[2]], nrow=2, ncol=2, 
                         legend = "top", common.legend = T, labels="auto")
  
}

ggsave("plot_pgs_scz.jpg", plot[[1]], width = 9, height = 9)
ggsave("plot_pgs_mdd.jpg", plot[[2]], width = 9, height = 9)



write.csv(assoc_pgs, "assoc_pgs_logistic.csv", row.names = F)
