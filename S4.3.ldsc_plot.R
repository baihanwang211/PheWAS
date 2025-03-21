rm(list = ls())

library(ggplot2)
library(dplyr)
library(tidyr)
library(corrplot)
library(grid)
library(gridExtra)
library(ggforce)
library(ggpubr)
library(ggtext)

setwd("")

######################## ldsc
rg <-  read.csv("all.rg", sep="")

# get ci
rg$lci <- rg$rg - 1.96*rg$se
rg$hci <- rg$rg + 1.96*rg$se

# remove _new
rg$p1 <- gsub("_new", "", rg$p1)
rg$p2 <- gsub("_new", "", rg$p2)

# remove suffix
rg$p1 <- gsub("\\.sumstats\\.gz", "", rg$p1)
rg$p2 <- gsub("(_EAS|_EUR)\\.sumstats\\.gz", "", rg$p2)

# duplicate rows between mdd and scz
rg_dup <- rg[rg$p2=="mdd",]
rg_dup$p1 <- c("mdd_EAS", "mdd_EUR")
rg_dup$p2 <- "scz"

rg <- rbind(rg, rg_dup)

rg$p_fdr <- p.adjust(rg$p, method = "fdr")

# split and reorder
rg_scz <- rg[grep("scz",rg$p1),]
rg_scz$p1 <- factor(rg_scz$p1, levels = c("scz_EAS", "scz_EUR"), labels=c("SCZ (EAS)", "SCZ (EUR)"))
# Convert p2 to a factor with specified levels
rg_scz$p2 <- factor(rg_scz$p2, levels = c("mdd", "bmi", "height", "ea", "smkinit", "cigday"), 
                    labels=c("MD", "BMI", "Height", "EA", "SmkInit", "CigDay"))
# Reorder the data frame based on the custom order of p2
rg_scz <- rg_scz %>% arrange(p2)

rg_mdd <- rg[grep("mdd",rg$p1),]
rg_mdd$p1 <- factor(rg_mdd$p1, levels = c("mdd_EAS", "mdd_EUR"), labels=c("MD (EAS)", "MD (EUR)"))
# Convert p2 to a factor with specified levels
rg_mdd$p2 <- factor(rg_mdd$p2, levels = c("scz", "bmi", "smkinit", "cigday", "cataract", "diabetes", "stroke", "cad", "lf", "asthma", "ra", "pud", "gall"),
                    labels=c("SCZ", "BMI", "SmkInit", "CigDay", "Cataract", "Diabetes", "Stroke", "CAD", "Lung function", "Asthma", "RA", "PUD", "Gallstone"))
# Reorder the data frame based on the custom order of p2
rg_mdd <- rg_mdd %>% arrange(p2)

# plot scz

rg_scz$logp <- -log10(rg_scz$p)

rg_scz$sig <- ifelse(rg_scz$p < 0.001, "***", 
                     ifelse(rg_scz$p < 0.01, "**",
                            ifelse(rg_scz$p < 0.05, "*", "") ))

within_scz <- ggplot(rg_scz, aes(x = p2, y = p1)) + 
  geom_point(aes(size = logp, fill = rg), alpha = 1, shape = 21) + 
  scale_size_continuous(limits = c(0, 120), range = c(6,15))  + 
  scale_fill_gradientn(breaks = c(-0.5, 0, 0.5), limits = c(-0.5, 0.5), colors = c("blue", "white", "red")) +
  labs( x= "", y = "", fill = expression(r[g]))  + 
  geom_text(aes(label = sig)) +
  guides(size = "none") +
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 12,  angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(colour = "black",  size = 11), 
        legend.text = element_text(size = 10,  colour ="black"), 
        legend.title = element_text(size = 12), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2))

# plot mdd

rg_mdd$logp <- -log10(rg_mdd$p)

rg_mdd$sig <- ifelse(rg_mdd$p < 0.001, "***", 
                     ifelse(rg_mdd$p < 0.01, "**",
                            ifelse(rg_mdd$p < 0.05, "*", "") ))

within_mdd <- ggplot(rg_mdd, aes(x = p2, y = p1)) + 
  geom_point(aes(size = logp, fill = rg), alpha = 1, shape = 21) + 
  scale_size_continuous(limits = c(0, 120), range = c(6,15))  + 
  scale_fill_gradientn(breaks = c(-0.5, 0, 0.5), limits = c(-0.5, 0.5), colors = c("blue", "white", "red")) +
  labs( x= "", y = "", fill = expression(r[g]))  + 
  geom_text(aes(label = sig)) +
  guides(size = "none") +
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 12,  angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(colour = "black",  size = 11), 
        legend.text = element_text(size = 10,  colour ="black"), 
        legend.title = element_text(size = 12), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2))

