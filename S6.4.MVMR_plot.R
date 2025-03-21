rm(list=ls())

setwd("")

library(ckbplotr)
library(grid)
library(gridExtra)
library(openxlsx)

load("mv_res_exp_scz.RData")
load("mv_res_exp_mdd.RData")
load("mv_res_scz_out.RData")
load("mv_res_mdd_out.RData")
mr_scz <- read.csv("mr_scz.csv")
mr_scz <- mr_scz[mr_scz$pop=="EUR" & mr_scz$method=="Inverse variance weighted",
                 c("exposure", "outcome", "b", "se", "pval", "nsnp", "F")]

mr_mdd <- read.csv("mr_mdd.csv")
mr_mdd <- mr_mdd[mr_mdd$pop=="EUR" & mr_mdd$method=="Inverse variance weighted",
                 c("exposure", "outcome", "b", "se", "pval", "nsnp", "F")]

mr_res <- rbind(mr_scz, mr_mdd)
mr_res$Q <- NA
mr_res$Q_pval <- NA
mr_res$b_robust <- NA
mr_res$covar <- "univariable"

pheno <- c("bmi", "smkinit", "cigday")

# pheno to psych

exp_psych_list <- list()

for (x in 1:length(pheno)){
  
  exp_psych <- data.frame()
  
  for (i in 1:length(mv_res_exp_scz[[1]])){
    df <- mv_res_exp_scz[[x]][[i]]
    df <- df[c("exposure", "outcome", "Estimate", "Std..Error", "Pr...t..", "nsnp", "F", "Qstat", "Qpval", "Estimates_robust")]
    names(df) <- c("exposure", "outcome", "b", "se", "pval", "nsnp", "F", "Q", "Q_pval", "b_robust")
    df$covar <- setdiff(df$exposure, pheno[x])
    exp_psych <- rbind(exp_psych, df[df$exposure==pheno[x],])
    
    df <- mv_res_exp_mdd[[x]][[i]]
    df <- df[c("exposure", "outcome", "Estimate", "Std..Error", "Pr...t..", "nsnp", "F", "Qstat", "Qpval", "Estimates_robust")]
    names(df) <- c("exposure", "outcome", "b", "se", "pval", "nsnp", "F", "Q", "Q_pval", "b_robust")
    df$covar <- setdiff(df$exposure, pheno[x])
    exp_psych <- rbind(exp_psych, df[df$exposure==pheno[x],])
  }
  
  exp_psych_list[[x]] <- rbind(mr_res[mr_res$exposure==pheno[x],], exp_psych)
  
  exp_psych_list[[x]]$covar <- factor(exp_psych_list[[x]]$covar,
                                      levels = c("univariable", "bmi", "smkinit", "cigday", "canud", "ea", "income"),
                                      labels = c("Univariable", "+ BMI", "+ SmkInit", "+ CigDay", "+ CanUD", "+ EA", "+ Income"))
  
}


# psych to pheno

psych_out_list <- list()

for (x in 1:length(pheno)){
  
  psych_out <- data.frame()
  
  for (i in 1:length(mv_res_scz_out[[1]])){
    df <- mv_res_scz_out[[x]][[i]]
    df <- df[c("exposure", "outcome", "Estimate", "Std..Error", "Pr...t..", "nsnp", "F", "Qstat", "Qpval", "Estimates_robust")]
    names(df) <- c("exposure", "outcome", "b", "se", "pval", "nsnp", "F", "Q", "Q_pval", "b_robust")
    df$covar <- setdiff(df$exposure, "scz")
    psych_out <- rbind(psych_out, df[df$exposure=="scz",])
    
    df <- mv_res_mdd_out[[x]][[i]]
    df <- df[c("exposure", "outcome", "Estimate", "Std..Error", "Pr...t..", "nsnp", "F", "Qstat", "Qpval", "Estimates_robust")]
    names(df) <- c("exposure", "outcome", "b", "se", "pval", "nsnp", "F", "Q", "Q_pval", "b_robust")
    df$covar <- setdiff(df$exposure, "mdd")
    psych_out <- rbind(psych_out, df[df$exposure=="mdd",])
  }
  
  psych_out_list[[x]] <- rbind(mr_res[mr_res$outcome==pheno[x],], psych_out)
  psych_out_list[[x]]$covar <- factor(psych_out_list[[x]]$covar,
                                      levels = c("univariable", "bmi", "smkinit", "cigday", "canud", "ea", "income"),
                                      labels = c("Univariable", "+ BMI", "+ SmkInit", "+ CigDay", "+ CanUD", "+ EA", "+ Income"))
  
}

### plot together

forest_list <- list()

pheno_label <- c("BMI", "SmkInit", "CigDay")

exp_psych_list_1 <- exp_psych_list
psych_out_list_1 <- psych_out_list

for (x in 1:length(pheno)){
  
  if (x !=1){
    exp_psych_list_1[[x]] <- exp_psych_list_1[[x]][exp_psych_list_1[[x]]$covar!="+ SmkInit" & exp_psych_list_1[[x]]$covar!="+ CigDay",]
    psych_out_list_1[[x]] <- psych_out_list_1[[x]][psych_out_list_1[[x]]$covar!="+ SmkInit" & psych_out_list_1[[x]]$covar!="+ CigDay",]
  }
  
  exp_psych_list_1[[x]]$pval <- formatC(exp_psych_list_1[[x]]$pval, digits = 3, format = "f")
  exp_psych_list_1[[x]]$pval[exp_psych_list_1[[x]]$pval < 0.001] <- "< 0.001"
  
  psych_out_list_1[[x]]$pval <- formatC(psych_out_list_1[[x]]$pval, digits = 3, format = "f")
  psych_out_list_1[[x]]$pval[psych_out_list_1[[x]]$pval < 0.001] <- "< 0.001"
  
  forest_1 <- forest_plot(exp_psych_list_1[[x]][exp_psych_list_1[[x]]$outcome=="scz",],
                          row.labels.heading	= "Exposure",
                          col.key = "covar",
                          col.estimate = "b",
                          col.stderr = "se",
                          estcolumn = FALSE,
                          col.right = "pval",
                          col.right.heading = "p",
                          xlab = "Beta (95% CI)",
                          xlim = c(-0.5,1.2),
                          exponentiate=FALSE,
                          nullval=0, pointsiz=2,
                          panel.headings = paste("From", pheno_label[x], "to SCZ"),
                          plot.margin = margin(5, 5, 5, 5, "mm"))
  
  
  forest_2 <- forest_plot(exp_psych_list_1[[x]][exp_psych_list_1[[x]]$outcome=="mdd",],
                          col.key = "covar",
                          row.labels = data.frame(covar = unique(exp_psych_list_1[[x]]$covar), label = rep(" ", length(unique(exp_psych_list_1[[x]]$covar)))),
                          col.estimate = "b",
                          col.stderr = "se",
                          col.left = NULL,
                          estcolumn = FALSE,
                          col.right = "pval",
                          col.right.heading = "p",
                          xlab = "Beta (95% CI)",
                          xlim = c(-0.2,0.6),
                          exponentiate=FALSE,
                          nullval=0, pointsiz=2,
                          panel.headings = paste("From", pheno_label[x], "to MD"),
                          plot.margin = margin(5, 5, 5, 5, "mm"))
  
  
  forest_3 <- forest_plot(psych_out_list_1[[x]][psych_out_list_1[[x]]$exposure=="scz",],
                          col.key = "covar",
                          row.labels = data.frame(covar = unique(psych_out_list_1[[x]]$covar), label = rep(" ", length(unique(psych_out_list_1[[x]]$covar)))),
                          col.estimate = "b",
                          col.stderr = "se",
                          estcolumn = FALSE,
                          col.right = "pval",
                          col.right.heading = "p",
                          xlab = "Beta (95% CI)",
                          xlim = c(-0.05,0.06),
                          exponentiate=FALSE,
                          nullval=0, pointsiz=2,
                          panel.headings = paste("From SCZ to", pheno_label[x]),
                          plot.margin = margin(5, 5, 5, 5, "mm"))
  
  
  
  forest_4 <- forest_plot(psych_out_list_1[[x]][psych_out_list_1[[x]]$exposure=="mdd",],
                          col.key = "covar",
                          row.labels = data.frame(covar = unique(psych_out_list_1[[x]]$covar), label = rep(" ", length(unique(psych_out_list_1[[x]]$covar)))),
                          col.estimate = "b",
                          col.stderr = "se",
                          col.left = NULL,
                          estcolumn = FALSE,
                          col.right = "pval",
                          col.right.heading = "p",
                          xlab = "Beta (95% CI)",
                          xlim = c(-0.1,0.23),
                          exponentiate=FALSE,
                          nullval=0, pointsiz=2,
                          panel.headings = paste("From MD to", pheno_label[x]),
                          plot.margin = margin(5, 5, 5, 5, "mm"))
  
  
  forest_list[[x]] <- grid.arrange(forest_1$plot, forest_2$plot, forest_3$plot, forest_4$plot, ncol = 4, widths = c(1.2, 1, 1, 1))
  
}

grid.newpage()

forest_combined <- grid.arrange(arrangeGrob(forest_list[[1]], left = textGrob(expression(bold("a")), x = unit(1, "npc"), y = unit(.95, "npc"))),
                                arrangeGrob(forest_list[[2]], left = textGrob(expression(bold("b")), x = unit(1, "npc"), y = unit(.95, "npc"))),
                                arrangeGrob(forest_list[[3]], left = textGrob(expression(bold("c")), x = unit(1, "npc"), y = unit(.95, "npc"))),
                                nrow = 3)

ggsave("J:/psych/results/mdd_new/mvmr_scale.png", forest_combined, width = 12, height = 8, bg = "white")
