rm(list=ls())

setwd("")

library(ckbplotr)
library(grid)
library(gridExtra)
library(openxlsx)

###### scz
mr_res <- read.csv("mr_scz.csv")
mr_res <- mr_res[mr_res$method=="Inverse variance weighted",]
mr_res <- mr_res[,-which(names(mr_res) %in% c("id.exposure", "id.outcome", "method"))]
mr_res$subgroup <- paste(mr_res$exposure, mr_res$outcome, mr_res$pop, sep="_")

mr_res$pval <- formatC(mr_res$pval, digits = 3, format = "f")
mr_res$pval[mr_res$pval < 0.001] <- "< 0.001"

row_labels <- data.frame(
  subgroup = c("bmi_scz_EAS", "bmi_scz_EUR",
               "smkinit_scz_EAS", "smkinit_scz_EUR",
               "cigday_scz_EAS", "cigday_scz_EUR"),
  group    = c("BMI","BMI",
               "SmkInit", "SmkInit", 
               "CigDay", "CigDay"),
  label    = c("EAS", "EUR",
               "EAS", "EUR",
               "EAS", "EUR")
)

forest_mr_exp_scz <- ckbplotr::forest_plot(mr_res[mr_res$outcome=="scz",],
                                           col.key = "subgroup",
                                           row.labels.heading	= "Exposure",
                                           row.labels = row_labels,
                                           col.estimate = "b",
                                           col.stderr = "se",
                                           col.right = "pval",
                                           col.right.heading = c("Beta (95% CI)","p"),
                                           col.right.hjust=0.5,
                                           xlab = "Beta (95% CI)",
                                           xlim = c(-0.4,1.0),
                                           panel.headings = "From exposure to SCZ",
                                           exponentiate=FALSE,
                                           nullval=0, pointsiz=2)
forest_mr_exp_scz


row_labels <- data.frame(
  subgroup = c("scz_bmi_EAS", "scz_bmi_EUR",
               "scz_smkinit_EAS", "scz_smkinit_EUR",
               "scz_cigday_EAS", "scz_cigday_EUR"),
  group    = c("BMI","BMI",
               "SmkInit", "SmkInit", 
               "CigDay", "CigDay"),
  label    = c("EAS", "EUR",
               "EAS", "EUR",
               "EAS", "EUR")
)

forest_mr_scz_out <- ckbplotr::forest_plot(mr_res[mr_res$exposure=="scz",],
                                           col.key = "subgroup",
                                           row.labels.heading	= "Outcome",
                                           row.labels = row_labels,
                                           col.estimate = "b",
                                           col.stderr = "se",
                                           col.right = "pval",
                                           col.right.heading = c("Beta (95% CI)","p"),
                                           col.right.hjust=0.5,
                                           xlab = "Beta (95% CI)",
                                           xlim = c(-0.05,0.23),
                                           panel.headings = "From SCZ to outcome",
                                           exponentiate=FALSE,
                                           nullval=0, pointsiz=2)
forest_mr_scz_out

########### mdd

mr_res <- read.csv("mr_mdd.csv")
mr_res <- mr_res[mr_res$method=="Inverse variance weighted",]
mr_res <- mr_res[,-which(names(mr_res) %in% c("id.exposure", "id.outcome", "method"))]
mr_res$subgroup <- paste(mr_res$exposure, mr_res$outcome, mr_res$pop, sep="_")

mr_res$pval <- formatC(mr_res$pval, digits = 3, format = "f")
mr_res$pval[mr_res$pval < 0.001] <- "< 0.001"

row_labels <- data.frame(
  subgroup = c("bmi_mdd_EAS", "bmi_mdd_EUR",
               "smkinit_mdd_EAS", "smkinit_mdd_EUR",
               "cigday_mdd_EAS", "cigday_mdd_EUR"),
  group    = c("BMI","BMI",
               "SmkInit", "SmkInit", 
               "CigDay", "CigDay"),
  label    = c("EAS", "EUR",
               "EAS", "EUR",
               "EAS", "EUR")
)

forest_mr_exp_mdd <- ckbplotr::forest_plot(mr_res[mr_res$outcome=="mdd",],
                                           col.key = "subgroup",
                                           row.labels.heading	= "Exposure",
                                           row.labels = row_labels,
                                           col.estimate = "b",
                                           col.stderr = "se",
                                           col.right = "pval",
                                           col.right.heading = c("Beta (95% CI)","p"),
                                           col.right.hjust=0.5,
                                           xlab = "Beta (95% CI)",
                                           xlim = c(-0.4,1.0),
                                           panel.headings = "From exposure to MD",
                                           exponentiate=FALSE,
                                           nullval=0, pointsiz=2)
forest_mr_exp_mdd


row_labels <- data.frame(
  subgroup = c("mdd_bmi_EAS", "mdd_bmi_EUR",
               "mdd_smkinit_EAS", "mdd_smkinit_EUR",
               "mdd_cigday_EAS", "mdd_cigday_EUR"),
  group    = c("BMI","BMI",
               "SmkInit", "SmkInit", 
               "CigDay", "CigDay"),
  label    = c("EAS", "EUR",
               "EAS", "EUR",
               "EAS", "EUR")
)

forest_mr_mdd_out <- ckbplotr::forest_plot(mr_res[mr_res$exposure=="mdd",],
                                           col.key = "subgroup",
                                           row.labels.heading	= "Outcome",
                                           row.labels = row_labels,
                                           col.estimate = "b",
                                           col.stderr = "se",
                                           col.right = "pval",
                                           col.right.heading = c("Beta (95% CI)","p"),
                                           col.right.hjust=0.5,
                                           xlab = "Beta (95% CI)",
                                           xlim = c(-0.05,0.23),
                                           panel.headings = "From MD to outcome",
                                           exponentiate=FALSE,
                                           nullval=0, pointsiz=2)
forest_mr_mdd_out

grid.newpage()
mr <- grid.arrange(arrangeGrob(forest_mr_exp_scz$plot, left = textGrob(expression(bold("a")), x = unit(1, "npc"), y = unit(.95, "npc"))),
                   arrangeGrob(forest_mr_scz_out$plot, left = textGrob(expression(bold("b")), x = unit(1, "npc"), y = unit(.95, "npc"))),
                   arrangeGrob(forest_mr_exp_mdd$plot, left = textGrob(expression(bold("c")), x = unit(1, "npc"), y = unit(.95, "npc"))),
                   arrangeGrob(forest_mr_mdd_out$plot, left = textGrob(expression(bold("d")), x = unit(1, "npc"), y = unit(.95, "npc"))),
                   ncol = 2, nrow = 2)

ggsave("mr_scale.png",  mr, width = 10, height = 6, bg = "white")

