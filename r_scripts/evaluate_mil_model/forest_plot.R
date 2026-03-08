library(survival)
library(dplyr)
library(survminer)
library(maxstat)
library(caret)
library(glmnet)
library(survcomp)
set.seed(42)

tcga_val <- read.csv("surv_clinical_vars.csv")

colnames(tcga_val)[colnames(tcga_val) == "haz"]                  <- "Hazard Score"
colnames(tcga_val)[colnames(tcga_val) == "c1p19q"]               <- "1p/19q Codeletion"
colnames(tcga_val)[colnames(tcga_val) == "age"]                  <- "Age (>= 40 years)"
colnames(tcga_val)[colnames(tcga_val) == "MGMT.promoter.status"] <- "MGMT Promoter Status"
colnames(tcga_val)[colnames(tcga_val) == "Grade"]                <- "Grade"
colnames(tcga_val)[colnames(tcga_val) == "CDKN2A_B_CNV"]        <- "CDKN2A/B CNV"
colnames(tcga_val)[colnames(tcga_val) == "RadiationT"]           <- "Radiation Therapy"
colnames(tcga_val)[colnames(tcga_val) == "Sex"]                  <- "Sex"

features_plot <- c("Hazard Score", "1p/19q Codeletion", "Age (>= 40 years)", "MGMT Promoter Status",
                   "Grade", "CDKN2A/B CNV", "Radiation Therapy", "Sex")
cox_formula <- as.formula(
  paste("Surv(time, status) ~", paste(paste0("`", features_plot, "`"), collapse = " + "))
)

multicox <- coxph(cox_formula, data = tcga_val)
summary(multicox)
multisum <- summary(multicox)
ggforest(multicox,main = "Forest Plot - OS", fontsize = 1.1)
