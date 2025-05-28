library(survival)
library(survminer)
library(readxl)
library(tidyverse)
library(survival)
library(survminer)
library(forestplot)

surv_data <- read.csv('tcga_val_multicoxph_cleaned.csv')

multicox <- coxph(Surv(time = time, event = status) ~ ., data = surv_data)
multisum <- summary(multicox)

Factor <- colnames(surv_data)[3:10]
HR <- multisum$coefficients[, 2]
L95CI <- multisum$conf.int[, 3]
H95CI <- multisum$conf.int[, 4]
pvalue <- multisum$coefficients[, 5]


multiresult <- data.frame(
  Factor = Factor,
  HR = round(HR, 6),
  L95CI = round(L95CI, 6),
  H95CI = round(H95CI, 6),
  pvalue = round(pvalue, 6)
)

multiresult <- multiresult[multiresult$pvalue < 1,]

multiresult$Factor <- c("Risk Score", "1p19q Status", "Age (>40ys)", "MGMTp Methyl","Histologic Grade", "CDKN2A/B Status","Radiation Therapy","Gender")

write.csv(multiresult, file = 'multi_cox_allfactor.csv', row.names = FALSE)

data <- read.csv('multi_cox_allfactor.csv')


data$log10HR <- log10(data$HR)
data$log10L95CI <- log10(data$L95CI)
data$log10H95CI <- log10(data$H95CI)



data$original_log10L95CI <- data$log10L95CI
data$original_log10H95CI <- data$log10H95CI


data$log10L95CI <- with(data, ifelse(
  (original_log10L95CI > 0 & original_log10H95CI > 0) | (original_log10L95CI < 0 & original_log10H95CI < 0), 
  ifelse(abs(original_log10L95CI) > abs(original_log10H95CI), 
         formatC(original_log10H95CI, format = "f", digits = 4), 
         formatC(original_log10L95CI, format = "f", digits = 4)),
  formatC(original_log10H95CI, format = "f", digits = 4)
))

data$log10H95CI <- with(data, ifelse(
  (original_log10L95CI > 0 & original_log10H95CI > 0) | (original_log10L95CI < 0 & original_log10H95CI < 0), 
  ifelse(abs(original_log10L95CI) > abs(original_log10H95CI), 
         formatC(original_log10L95CI, format = "f", digits = 4), 
         formatC(original_log10H95CI, format = "f", digits = 4)),
  formatC(original_log10L95CI, format = "f", digits = 4)
))


data$log10HR <- formatC(data$log10HR, format = "f", digits = 4)


data$`log10HR(95% CI)` <- with(data, 
                               ifelse(
                                 log10L95CI < 0 & log10H95CI < 0,
                                 paste0(
                                   formatC(log10HR, format = "f", digits = 4),
                                   " (", 
                                   formatC(log10H95CI, format = "f", digits = 4),
                                   ", ", 
                                   formatC(log10L95CI, format = "f", digits = 4),
                                   ")"
                                 ),
                                 paste0(
                                   formatC(log10HR, format = "f", digits = 4),
                                   " (", 
                                   ifelse(log10L95CI < log10H95CI, 
                                          formatC(log10L95CI, format = "f", digits = 4),
                                          formatC(log10H95CI, format = "f", digits = 4)),
                                   ", ", 
                                   ifelse(log10L95CI < log10H95CI, 
                                          formatC(log10H95CI, format = "f", digits = 4),
                                          formatC(log10L95CI, format = "f", digits = 4)),
                                   ")"
                                 )
                               )
)


print(data$`log10HR(95% CI)`)


data$pvalue <- ifelse(data$pvalue < 0.0001, 
                      '<0.0001', 
                      sprintf("%.4f", data$pvalue))

column_names <- colnames(data)


column_names_row <- matrix(column_names, nrow = 1)
colnames(column_names_row) <- NULL

data_matrix <- as.matrix(data)


data_with_colnames <- rbind(column_names_row, data_matrix)

result1 <- data_with_colnames


fig5 <- forestplot(result1[,c(1,11,5)], 
                   mean=as.numeric(result1[,6]),   
                   lower=as.numeric(result1[,9]),  
                   upper=as.numeric(result1[,10]), 
                   zero=0,           
                   boxsize=0.3,      
                   graph.pos="right",
                   hrzl_lines=list("1" = gpar(lty=1,lwd=2),
                                   "2" = gpar(lty=2),
                                   "10"= gpar(lwd=2,lty=1,columns=c(1:4)) ),
                   graphwidth=unit(0.25, "npc"),
                   xlab="log10 Hazard Ratio",
                   xticks=c(-1.0, 0, 1.0),
                   is.summary=c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
                   txt_gp=fpTxtGp(
                     label = gpar(cex = 1),
                     ticks = gpar(cex = 1), 
                     xlab = gpar(cex = 1.2), 
                     title = gpar(cex = 2)
                   ),
                   lwd.zero=2,
                   lwd.ci=2,
                   lwd.xaxis=3, 
                   lty.ci=1.5,
                   ci.vertices=TRUE,
                   ci.vertices.height=0.25, 
                   clip=c(-0.7, 1.0),
                   lineheight=unit(10, 'mm'), 
                   line.margin=unit(6, 'mm'),
                   colgap=unit(14, 'mm'),
                   fn.ci_norm="fpDrawDiamondCI", 
                   col=fpColors(box='#ca7465', 
                                lines='black', 
                                zero="#7ac6ce")) 

print(fig5)

pdf('1_multicox_forestplot.pdf', width = 9, height = 4.5) 
print(fig5)
dev.off()