library(tidyverse)
library(dplyr)
library(ggplot2)
library(survival)
library(survminer)
set.seed(42)

#read in csv
surv_data <- read.csv("surv_data_wrisk.csv")
surv_obj <- Surv(time = surv_data$Duration.Months, event = surv_data$Event)
fit <- survfit(surv_obj ~ COHORT, data = surv_data)

n_high <- sum(surv_data$COHORT == "HighRisk")
n_low <- sum(surv_data$COHORT == "LowRisk")
n_total <- nrow(surv_data)

plot_title <- paste0("TEST dataset (n=", n_total, ")")
legend_labels <- c(paste0("High risk (n=", n_high, ")"), paste0("Low risk (n=", n_low, ")"))

surv_plot <- ggsurvplot(
  fit,
  data = surv_data,
  pval = TRUE, # Display p-value
  pval.size = 2.8,
  pval.coord = c(4, 0.2), #use to change location of p-value
  conf.int = TRUE, # Display confidence intervals
  palette = c("#E41A1C", "#377EB8"), # Colors for high and low risk
  surv.median.line = 'hv',
  legend.title = "Risk", # Custom legend title
  title = plot_title,
  legend.labs = legend_labels, # Legend labels
  xlab = "Time (months)", # X-axis label
  ylab = "Survival probability", # Y-axis label
  size=0.7,
  ggtheme = theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
)
med_surv <- summary(fit)$table[, "median"]

#use to place the median survival of HR and LR
med_text <- paste0("Median Survival HR: ", round(med_surv[[1]], 2), " months")
print(med_text)
surv_plot$plot <- surv_plot$plot +
  annotate("text", x = 4, y = 0.15, label = med_text, size = 2.7, hjust = 0,lineheight = 0.9)
med_text <- paste0("Median Survival LR: ", round(med_surv[[2]], 2), " months")
print(med_text)
surv_plot$plot <- surv_plot$plot +annotate("text", x = 4, y = 0.1, label = med_text, size = 2.7, hjust = 0,lineheight = 0.9)

print(surv_plot)


# Save to file with high resolution
ggsave_workaround <- function(g){survminer:::.build_ggsurvplot(x = g,
                                                               surv.plot.height = NULL,
                                                               risk.table.height = NULL,
                                                               ncensor.plot.height = NULL)}

g_to_save <- ggsave_workaround(surv_plot)

ggsave(filename = "km_curve.pdf", plot = g_to_save,
       width = 16, height = 12, dpi = 1000, units = "cm")
