install.packages("survival")
install.packages("survminer")
install.packages("maxstat")
install.packages("tidyverse")
install.packages("ggplot2")
library(tidyverse)
library(dplyr)
library(ggplot2)
library(survival)
library(survminer)

set.seed(42)

surv_data <- read.csv("dbta.csv")
#surv_data <- read.csv("csv_512_256/tcga_wt_val_amft.csv")
#surv_data <- read.csv("vgh.csv")
#surv_data <- surv_data %>% mutate(Duration = round(Duration/30.417, digit=0))
#surv_data <- read.csv("dbta.csv")
surv_obj <- Surv(time = surv_data$Duration, event = surv_data$Event)
fit <- survfit(surv_obj ~ COHORT, data = surv_data)

n_high <- sum(surv_data$COHORT == "HighRisk")
n_low <- sum(surv_data$COHORT == "LowRisk")
n_total <- nrow(surv_data)

plot_title <- paste0("TCGA Validation dataset (IDHWT) (n=", n_total, ")")
legend_labels <- c(paste0("High risk (n=", n_high, ")"), paste0("Low risk (n=", n_low, ")"))
# pval.coord = c(65, 0.9),
surv_plot <- ggsurvplot(
  fit,
  data = surv_data,
  pval = TRUE, # Display p-value
  pval.size = 2.8,
  conf.int = TRUE, # Display confidence intervals
  palette = c("#E41A1C", "#377EB8"), # Colors for high and low risk
  surv.median.line = 'hv',
  legend.title = "Risk", # Custom legend title
  title = plot_title,
  legend.labs = legend_labels, # Legend labels
  xlab = "Time (months)", # X-axis label
  ylab = "Survival probability", # Y-axis label
  size=0.7,
  ggtheme = theme_minimal() + theme(plot.title = element_text(hjust = 0.5)) # Apply minimal theme
)
med_surv <- summary(fit)$table[, "median"]
med_text <- paste0("Median Survival HR: ", round(med_surv[[1]], 2), " months")
print(med_text)
surv_plot$plot <- surv_plot$plot +
  annotate("text", x = 5, y = 0.15, label = med_text, size = 2.7, hjust = 0,lineheight = 0.9)
#lineheight = 0.8
med_text <- paste0("Median Survival LR: ", round(med_surv[[2]], 2), " months")
print(med_text)
surv_plot$plot <- surv_plot$plot +annotate("text", x = 5, y = 0.1, label = med_text, size = 2.7, hjust = 0,lineheight = 0.9)

print(surv_plot)


# Save to file with high resolution
ggsave_workaround <- function(g){survminer:::.build_ggsurvplot(x = g,
                                                               surv.plot.height = NULL,
                                                               risk.table.height = NULL,
                                                               ncensor.plot.height = NULL)}

g_to_save <- ggsave_workaround(surv_plot)

ggsave(filename = "tcga_val_wt_80perc.pdf", plot = g_to_save,
       width = 16, height = 12, dpi = 1000, units = "cm")

