library(survcomp)

surv_data <- read.csv("csv_files/csv_ablation/5_cut_DBTA.csv")
c_index_result <- concordance.index(surv_data$haz, surv_data$Duration, surv_data$Event)
print(c_index_result$c.index)


