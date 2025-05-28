install.packages("survival")
install.packages("survminer")
install.packages("maxstat")

library(survival)
library(survminer)
library(maxstat)
set.seed(42)
surv_data <- read.csv("csv_files/dbta_gigapath.csv")

res.cut <- surv_cutpoint(surv_data, time = "Duration", event = "Event",
                         variables ="haz",minprop=0.2)

summary(res.cut)


surv_data_filtered <- surv_data[surv_data$haz < -1.754976, ]
res.cut <- surv_cutpoint(surv_data_filtered, time = "Duration", event = "Event",
                         variables ="haz",minprop=0.15)

summary(res.cut)
