library(survcomp)

#read the CSV file
data <- read.csv("surv_data_wrisk.csv")

duration <- data$Duration
event <- data$Event

#store c-index values for each hazard score
c_index_values <- numeric(5)

#loop through 5 folds
for (i in 1:5) {
  hazard_score <- data[[paste0("Hazard", i)]]  # Assuming hazard scores are named 'Hazard1', 'Hazard2', ..., 'Hazard5'

  #calculate the c-index
  c_index <- concordance.index(hazard_score, surv.time = duration, surv.event = event)
  c_index_values[i] <- c_index$c.index
}

#compute the average and standard deviation of the c-index values
mean_c_index <- mean(c_index_values)
std_c_index <- sd(c_index_values)

# Output the results
cat("Average c-index: ", mean_c_index , std_c_index)

#get c-index for average hazard score
surv_data <- read.csv("risk_scores.csv")
c_index_result <- concordance.index(surv_data$haz, surv_data$Duration, surv_data$Event)
print(c_index_result$c.index)
