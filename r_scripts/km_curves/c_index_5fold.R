BiocManager::install("survcomp") 
library(survcomp)

# Step 2: Read the CSV file
data <- read.csv("csv_files/dbta_gigapath.csv")

duration <- data$Duration   # Assumingqthe column name is 'duration'
event <- data$Event         # Assuming the column name is 'event'

# Step 4: Create a vector to store c-index values for each hazard score
c_index_values <- numeric(5)

# Step 5: Calculate c-index for each hazard score
for (i in 1:5) {
  hazard_score <- data[[paste0("Hazard", i)]]  # Assuming hazard scores are named 'hazard_score1', 'hazard_score2', ..., 'hazard_score5'
  
  # Calculate the c-index
  c_index <- concordance.index(hazard_score, surv.time = duration, surv.event = event)
  
  # Store the c-index value
  c_index_values[i] <- c_index$c.index
}
c_index_values
# Step 6: Compute the average and standard deviation of the c-index values
mean_c_index <- mean(c_index_values)
std_c_index <- sd(c_index_values)

# Output the results
cat("Average c-index: ", mean_c_index , std_c_index)