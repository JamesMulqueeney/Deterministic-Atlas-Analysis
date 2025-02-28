

library(gtools)  # for combinations function
library(dplyr)

# Filter to keep only Orders with more than 10 specimens
order_counts <- pcscores %>%
  group_by(Order) %>%
  summarise(specimen_count = n()) %>%
  filter(specimen_count > 10)  # Keep only orders with more than 10 specimens

# Get the filtered list of Orders
filtered_orders <- order_counts$Order

# Create an empty data frame to store results
confusion_matrix_results <- data.frame()

# Loop through all combinations of Orders (based on included orders)
for (i in 1:length(filtered_orders)) {
  combinations <- combinations(length(filtered_orders), i, filtered_orders)  # Generate all combinations of 'i' Orders
  
  for (combo in 1:nrow(combinations)) {
    included_orders <- combinations[combo, ]
    
    # Filter data by including only the selected orders
    data_filtered <- merged_data_filtered4 %>%
      filter(Order %in% included_orders)
    
    # Fit linear model and extract R-squared and p-value
    lm_model <- lm(PC_distance ~ Comp_distance, data = data_filtered)
    lm_summary <- summary(lm_model)
    
    r_squared <- lm_summary$r.squared
    p_value <- lm_summary$coefficients["Comp_distance", "Pr(>|t|)"]
    
    # Store results in the confusion matrix
    confusion_matrix_results <- rbind(
      confusion_matrix_results,
      data.frame(
        included_orders = paste(included_orders, collapse = ", "),  # Name of included orders
        r_squared = r_squared,
        p_value = p_value
      )
    )
  }
}

# View the results
print(confusion_matrix_results)
