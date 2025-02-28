# Paper - Assessing the application of landmark-free morphometrics to macroevolutionary analyses

# Author: James M. Mulqueeney 

# Date Last Modified: 28/02/2025

# Within Order Euclidean Distance Statistics 

################################################################################

# Load required libraries
library(dplyr)

################################################################################

# Euclidean Distance Order Correlations 

################################################################################

# Define input directory (change accordingly)
input_dir <- "C:/Users/jmm1e21/OneDrive - University of Southampton/Documents/Main PhD Work/Chapter 3 - Comparison of Automated Methods/Final/BMC Ecology and Evolution/Data & Code/Data/Data"

# Read the CSV file into a data frame
data <- read.csv(file.path(input_dir, "Data_A10-Order_Euclidean_Distance_Correlations.csv"))

# Show the Data 
print(data)

# Ensure column names are correct
colnames(data)

# Calculate mean and standard deviation for Aligned_Only and Poisson columns by Order
order_stats <- data %>%
  group_by(Order) %>%
  summarise(
    Aligned_Mean = mean(c(Aligned_Only_40, Aligned_Only_20, Aligned_Only_10), na.rm = TRUE),
    Aligned_SD = sd(c(Aligned_Only_40, Aligned_Only_20, Aligned_Only_10), na.rm = TRUE),
    Poisson_Mean = mean(c(Poisson_40, Poisson_20, Poisson_10), na.rm = TRUE),
    Poisson_SD = sd(c(Poisson_40, Poisson_20, Poisson_10), na.rm = TRUE),
    .groups = "drop"
  )

# View the results
print(order_stats)

#########################################################################################

# Compare affects of Mixed Mesh Modalities on Each Group 

#########################################################################################

# Perform paired t-test for each Order (Poisson compared to Aligned-Only)
ttest_results <- data %>%
  group_by(Order) %>%
  summarise(
    t_statistic = t.test(c(Poisson_40, Poisson_20, Poisson_10),  # Poisson first
                         c(Aligned_Only_40, Aligned_Only_20, Aligned_Only_10),
                         paired = TRUE)$statistic,
    p_value = t.test(c(Poisson_40, Poisson_20, Poisson_10),
                     c(Aligned_Only_40, Aligned_Only_20, Aligned_Only_10),
                     paired = TRUE)$p.value,
    mean_diff = mean(c(Poisson_40, Poisson_20, Poisson_10)) -  # Poisson - Aligned
      mean(c(Aligned_Only_40, Aligned_Only_20, Aligned_Only_10)),
    .groups = "drop"
  )

# View the reversed t-test results
print(ttest_results)





